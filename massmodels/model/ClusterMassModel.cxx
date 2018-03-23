// Mixing model for determining the gravitating mass that confines
// a spherically symmetric, hydrostatic atmosphere.
//
// ClusterMassModel and Hydrostatic::Atmosphere are base classes
// from which models for the gravitational potential are derived.


#include <XSModel/GlobalContainer/DataContainer.h>
using namespace XSContainer;
#include <XSModel/Data/SpectralData.h>
#include <XSUtil/Numerics/LinearInterp.h>
#include <cmath>
#include <limits>
#include <gsl/gsl_sf_gamma.h>
#include "XSContainer.h"
#include "XSstreams.h"
#include "Potential.h"


const int ClusterMassModel::chatterInfo = 20;
const int ClusterMassModel::chatterDebug = 30;

Real ClusterMassModel::Sector::FLAG = -999.0;

ClusterMassModel::Sector::Sector (Real r_inner, Real r_outer, Real delta_angle)
  : rin (r_inner), rout (r_outer), angle (delta_angle) {
}

int ClusterMassModel::Sector::operator== (const ClusterMassModel::Sector &right)
  const {
  // For comparing spectral keywords, require radial bounds to match
  return (rin == right.rin) && (rout == right.rout);
}

int ClusterMassModel::Sector::operator!= (const ClusterMassModel::Sector &right)
  const {
  return !operator==(right);
}

ClusterMassModel::ClusterMassModel (const string& name)
  : MixUtility (name) {
}

ClusterMassModel::~ClusterMassModel () {
}

ClusterMassModel::Sector ClusterMassModel::testSector 
(const std::map<size_t, Real>& xflt, size_t spectrumNumber) {
  // At least one XFLT must be specified
  if (xflt.size() >= 1) {
    return Sector (0.0, xflt.find(1)->second, 360.0);
  } else {
    std::ostringstream diag;
    diag << "\n***  Spectrum number " << spectrumNumber 
	 << " has no annulus definition\n*** model removed";
    throw IncompatibleData (diag.str());                
  }
}

void ClusterMassModel::verifyData ()  {
  // Checks for:
  // data present
  // consistent definitions of annuli within data groups
   
  // Adapted to Xspec 12.4.0, which allows selective inclusion of spectra
  // in a model.

  if (datasets->dataArray ().size () == 0) {
    throw NoDataPresent ("");
  }

  // Mixing model assumed only to apply to source number 1
  const size_t SOURCENUM = 1;
  
  typedef std::map<size_t, size_t> dgMapType;
  typedef std::map<size_t, dgMapType> sourceDgMapType;
  sourceDgMapType::const_iterator itSourceDgMap 
    = datasets->sourceToDgs ().find (SOURCENUM);
  if (itSourceDgMap == datasets->sourceToDgs ().end ()) {
    throw YellowAlert ("No data sets are assigned to mixing model.\n");
  }
  const dgMapType& dgMap = itSourceDgMap->second;

  // Prepare to check consistency of sector definitions for data groups
  std::map<size_t, Sector> sectorTest;
  dgMapType::const_iterator itDgMap = dgMap.begin ();
  dgMapType::const_iterator itDgMapEnd = dgMap.end ();
  while (itDgMap != itDgMapEnd) {
    sectorTest [itDgMap->first] = Sector ();
    ++itDgMap;
  }

  // Check consistency of sector definitions
  size_t nTotalSpec = datasets->numberOfSpectra ();
  for (size_t iSpec = 1; iSpec <= nTotalSpec; ++iSpec) {
    const SpectralData *spec = datasets->lookup (iSpec);
    size_t groupNum = spec->parent ()->dataGroup ();
    // Spectrum active for the mixing model source?
    if (spec->detector (SOURCENUM - 1)) {
      if (sectorTest [groupNum].rout == Sector::FLAG) {
	sectorTest [groupNum] = testSector (convertToOldXfltMap(spec->xflt ()), iSpec);
      } else {
	Sector currentSector (testSector (convertToOldXfltMap(spec->xflt ()), iSpec));
	if (currentSector != sectorTest[groupNum]) {
	  std::ostringstream msg;
	  msg << "\n***  Spectrum: " << iSpec << " in data group " << groupNum
	      << "\n***  contains spheroid definition inconsistent "
	      << "with the group (check data entry order)";
	  throw DataOrderingError (msg.str ());
	}
      }
    }
  }
}

// Make a Sector for each spectrum, filling in properties that
// can be determined from XFLT keywords for this spectrum.  Does
// not include rin.
void ClusterMassModel::createSectorDescription (size_t observation, 
						size_t spectrumNumber,
						size_t iShell) {
  SpectralData* spectrum (datasets->lookup (spectrumNumber));      
  const std::map<size_t, Real> filterKeys = convertToOldXfltMap(spectrum->xflt ());  
  size_t N (filterKeys.size ());  
  
  if (N < 1) {
    string diag(" 1 keyword (outer radius) required");
    throw IncompatibleData (diag);
  }

  Sector& currShell = m_sector [iShell][observation];  
  currShell.rout = filterKeys.find (1)->second;
  // Fill rbounds, unsorted
  rbounds [iShell + 1] = currShell.rout;

  if (N == 1) {
    currShell.angle = 360.0;  
  } else {
    // "Ellipse" must be circular
    Real b = filterKeys.find (2)->second;
    if (b != currShell.rout) {
      string diag(" mass models requires circular annuli");
      throw IncompatibleData (diag);
    }

    if (N <= 4) {
      currShell.angle = 360.0;
    } else {
      currShell.angle = 0.0;
      // Allow multiple angle ranges (compatible with projct)
      for (size_t j = 4; j < N; j += 2) {
	Real deg = filterKeys.find (j + 1)->second 
	  - filterKeys.find (j)->second;
	while (deg < 0.0) {
	  deg += 360.0;
	}
	currShell.angle += deg;
      }
    }
  }
}

void ClusterMassModel::makeBetaEM (const std::vector<Real>& params, 
				   bool force) {

  // Computes emission measure of gas outside the inner boundary of
  // the outermost shell in each annulus, assuming unit density at the
  // inner edge of that shell.  
  // Results are divided by 4 pi, to be consistent with emission measure
  // computation in Atmosphere.

  // Allows acore == 0.
  // Must have beta > 1/6 and not very close, but unequal to 0.5.
  // This is potentially a trap if the model is used to fit beta.

  // Try to avoid calculation.  A change to rinner or turning switch on
  // triggers initialize(), which sets rbounds[] and useBetaModel, then
  // calls this with force true, setting acore, beta and betaModelEM. 
  // If force is false, only need to recompute betaModelEM if acore or
  // beta have changed.
  if (!useBetaModel || !force && acore == params[1] && beta == params[2])
    return;
  acore = params[1];
  beta = params[2];

  if (beta <= (1.0 / 6.0)) {
    throw clmassError ("Must have beta > 1/6 in mass models.");
  }

  betaModelEM.resize (m_numberOfShells);
  Real a2 = acore * acore;
  // Inner edge of outermost shell.
  size_t k = m_numberOfShells - 1;
  Real& R = rbounds[k];
  Real R2 = R * R;
  Real betaFunc = gsl_sf_beta (3.0 * beta - 0.5, 0.5);

  // Outermost annulus must be treated separately
  {
    Real rin = R;
    Real rout = rbounds [m_numberOfShells];
    Real c = 0.25 * pow (a2 + R2, 3.0 * beta) * betaFunc;
    if (beta == 0.5) {
      betaModelEM[k] = c * log ((a2 + rout * rout) / (a2 + rin * rin));
    } else {
      Real p = 1.5 - 3.0 * beta;
      betaModelEM[k] = c 
	* (pow (a2 + rout * rout, p) - pow (a2 + rin * rin, p)) / p;
    }
  }

  // Remaining annuli
  while (k != 0) {
    --k;
    Real rin = rbounds[k];
    Real rout = rbounds[k + 1];
    Real ri2 = rin * rin;
    Real ro2 = rout * rout;
    if (beta == 0.5) {
      Real t = sqrt (a2 + R2);
      Real ei = sqrt (R2 - ri2);
      Real eo = sqrt (R2 - ro2);
      betaModelEM[k] = t * t * t * (log ((t + eo) / (t + ei)) - (eo - ei) / t);
    } else {
      Real p = 1.5 - 3.0 * beta;
      Real t = a2 + R2;
      Real aigo = pow (a2 + ro2, p)
	* gsl_sf_beta_inc (3.0 * beta - 0.5, 0.5, (a2 + ro2) / t);
      Real aigi;
      // Can have a2 + ri2 == 0 for the innermost annulus
      if (a2 + ri2 > 0.0) {
	aigi = pow (a2 + ri2, p)
	  * gsl_sf_beta_inc (3.0 * beta - 0.5, 0.5, (a2 + ri2) / t);
      } else {
	// The appropriate limit
	aigi = 0.0;
      }
      // This incomplete beta function is normalized
      Real t1 = 0.5 * pow (t, 3.0 * beta) * (aigo - aigi) * betaFunc;
      t1 += t * (sqrt (R2 - ro2) - sqrt (R2 - ri2));
      betaModelEM[k] = 0.5 * t1 / p;
    }
  }

  if (tpout.maxChatter() >= chatterInfo) {
    tcout << "Beta model emission measures:";
    for (size_t i = 0; i < m_numberOfShells; ++i) {
      tcout << "\n" << rbounds[i] << " " << rbounds[i + 1] << " "
	    << betaModelEM[i];
    }
    tcout << std::endl;
  }

}

void ClusterMassModel::initialize (const std::vector<Real>& params,
                                   const IntegerArray& specNums,
				   const std::string& modelName) {

  // Check data keywords and creates Sectors for all spectra.
  // Resize, but do not compute the weight matrix.
  verifyData ();

  // Mixing model is assumed to apply only to source number 1
  static const size_t SOURCENUM = 1;
  // Data group map checked in verifyData () 
  typedef std::map<size_t, size_t> dgMapType;
  const dgMapType& dgMap = datasets->sourceToDgs ().find (SOURCENUM)->second;
  m_numberOfShells = dgMap.size ();

  // Count active spectra in each data group
  m_groupSpectra.resize (m_numberOfShells);
  for (size_t i = 0; i < m_numberOfShells; ++i) {
    m_groupSpectra [i] = 0;
  }
  size_t nTotalSpec = datasets->numberOfSpectra ();
  for (size_t iSpec = 1; iSpec <= nTotalSpec; ++iSpec) {
    const SpectralData *spec = datasets->lookup (iSpec);
    size_t groupNum = spec->parent ()->dataGroup ();
    // Spectrum active for the mixing model source?
    if (spec->detector (SOURCENUM - 1)) {
      size_t groupPosNum = dgMap.find (groupNum)->second - 1;
      m_groupSpectra [groupPosNum] += 1;
    }
  }

  // Make room for sector descriptions
  m_sector.resize (m_numberOfShells);
  m_frac.resize (m_numberOfShells);
  for (size_t grPos = 0; grPos < m_numberOfShells; ++grPos) {
    if (m_groupSpectra [grPos] == 0) {
      string diag (" no active spectra in group");
      throw IncompatibleData (diag);
    }
    m_sector [grPos].resize (m_groupSpectra [grPos]);
    m_frac [grPos].resize (m_groupSpectra [grPos]);
  }

  // If the observation 1 files are o1a1, o1a2, o1a3 and the 
  // observation 2 files are o2a1, o2a2, o2a3 then they should
  // be read in using: 
  // XSPEC> data 1:1 o1a1 1:2 o2a1 2:3 o1a2 2:4 o2a2 3:5 o1a3 3:6 o2a3
  //
  // ...but, since we're now allowing interspersed spectra that don't use
  // this model, we can no longer so easily obtain the obs number from
  // the spectrum number.  We CAN still assume obs number increases
  // consecutively within each data group for the spectra that ARE using
  // this model, so keep track of things with a counter for each 
  // data group (or annulus).

  // Create a Sector for every active spectrum for this model.
  // Also puts entries in rbounds.
  rbounds.resize (m_numberOfShells + 1);
  std::vector<size_t> obsCounter (m_numberOfShells, 0);
  for (size_t iSpec = 1; iSpec <= nTotalSpec; ++iSpec) {
    const SpectralData *spec = datasets->lookup (iSpec);
    if (spec->detector (SOURCENUM - 1)) {
      size_t groupNum = spec->parent ()->dataGroup ();
      size_t groupPosNum = dgMap.find (groupNum)->second - 1;
      size_t iObs = obsCounter [groupPosNum]++;
      createSectorDescription (iObs, iSpec, groupPosNum);
    }
  }

  // Make map from radial order to groupPosNum, as above
  radialToShell.resize (m_numberOfShells);
  for (size_t kRadial = 0; kRadial < m_numberOfShells; ++kRadial) {
    radialToShell [kRadial] = kRadial;
  }
  std::sort (radialToShell.begin (), radialToShell.end (), 
	     radialOrder (rbounds));
  
  // Complete the radial bounds and sort them
  rbounds [0] = params [0];
  std::sort (rbounds.begin (), rbounds.end ());
  if (rbounds[0] != params[0]) {
    string diag (" innermost shell is inside rinner");
    throw IncompatibleData (diag);
  }
  // Shell boundaries must be strictly increasing
  for (size_t kRadial = 0; kRadial < m_numberOfShells; ++kRadial) {
    if (rbounds[kRadial] >= rbounds[kRadial + 1]) {
      string diag (" shell boundaries overlap");
      throw IncompatibleData (diag);
    }
  }

  // Finish setting up Sectors and m_frac for active spectra by
  // setting rin from rout for the shell just inside it
  // (or rinner).
  for (size_t kRadial = 0; kRadial < m_numberOfShells; ++kRadial) {
    size_t grPos = radialToShell [kRadial];
    std::vector<Sector>& shell = m_sector [grPos];
    std::vector<Real>& fracs = m_frac [grPos];
    for (size_t iObs = 0; iObs < m_groupSpectra [grPos]; ++iObs) {
      shell [iObs].rin = rbounds [kRadial];
      // Redundant
      if (shell [iObs].rin >= shell [iObs].rout) {
	throw RedAlert ("Coding error in mass model initialize.");
      }
      fracs [iObs] = shell [iObs].angle * (1.0 / 360.0);
    }
    if (tpout.maxChatter() >= chatterDebug) {
      tcout << "Sectors for shell " << kRadial << ", group " << grPos << ":";
      for (size_t iObs = 0; iObs < m_groupSpectra [grPos]; ++iObs) {
	Sector& tmp = shell [iObs];
	tcout << "\n" << iObs << " " << rbounds [kRadial] << " " 
	      << rbounds[kRadial + 1] << " " << tmp.rin << " " << tmp.rout 
	      << " " << tmp.angle << " " << fracs [iObs];
      }
      tcout << std::endl;
    }
  }
 
  // Get storage for the mixing weight matrix
  m_matrix.resize (m_numberOfShells * m_numberOfShells, 0.0);

  // Start beta model (requires rbounds)
  useBetaModel = params[3];
  if (useBetaModel) makeBetaEM (params, true); // NB: calculation forced
}

// Accessor for matrix of mixing weights
void ClusterMassModel::setElement (size_t nshell, size_t nannulus, Real value) {
  m_matrix [m_numberOfShells * nannulus + nshell] = value;
}

// Accessor for matrix of mixing weights
Real ClusterMassModel::getElement (size_t nshell, size_t nannulus) const { 
  return m_matrix [m_numberOfShells * nannulus + nshell];
}

void ClusterMassModel::makeMatrix (const std::vector<Real>& params) {
  // Called by perform to compute the matrix of mixing weights,
  // m_matrix.  A matrix entry gives the ratio of the emission measure
  // from a shell viewed in an annulus to the emission measure of the
  // innermost shell seen in the innermost annulus.
  //
  // The distribution of gravitating mass is determined by the potential
  // model for the class derived from Hydrostatic::Atmosphere that is
  // constructed in the call to the virtual function here.  

  // Compute emission measures and mixing weights.
  Hydrostatic::Atmosphere *atmos = makeAtmos (params);
  // Finish setting up atmosphere and calculate mixing weights
  atmos->atmosProps ();

  // Use them to fill the matrix of mixing weights
  atmos->fillMatrix (*this);
  delete atmos;
  if (tpout.maxChatter() >= chatterInfo) {
    tcout << "Weight matrix:";
    for (size_t jGrp = 0; jGrp < m_numberOfShells; jGrp++) {
      tcout << "\nAnnular group " << jGrp;
      for (size_t iGrp = 0; iGrp < m_numberOfShells; iGrp++) {
	tcout << " " << getElement (iGrp, jGrp);
      }
    }
    tcout << std::endl;
  }
}

void ClusterMassModel::doPerform (const EnergyPointer& energy, 
				  GroupFluxContainer& fluxes) {

  // Weight matrix must be defined on entry.  Spectra are mixed here.
  GroupFluxContainer tmpFlux (fluxes);  
  GroupFluxContainer::iterator t (tmpFlux.begin());
  GroupFluxContainer::iterator tEnd (tmpFlux.end());
  while (t != tEnd) {
    ArrayContainer::iterator s (t->second.begin());
    ArrayContainer::iterator sEnd (t->second.end());
    while (s != sEnd) {
      // ArrayContainer == std::map<size_t,std::valarray<Real> >, so that
      // s->second has a vectorized assignment operator that takes a scalar
      s->second = 0.0;
      ++s;
    }
    ++t;         
  }

  size_t Nregions (m_numberOfShells);  
  if (tpout.maxChatter() >= chatterDebug) {
    tcout << "Annulus group, shell group, observation, mixing weight:"
	  << std::endl;
  }  
  // roi = region of interest - the data group (annulus) for which a 
  // spectral model is currently being computed.
  for (size_t roi = 1; roi <= Nregions; ++roi) {
    const EnergyPointer::const_iterator en_roi = energy.find (roi);
    const ArrayContainer& energyContainer = *(en_roi->second); 
    ArrayContainer::const_iterator e (energyContainer.begin());
    ArrayContainer& flux_roi = tmpFlux[roi];
    ArrayContainer::iterator tf (flux_roi.begin());
    ArrayContainer::iterator tfEnd (flux_roi.end());
    // Counter for observations
    size_t obs (0);
    while (tf != tfEnd) {
      // Output spectrum being modified on this pass
      // and the corresponding energies
      const RealArray& en_obs = e->second;
      RealArray& tf_obs = tf->second;
      RealArray tmp (tf_obs.size());

      // region indexes shells contributing to the roi flux
      for (size_t region = 0; region < m_numberOfShells; ++region) {
	size_t r0 (roi - 1);
	// NB: 
	Real volumeFactor = getElement (region, r0) * m_frac [r0][obs];
	if (tpout.maxChatter() >= chatterDebug) { 
	    tcout << roi << " " << region + 1 << " " << obs << " " 
		  << volumeFactor << std::endl;
	}
	if (volumeFactor > 0) {
	  // Locate the input spectrum and energy arrays for the
	  // current observation of region
	  size_t r1 (region + 1);
	  ArrayContainer& fluxr = fluxes [r1];
	  ArrayContainer::iterator ir (fluxr.begin ());
	  ArrayContainer::iterator irEnd (fluxr.end ());
	  const EnergyPointer::const_iterator en = energy.find (r1);
	  const ArrayContainer& en_calc = *(en->second);
	  ArrayContainer::const_iterator er (en_calc.begin ());
	  for (size_t i = 0; i < obs && ir != irEnd; ++i, ++ir, ++er);
	  if (ir == irEnd) {
	    // This data group has no spectrum corresponding to obs
	    // in region.  Use the spectrum for the first observation
	    // instead.
	    ir = fluxr.begin ();
	    er = en_calc.begin ();
	  }
	  RealArray& fluxir = ir->second;
	  const RealArray& en_r = er->second;
	  // Add contribution of the outer shell to the current annulus
	  if (XSutility::equalVector (en_obs, en_r)) {
	    tf_obs  += fluxir * volumeFactor;               
	  } else {
	    using namespace Numerics::Rebin;
	    size_t inputBin (0);
	    size_t outputBin (0);
	    static const Real FUZZY (1.e-06);

	    initializeBins (en_obs, en_r, FUZZY, inputBin, outputBin,
			    m_startBin, m_endBin, m_startWeight, m_endWeight); 
	    rebin (fluxir, m_startBin, m_endBin, m_startWeight, m_endWeight,
		   tmp);
	    tf_obs += volumeFactor * tmp;                                
	  }
                               
	}             
      }
      ++obs;
      ++tf;
      ++e;
    }
  }  

  if (tpout.maxChatter() >= chatterDebug) { 
    t = tmpFlux.begin();
    size_t j (1);
    using namespace std;
    ios_base::fmtflags save (tcout.flags());
    streamsize p (tcout.precision(6));
    GroupFluxContainer::iterator f (fluxes.begin());
    while (t != tEnd) {
      ArrayContainer::iterator s (t->second.begin());
      ArrayContainer::iterator ff (f->second.begin());
      ArrayContainer::iterator sEnd (t->second.end());
      const EnergyPointer::const_iterator en_j = energy.find(j);
      const ArrayContainer& energyContainer = *(en_j->second);                 
      ArrayContainer::const_iterator e (energyContainer.begin());
      tcout << " DataGroup " << t->first;
      for (;s != sEnd; ++s, ++e) {
	const RealArray& en = e->second;
	const RealArray& projFlux = s->second;
	const RealArray& unprojFlux = ff->second;
	size_t Nvec = projFlux.size();
	if (Nvec > 10) Nvec = 10;
	for (size_t k = 0; k < Nvec; ++k) {
	  tcout << "\n" << setw(4) << k << "   "  << setw(9)  << en[k+1]
		<< "   " << setw(9) << unprojFlux[k] << "   " 
		<< setw(9) << projFlux[k];
	} 
      }
      tcout << endl;
      ++t, ++f, ++j;         
    }
    tcout.flags(save);
    tcout.precision(p);
  }
  
  // This will work as a copy assignment, since the dimensions of the
  // arrays under each are identical by definition.
  fluxes = tmpFlux;   

}

void ClusterMassModel::perform (const EnergyPointer& energy, 
				const std::vector<Real>& params, 
				GroupFluxContainer& flux, 
				GroupFluxContainer& fluxError) {
  
  // The map of flux is indexed by data group (although the index is
  // groupPosNum - see initialize).  The inner map contains model
  // spectra for the different observations (numbers may differ from
  // group to group).

  // XSPEC handles throws during model calculation poorly, so
  // catch them here.
  try {
    // Compute emission measures for beta model, as needed
    makeBetaEM (params);

    // Compute emission measures and mixing weights
    makeMatrix (params);
  }
  catch (clmassError e) {
    // Inform the user
    tcout << e.message << std::endl;
    // Fill flux with NaN's to reinforce the message
    GroupFluxContainer::iterator t (flux.begin());
    GroupFluxContainer::iterator tEnd (flux.end());
    while (t != tEnd) {
      ArrayContainer::iterator s (t->second.begin());
      ArrayContainer::iterator sEnd (t->second.end());
      while (s != sEnd) {
	// Vector assigned from scalar - see comment near start of doPerform
	s->second = std::numeric_limits<Real>::quiet_NaN ();
	++s;
      }
      ++t;         
    }
    return;
  }
  
  // Mix the model spectra
  doPerform (energy, flux);
  
  if (!fluxError.empty()) doPerform (energy, fluxError);

}

void ClusterMassModel::initializeForFit (const std::vector<Real>& params, 
					 bool paramsAreFrozen) {
  // Use initialize to update rbounds or betalModelEM, if switch is
  // turned on
  IntegerArray dummy;
  if (params[0] != rbounds[0] || params[3] && !useBetaModel)
    initialize (params,dummy,string(""));
  // Updates useBetaModel in case it is turned off
  useBetaModel = params[3];
}

std::map<size_t, Real> ClusterMassModel::convertToOldXfltMap(const std::map<string,Real>& xfltMap)
{
   std::map<size_t,Real> oldStyleMap;
   std::map<string,Real>::const_iterator itXfltMap = xfltMap.begin();
   std::map<string,Real>::const_iterator itXfltMapEnd = xfltMap.end();
   while (itXfltMap != itXfltMapEnd)
   {
      // Remove "key" from front of  key string, convert the remainder to a size_t.
      string keyString = itXfltMap->first;
      keyString.erase(0,3);
      std::istringstream iss(keyString);
      size_t val=string::npos;
      iss >> val;
      oldStyleMap[val] = itXfltMap->second;
      ++itXfltMap;
   }
   return oldStyleMap;
}
