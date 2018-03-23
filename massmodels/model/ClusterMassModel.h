// Mixing model for cluster mass determinations.
//
// ClusterMassModel and Hydrostatic::Atmosphere are abstract base classes
// for models of spherical, hydrostatic atmospheres.  They implement the
// bulk of the code to mix spectra for a hydrostatic atmosphere composed
// of isothermal, spherical shells.  Classes derived from these define
// the model for the gravitational potential.


#ifndef CLUSTERMASSMODEL_H
#define CLUSTERMASSMODEL_H 1

#include "xsTypes.h"
#include <utility>
#include <XSModel/Model/Component/MixUtility.h>


namespace XSContainer {
  class DataContainer;
}

namespace Hydrostatic {
  class Atmosphere;
}


// Interface to XSPEC that applies the mixing model
class ClusterMassModel : public MixUtility {
 private:
  static const int chatterInfo;
  static const int chatterDebug;

  // A shell is a 3d region between two concentric spheres.
  // An annulus is a 2d region on the sky bounded by concentric circles.
  // A Sector defines the part of an annulus from which a spectrum was
  // extracted.
  // The radial limits of a Sector defines a corresponding shell
  // and cylinder.
  struct Sector {
    Sector (Real r_inner = ClusterMassModel::Sector::FLAG, 
	    Real r_outer = ClusterMassModel::Sector::FLAG, 
	    Real delta_angle = ClusterMassModel::Sector::FLAG);
    int operator== (const Sector &right) const;
    int operator!= (const Sector &right) const;
    Real rin;
    Real rout;
    Real angle;
    static Real FLAG;
  };

  // Number of data groups is the number of spherical shells in the model
  size_t m_numberOfShells;
  // Number of spectra for each data group (shell)
  std::vector<size_t> m_groupSpectra;
  // Shell boundaries
  std::vector<Real> rbounds;
  // Translates from radial order to shell order
  std::vector<size_t> radialToShell;
  // One vector per data group (shell) and one sector per spectrum
  std::vector<std::vector<Sector> > m_sector;
  // Fraction of the annulus observed for a spectrum - parallels m_sector
  std::vector<std::vector<Real> > m_frac;
  // Mixing weights.  Access through setElement and getElement.
  // NB: One copy for all observations with sector corrections in m_frac.
  std::vector<Real> m_matrix;

  // Used for rebinning spectra
  RealArray m_startWeight;
  RealArray m_endWeight;
  IntegerArray m_startBin;
  IntegerArray m_endBin;

  // Beta model emission measures and parameters
  RealArray betaModelEM;
  Real beta;
  Real acore;
  bool useBetaModel;

  // Used to sort shells into radial order
  class radialOrder {
  private:
    const std::vector<Real>& rb;
  public:
    radialOrder (const std::vector<Real>& rbin) : rb (rbin) {}
    bool operator() (size_t i, size_t j) const {return rb [i + 1] < rb [j + 1];}
  };

  ClusterMassModel(const ClusterMassModel &right);
  ClusterMassModel & operator= (const ClusterMassModel &right);

  // Mixes the model spectra
  void doPerform (const EnergyPointer& energy, GroupFluxContainer& fluxes);

  // Converts shell emission measures to mixing weights
  void makeMatrix (const std::vector<Real>& params);
  // Defines the gravitational potential
  virtual Hydrostatic::Atmosphere *makeAtmos (const std::vector<Real>& params)
    = 0;
  // Computes emission measures for a beta model
  void makeBetaEM (const std::vector<Real>& params, bool force = false);
  // Fills in details for one, but not rin
  void createSectorDescription (size_t iObs, size_t iSpec, size_t iShell);
  // Used by verifyData
  static Sector testSector (const std::map<size_t, Real>& xflt, 
			    size_t spectrumNumber);

  // Fortran style accessors for projection matrix
  Real getElement (size_t nshell, size_t nannulus) const;
  void setElement (size_t nshell, size_t nannulus, Real value);
  
  static std::map<size_t, Real> convertToOldXfltMap(const std::map<string,Real>& xfltMap);

 public:
  // For issues encountered during model calculation
  struct clmassError {
    string message;
    clmassError (string wha = "") : message (wha) {}
  };

  ClusterMassModel (const string& name);
  virtual ~ClusterMassModel();

  // Member functions required for a mixing model.
  // Initialize reads the input keys and uses them to create the sectors.
  virtual void initialize (const std::vector<Real>& params,
                           const IntegerArray& specNums, 
			   const std::string& modelName);
  // The matrix of miximg weights depends on model parameters, so it
  // has to be computed by perform.
  virtual void perform (const EnergyPointer& energy, 
			const std::vector<Real>& params, 
			GroupFluxContainer& flux, 
			GroupFluxContainer& fluxError);
  // In particular, called when model parameters are changed
  virtual void initializeForFit (const std::vector<Real>& params, 
				 bool paramsAreFrozen);

  protected:
     virtual void verifyData();
  friend class Hydrostatic::Atmosphere;
};

#endif
