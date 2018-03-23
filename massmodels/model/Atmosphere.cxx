// Models a spherically symmetric, hydrostatic atmosphere.
//
// ClusterMassModel and Hydrostatic::Atmosphere are base classes
// from which models for the gravitational potential are derived.


#include <XSModel/Data/SpectralData.h>
#include <cmath>
#include "XSstreams.h"
#include "Potential.h"


const int Hydrostatic::Atmosphere::chatterInfo = 20;
const int Hydrostatic::Atmosphere::chatterDebug = 30;

// Takes array of values for the squared density and computes emission
// measure integrals for one shell in every annulus, with error estimate
Real Hydrostatic::Atmosphere::tpinteg (const std::vector<Real>& r, 
				       size_t i, size_t ny, 
				       const RealArray& yint, 
				       const RealArray& fint, 
				       RealArray& res) {

  // Integrator: does EM integrals for one spherical shell. 
  // Returns the emission measure integrals for shell i for all annuli.
  // For annulus k this is
  // $$\int_0^\sqrt{r_{i+1}^2 - r_i^2} f_i(\sqrt{r_i^2 + y^2})
  //   (\sqrt{y^2 + r_i^2 - r_k^2} - \sqrt{y^2 + r_i^2 - r_{k+1}^2}) y \, dy.$$
  // Result for annulus k is in res[k].

  RealArray last (i + 1), d2r (i + 1), rhalf (i + 1);
  // Initializations
  for (size_t k = 0; k <= i; k++) {
    res[k] = 0.0; // Integrals
    rhalf[k] = 0.0; // Integrals for half as many points
    // Integrand is always zero at y = 0.
    last[k] = 0.0; // Last evaluation of integrand k
    d2r[k] = (r[i] - r[k]) * (r[i] + r[k]);
  }

  // Three point, fourth order, integrator
  for (size_t j = 0; j < ny; j += 2) {
    // Separate off the integral for the outermost overlapping annulus
    // (k = i) to avoid redundant square roots and a negative argument.
    res[i] += last[i] + 4.0 * fint[j + 1] * yint[j + 1] * yint[j + 1];
    Real qrtu[2];
    qrtu[0] = yint[j + 1]; // Save for the next annulus inward
    Real hl = last[i]; // Save for rhalf
    last[i] = fint[j + 2] * yint[j + 2] * yint[j + 2];
    qrtu[1] = yint[j+2]; // Save for the next annulus inward
    res[i] += last[i];
    // Integral with half as many points (the test is j/2 odd)
    if (j & 2) {
      rhalf[i] += last[i];
    } else {
      rhalf[i] += hl + 4.0 * last[i];
    }

    // Integrals for all inner annuli (caution: counting down with size_t's)
    if (i != 0) {
      size_t k = i;
      do {
	--k;
	Real t = sqrt (yint[j + 1] * yint[j + 1] + d2r[k]);
	res[k] += last[k] + 4.0 * fint[j + 1] * (t - qrtu[0]) * yint[j + 1];
	qrtu[0] = t; // Saved for next annulus (k - 1), term j + 1
	t = sqrt (yint[j + 2] * yint[j + 2] + d2r[k]);
	hl = last[k]; // Save for rhalf
	// Save for the next y interval
	last[k] = fint[j + 2] * (t - qrtu[1]) * yint[j + 2];
	res[k] += last[k];
	qrtu[1] = t; // For next lower k, term j + 2
	// Integrate with half as many points for error estimate.
	// ASSUMES ny is a multiple of 4.
	if (j & 2) {
	  rhalf[k] += last[k];
	} else {
	  rhalf[k] += hl + 4.0 * last[k];
	}
      } while (k != 0);
    }
  }  

  // Normalize integrals
  Real nfac = yint[ny] / (3.0 * ny);
  Real emax = 0.0;
  for (size_t k = 0; k <= i; k++) {
    res[k] *= nfac;
    rhalf[k] *= 2.0 * nfac;
    Real eest = fabs (1.0 - rhalf[k] / res[k]);
    emax = (eest > emax)? eest: emax;
  }

  // Zero results for the non-overlapping annuli
  size_t n = res.size ();
  for (size_t k = i + 1; k < n; k++) {
    res[k] = 0.0;
  }
  return emax;
}

// Makes array of squared densities and computes EM's for one shell
Real Hydrostatic::Atmosphere::EMint (const std::vector<Real>& r, size_t i, 
				     size_t ny, RealArray& res) {
  Real a = r[i], n2 = ng2[i];
  Real dy = sqrt ((r[i + 1] - a) * (r[i + 1] + a)) / ny;
  RealArray yint (ny + 1);
  RealArray fint (ny + 1);

  // Integrate in $y = sqrt{r^2 - r_i^2}$ in order to make the
  // integrand regular (and the numerical integral fourth order).
  // Tabulate y values and squared density for integration.
  yint[0] = 0.0;
  fint[0] = n2;
  Real twoOnkT = 2.0 / kT[i];
  for (size_t j = 1; j <= ny; j++) {
    Real y = j * dy;
    Real rwk = sqrt (a * a + y * y);
    // Factor of 2 to get density squared
    Real tdpokt = twoOnkT * dphi (i, a, rwk);
    yint[j] = y;
    fint[j] = n2 * exp (- tdpokt);
  }

  return tpinteg (r, i, ny, yint, fint, res);
}

// Gets emission measures for every intersection between a shell and
// annulus
void Hydrostatic::Atmosphere::mixWeight (const ClusterMassModel& mc) {
  // Default number of points to use in integration.
  // These three are a compromise between speed and accuracy.
  // May require tuning ***
  static const size_t Nint = 32;
  // Maximum number of points to allow in integration
  static const size_t NintMax = 256;
  // Limit on acceptable error.  Assuming fourth order behaviour,
  // the actual error limit should be ~ 1/16 of this (should generally
  // do a lot better).
  static const Real Errlim = 0.0016;

  size_t n = mc.m_numberOfShells;

  for (size_t i = 0; i < n; i++) {
    // Emission measure integrals for shell i in all annuli
    RealArray& res = emissionMeasure[i];
    res.resize (n);
    if (i < n - 1 || !mc.useBetaModel) {
#if 0
      Real err = EMint (mc.rbounds, i, Nint, res);
      if (err > Errlim) {
	// Should be rare.  Estimate number of points needed to
	// reach error goal, assuming fourth order behaviour.  
	// NB: nint must be a multiple of 4.
	size_t nint = 4 * (int) ceil (0.25 * Nint * pow (err / Errlim, 0.25));
	nint = nint > NintMax ? NintMax : nint;
	err = EMint (mc.rbounds, i, nint, res);
	if (err > Errlim) {
	  pdump ();
	  throw YellowAlert ("Convergence problem in "
			     "Hydrostatic::Atmosphere::mixWeight."
			     "\nMay be due to bad model parameters.");
	}
      }
#else
      // Using a fixed number of points for the integrals means
      // no numerical glitches as Nint changes, which can upset
      // the minimization algorithm
      Real err = EMint (mc.rbounds, i, NintMax, res);
      if (err > Errlim) {
	pdump ();
	throw ClusterMassModel::clmassError (
	    "Convergence problem in Hydrostatic::Atmosphere::mixWeight."
	    "\nMay be due to bad model parameters.");
      }
#endif
    } else {
      // Beta model gives emission measures for outermost shell.
      // Note: vector times scalar.
      res = mc.betaModelEM * ng2 [n - 1];
    }
  }
  if (tpout.maxChatter() >= ClusterMassModel::chatterInfo) {
    tcout << "Emission measures (radial order)";
    for (size_t i = 0; i < n; i++) {
      tcout << "\nShell " << i;
      for (size_t j = 0; j < n; j++) {
	tcout << " " << emissionMeasure[i][j];
      }
    }
    tcout << std::endl;
  }
}

Hydrostatic::Atmosphere::Atmosphere (const std::vector<Real>& params,
				     const ClusterMassModel& mc) 
  : parent (mc) {
  // Full initialization must wait until the derived instance is ready
  size_t nShells = parent.m_numberOfShells;
  kT.resize (nShells);
  ng2.resize (nShells);
  emissionMeasure.resize (nShells);
  for (size_t i = 0; i < nShells; ++i) {
    kT [i] = params [kTBase + shellOfRadius (i)];
  }
}

Hydrostatic::Atmosphere::~Atmosphere () {
}

// Finish calculating mixing weights, after derived instance is built
void Hydrostatic::Atmosphere::atmosProps () {
  // Squared density at inner edge of each shell
  size_t n = numShells ();
  const std::vector<Real>& r = radii ();
  ng2[0] = 1.0; // Arbitrary in fact
  for (size_t i = 0; i < n - 1; ++i) {
    // 2 \Delta \phi / (kT)
    Real tdpokt = 2.0 * dphi (i, r[i], r[i + 1]) / kT[i];
    Real x = kT[i] / kT[i + 1];
    // Density variation across shell i and pressure jump to shell i + 1
     ng2[i + 1] = ng2[i] * exp (- tdpokt) * x * x;
  }
  if (tpout.maxChatter() >= ClusterMassModel::chatterInfo) {
    tcout << "Squared densities:";
    for (size_t i = 0; i < n; i++) {
      tcout << " " << ng2[i];
    }
    tcout << std::endl;
  }
  // Emission measure integrals
  mixWeight (parent); 
}

void Hydrostatic::Atmosphere::fillMatrix (ClusterMassModel& mc) {

  // Normalize to the emission measure in the intersection between
  // the innermost shell and the whole of the innermost annulus.
  // If rinner == 0, this is the emission measure for the
  // the innermost sphere.
  // Note that corrections for the fraction of an annulus that is
  // observed are in m_frac.
  Real globalNorm = 1.0 / emissionMeasure [0][0];

  size_t nshells (mc.m_numberOfShells);
  const std::vector<size_t>& radord = mc.radialToShell;
  // jAnnulus counts annuli in radial order
  for (size_t jAnnulus = 0; jAnnulus < nshells; jAnnulus++) {
    // jGrp corresponds to jAnnulus in group order
    size_t jGrp = radord[jAnnulus];
    // iShell counts shells in radial order
    for (size_t iShell = 0; iShell < nshells; iShell++) {
      // iGrp corresponds to iShell in group order
      size_t iGrp = radord[iShell];
      mc.setElement (iGrp, jGrp, globalNorm 
			     * emissionMeasure [iShell][jAnnulus]);
    }
  }
}
