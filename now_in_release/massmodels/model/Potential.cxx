// Implement some models for the gravitational potential.  Each potential
// model requires a class derived from Hydrostatic::Atmosphere.


#include <XSModel/Data/SpectralData.h>
#include <cmath>
#include "XSstreams.h"
#include "Potential.h"


// Completes construction of an MIAtmosphere by filling in the
// parameters used to compute the potential.
void Hydrostatic::MIAtmosphere::shellPars (const std::vector<Real>& r) {
  size_t n = r.size() - 1;
  mass[0] = 0.0;
  for (size_t i = 0; i < n; i++) {
    Real mca = (4.0 * M_PI / 3.0) * rho[i];
    Real r2 = r[i] * r[i];

    // Total mass within r[i + 1]
    mass[i + 1] = mass[i] + mca * (r[i + 1] - r[i]) 
      * (r[i + 1] * (r[i + 1] + r[i]) + r2);
    // Coefficients for gravitational potential.
    // In physical units, this coefficient is
    // $ G [M(r_i) - 4 \pi \rho_i r_i^3 / 3] / r_i$.
    if (r[i] > 0.0) {
      pca[i] = (mass[i] - mca * r2 * r[i]) / r[i];
    } else {
      pca[i] = 0.0;
    }
    // In physical units, this is $ 2 \pi G \rho_i / 3$
    pcb[i] = 0.5 * mca;
  }
  if (tpout.maxChatter() >= chatterInfo) {
    tcout << "Masses:";
    for (size_t i = 0; i < n; i++) {
      tcout << " " << mass[i];
    }
    tcout << std::endl;
  }
}

// Implements gravitational potential for model-independent mass model
Hydrostatic::MIAtmosphere::MIAtmosphere (const std::vector<Real>& params,
					 const ClusterMassModel& mc,
					 bool partial)
  : Atmosphere (params, mc) {
  // Model instance understands parameter layout, so do checks here
  size_t perShellParams = params.size() - kTBase;
  if (perShellParams % 2) {
    throw ClusterMassModel::clmassError (
	    "Odd number of per shell parameters - check model.dat.");
  }
  size_t nShells = numShells ();
  if (perShellParams < 2 * nShells) {
    throw ClusterMassModel::clmassError (
            "Insufficient model shells for data.  Edit model.dat.");
  }
  rho.resize (nShells);
  mass.resize (nShells + 1);
  pca.resize (nShells);
  pcb.resize (nShells);

  size_t rhoBase = kTBase + perShellParams / 2;
  for (size_t i = 0; i < nShells; ++i) {
    rho [i] = params [rhoBase + shellOfRadius (i)];
  }

  // Finish collecting info for computing the potential
  if (!partial)
    shellPars (radii ());
}

// Gravitational potential for model-independent mass distribution
double Hydrostatic::MIAtmosphere::dphi (size_t i, const double& ri, 
					const double& r) const {
  return (pca[i] / r + pcb[i] * (ri + r)) * (r - ri);
}

void Hydrostatic::MIAtmosphere::pdump () const {
  size_t nShells = numShells ();
  tcout << "Densities:";
  for (size_t i = 0; i < nShells; ++i) {
    tcout << " " << rho[i];
  }
  tcout << "\nMasses:";
  for (size_t i = 0; i <= nShells; ++i) {
    tcout << " " << mass[i];
  }
  tcout << "\n";
}

// Translate MonoAtmosphere parameters to those for MIAtmosphere.
Hydrostatic::MonoAtmosphere::MonoAtmosphere (const std::vector<Real>& params, 
					     const ClusterMassModel& mc)
  : MIAtmosphere (params, mc, true) {
  size_t nShells = numShells ();
  // Density difference for outermost shell is ignored if using the beta model
  if (params[3] != 0.0)
    --nShells;
  // Convert density differences to densities
  double denlast = 0.0;
  for (size_t i = nShells; i != 0; --i) {
    rho [i - 1] += denlast;
    denlast = rho [i - 1];
  }
  // Now finish setting parameters
  shellPars (radii ());
}


// Implements gravitational potential for model-independent mass model
Hydrostatic::NFWAtmosphere::NFWAtmosphere (const std::vector<Real>& params,
					   const ClusterMassModel& mc)
  : Atmosphere (params, mc) {
  // Model instance understands parameter layout, so do checks here
  static const size_t NFWPars = 2;
  size_t perShellParams = params.size() - kTBase - NFWPars;
  size_t nShells = numShells ();
  if (perShellParams < nShells) {
    throw ClusterMassModel::clmassError (
            "Insufficient model shells for data.  Edit model.dat.");
  }
  NFWScale = params [kTBase + perShellParams];
  NFWNorm = params [kTBase + perShellParams + 1];
}

// NFW helper function
Real Hydrostatic::NFWAtmosphere::dpaux (const Real& x) const {
  static const Real x_break = 0.001;
  if (x < x_break) {
    return 1.0 - x * (0.5 - x * ((1.0 / 3.0) - x * (0.25 - 0.2 * x)));
  }
  return log (1.0 + x) / x;
}

// NFW potential version of dphi
Real Hydrostatic::NFWAtmosphere::dphi (size_t i, const Real& ri, 
				       const Real& r) const {
  return NFWNorm * (dpaux (ri / NFWScale) - dpaux (r / NFWScale));
}

void Hydrostatic::NFWAtmosphere::pdump () const {
  tcout << "nfwa: " << NFWScale << ", nfwnorm: " << NFWNorm << "\n";
}
