// Turn parameters from the NFW model fit into critical radii, masses, etc

#include <iostream>
#include <cmath>
#include "phyconst.h"
#include "cosbits.h"
#include "newton.h"


// Find the value of x = r/a that maximises the mass wrt a, for fixed r
// and norm for the NFW potential (also maximises the rotation curve)
class Maxmass : public Newton {
private:
  double xmax;

public:
  Maxmass (bool noise = false) : Newton (noise), xmax (0.0) {}
  double funct (const double& x) const;
  double dfunct (const double& x) const;
  bool converged (const double& p, const double& c, const double& fp,
		  const double& dfp) const;
  double tellmax ();
};

double Maxmass::funct (const double& x) const {
  double xp= 1.0 + x, t = x / xp;
  return log (xp) - t * (1.0 + t);
}

double Maxmass::dfunct (const double& x) const {
  double onxp = 1.0 / (1.0 + x), t = x * onxp;
  return (2.0 * t - 1.0) * t * onxp;
}

bool Maxmass::converged (const double& p, const double& c, const double& fp,
			 const double& dfp) const {
  static const double eps = 1e-15;
  return fabs (fp) < eps;
}

double Maxmass::tellmax () {
  if (xmax == 0.0) {
    xmax = root (2.0);
  }
  return xmax;
}


class NFWconc : public Newton {
private:
  const double norm;
  const double target;

public:
  NFWconc (const double& rbn, const double& t, bool noise =false) 
    : Newton (noise), norm (rbn), target (t) {}
  double mform (const double& x) const;
  double dmform (const double& x) const;
  double funct (const double& x) const {
    return norm * mform (x) - target * x * x * x;
  }
  double dfunct (const double& x) const {
    return norm * dmform (x) - 3.0 * target * x * x;
  }
  bool converged (const double& p, const double& c, const double& fp,
		  const double& dfp) const;
  double dadb (const double& x) const;
};

// NFW mass distribution
double NFWconc::mform (const double& x) const {
  if (x < 1.8e-3) {
    return x * x * (0.5 - x * ((2.0 / 3.0) - x * (0.75 - 0.8 * x)));
  } else {
    return log (1.0 + x) - x / (1.0 + x);
  }
}

// Derivative of NFW mass distribution
double NFWconc::dmform (const double& x) const {
  // Should really move break
  if (x < 1.8e-3) {
    return x * (1.0 - x * (2.0 - x * (3.0 - 4.0 * x)));
  } else {
    double t = 1.0 + x;
    return x / (t * t);
  }
}

bool NFWconc::converged (const double& p, const double& c, const double& fp,
			 const double& dfp) const {
  static const double eps = 1e-14;
  return fabs (fp) < eps * target * p * p * p;
}

double NFWconc::dadb (const double& x) const {
  double t = x / (1.0 + x);
  return 1.0 - t * t / mform (x);
}


void nfwprops_usage (char **argv) {
  std::cerr << "Usage: " << argv[0] 
	    << " <z> <arcsec/scale unit> <NFW scale> <NFW norm>"
	    << " <overdensity>\n";
  exit (EXIT_FAILURE);
}

int main (int argc, char **argv) {
  enum {az = 1, aunit, asc, anorm, aover, nargs};
  if (argc != nargs)
    nfwprops_usage (argv);
  char *p;
  double z = strtod (argv[az], &p);
  if (p == argv[az])
    nfwprops_usage (argv);
  double secperunit = strtod (argv[aunit], &p);
  if (p == argv[aunit])
    nfwprops_usage (argv);
  double fitscale = strtod (argv[asc], &p);
  if (p == argv[asc])
    nfwprops_usage (argv);
  double fitnorm = strtod (argv[anorm], &p);
  if (p == argv[anorm])
    nfwprops_usage (argv);
  double overden = strtod (argv[aover], &p);
  if (p == argv[aover])
    nfwprops_usage (argv);

  Cosbits cosmo (z);

  // NFW scale in SI
  double ascale = fitscale * secperunit * cosmo.metresperarcsec ();
  std::cout << "NFW_scale: " << ascale / phyconst::kiloparsec << " kpc"
	    << std::endl;

  // NFW potential norm ($4 \pi G \rho_0 a^2$) in SI
  double potNorm = fitnorm * phyconst::keV / (phyconst::mu * phyconst::m_H);
  std::cout << "NFW_potential_norm: " << potNorm << " m^2_s^{-2}" << std::endl;

  double tmp = ascale * cosmo.H0 ();
  // Normalization for mean density / critical density
  double obarscale = 2.0 * potNorm / (tmp * tmp);
  //  std::cout << "Norm for rho_av / rho_crit: " << obarscale << std::endl;

  NFWconc conc (obarscale, overden);
  double cstart;
  if (obarscale < overden) {
    cstart = 0.5 * obarscale / overden;
  } else {
    // Crude
    cstart = pow (obarscale / overden, 1.0 / 3.0);
  }
  double cpar = conc.root (cstart);
  std::cout << "Concentration_parameter: " << cpar << std::endl;

  double rover = ascale * cpar;
  double mass = (4.0 * M_PI / 3.0) * overden * cosmo.rhocrit0 () 
    * rover * rover * rover;
  std::cout << "Enclosed_mass: " << mass / phyconst::Msun << " Msun" 
	    << std::endl;

  // Find where the mass peaks (as a function of ascale)
  //  Maxmass mp;
  //  double xmassmax = mp.tellmax ();
  //  std::cout << "Value x that maximises the mass wrt a, for fixed r: "
  //	    << xmassmax << std::endl;

  // Gradient of increasing mass from here (as da/dB)
  double grad = conc.dadb (cpar) * fitnorm / fitscale;
  std::cout << "da/dB: " << grad << std::endl;

  return EXIT_SUCCESS;
}

