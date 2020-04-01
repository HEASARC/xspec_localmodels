// Cosmic time and distances

#include <sstream>
#include "cospar.h"
#include "integ.h"
#include "cosbits.h"


// Integrand for cosmic distance
class Distint : public Integrable {
private:
  const double orad;
  const double omat;
  const double ocom;
  const double ovac;

public:
  Distint (const Cospar& pcos) : orad (pcos.omega_rad ()), 
				 omat (pcos.omega_mat ()),
				 ocom (pcos.omega_curv ()),
				 ovac (pcos.omega_vac ()) {}
  double integrand (const double& a) const;
};

double Distint::integrand (const double& a) const {
  double arg = orad + a * (omat + a * (ocom + a * a * ovac));
  if (arg <= 0.0) {
    std::cerr << "Distint::integrand: should have turned: " << a << "\n";
    exit (EXIT_FAILURE);
  }
  return 1.0 / sqrt (arg);
}


// Integrnd for cosmic time
class Timeint : public Integrable {
private:
  const double orad;
  const double omat;
  const double ocom;
  const double ovac;

public:
  Timeint (const Cospar& pcos) : orad (pcos.omega_rad ()), 
				 omat (pcos.omega_mat ()),
				 ocom (pcos.omega_curv ()),
				 ovac (pcos.omega_vac ()) {}
  double integrand (const double& a) const;
};

double Timeint::integrand (const double& a) const {
  double arg = orad + a * (omat + a * (ocom + a * a * ovac));
  if (arg <= 0.0) {
    std::cerr << "Timeint::integrand: should have turned: " << a << "\n";
    exit (EXIT_FAILURE);
  }
  return a / sqrt (arg);
}


// Dimensionless angular diameter distance from comoving distance
// and scale factor
double Cosbits::angdd (const double& ocom, const double& ell, 
		       const double& a) {
  if (ocom > 0.0) {
    // Open
    double s = sqrt (ocom), arg = s * ell;
    double t = exp (arg);
    return 0.5 * (t - 1.0 / t) * a / s;
  } else if (ocom < 0.0) {
    // Closed
    double s = sqrt (-ocom), arg = s * ell;
    return sin (arg) * a / s;
  }
  // Flat
  return a * ell;
}

// Varying part of the Hubble constant
double Cosbits::hvar (const Cospar& pcos, const double& z) {
  double ona = 1.0 + z;
  return sqrt (pcos.omega_vac () 
	       + ona * ona * (pcos.omega_curv () 
			      + ona * (pcos.omega_mat () 
				       + ona * pcos.omega_rad ())));
}

// Compute times and distances from redshift (returns SI units)
void Cosbits::costd (const Cospar& pcos, const double& z, double& time, 
		     double& tbp, double& dcom, double& dang, double& dlum,
		     double& H, double& rhocrit) {
  static const double e = 0.0, er = 1e-10;
  Timeint tint (pcos);
  double a = 1.0 / (1.0 + z), b = 1.0, orig = 0.0;
  // Time after big bang
  double tHubble = 1.0 / pcos.H_0 ();
  time = tint.ecint (orig, a, e, er) * tHubble;
  // Time before present (sec)
  tbp = tint.ecint (a, b, e, er) * tHubble;
  // Comoving distance
  Distint dint (pcos);
  double dHubble = phyconst::c_light * tHubble;
  double ell = dint.ecint (a, b, e, er);
  dcom = ell * dHubble;
  dang = angdd (pcos.omega_curv (), ell, a) * dHubble;
  dlum = dang / (a * a);
  // Hubble constant and critical density
  double hv = hvar (pcos, z);
  H = pcos.H_0 () * hv;
  rhocrit = pcos.rhocrit_0 () * hv * hv;
}

// Enables the private data to be made constant by providing the values
Cosbits::Dbits::Dbits (const double& zz, const char *parfile)
  : cp (new Cospar (parfile)) {
  costd (*cp, zz, time, tbp, dcm, dang, dlum, H, rhocrit);
}

const double Cosbits::arcsec = (M_PI / (180.0 * 3600.0));

// Test main
//int main (int argc, char **argv) {
int cosbits_main (int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <redshift>\n";
    exit (EXIT_FAILURE);
  }
  double z;
  {
    std::istringstream a (argv[1]);
    a >> z;
  }
  Cosbits cb (z);
  double Gy = 1e9 * phyconst::year_secs;
  std::cout << "Scale: " << cb.metresperarcsec () / phyconst::kiloparsec 
	    << " kpc/arcsec\nCosmic time: " << cb.cosmictime () / Gy
	    << " Gy\nTime before present: " << cb.timebefore () / Gy
	    << " Gy\nComoving distance: " 
	    << cb.comovingdist () / phyconst::megaparsec 
	    << " Mpc\nAngular diameter distance: " 
	    << cb.angulardiamd () / phyconst::megaparsec 
	    << " Mpc\nLuminosity distance: " 
	    << cb.luminosityd () / phyconst::megaparsec << " Mpc\nH_0: " 
	    << cb.H0 () * phyconst::megaparsec * 1e-3 
	    << " km s^{-1} Mpc^{-1}\nH(z): " 
	    << cb.H () * phyconst::megaparsec * 1e-3
	    << " km s^{-1} Mpc^{-1}\nrho_crit(z): " << cb.rhocrit ()
	    << " kg m^{-3}" << std::endl;
  return 0;
}
