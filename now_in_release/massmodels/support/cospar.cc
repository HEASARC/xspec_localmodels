// Provide parameters for LCDM cosmology

#include "cospar.h"
#include "labelval.h"


std::ifstream* Cospar::opendata (const char* fname) {
  std::ifstream *f = new std::ifstream (fname);
  if (!f->is_open ()) {
    std::cerr << "Cospar::opendata: failed to open " << fname << "\n";
    exit (EXIT_FAILURE);
  }
  return f;
}

double Cospar::getnext (const char *lab) {
  Labelval x;
  *dfp >> x;
  if (!x.match (lab)) {
    std::cerr << "Expecting " << lab << ", got " << x << "\n";
    exit (EXIT_FAILURE);
  }
  return x.val ();
}

// Compute omega_rad from the radiation temperature.
// Plenty of assumptions here.
double Cospar::oradfn () {
  // Mass density of CMB
  double t2 = T_rad * T_rad, T4 = t2 * t2;
  double radden = phyconst::SB_density * T4 
    / (phyconst::c_light * phyconst::c_light);
  // CMB contribution to Omega_0
  double emrad = radden / rhocrit_0 ();
  // Neutrino boost to Omega_radiation
  // 3 neutrino species
  static const double nu_species = 3.0;
  // 7/8 thermal energy density of neutrinos / photons
  static const double fermion_boson = 7.0 / 8.0;
  // (4/11)^{4/3} from photon boost due to electron-positron annihilation
  static const double annihilation = pow (4.0 / 11.0, 4.0 / 3.0);
  static const double nuboost = nu_species * fermion_boson * annihilation;
  return emrad * (1.0 + nuboost);
}

Cospar::Cospar (const char* fname)
  : dfp (opendata (fname)), 
    h (getnext ("h")),
    T_rad (getnext ("T_CMB")),
    om_rad (oradfn ()),
    om_mat (getnext ("Omega_matter")),
    om_curv (getnext ("Omega_curvature")),
    om_vac (1.0 - om_curv - om_mat - om_rad)
{
  delete dfp;
}


//int main (int argc, char **argv) {
int cospar_main (int argc, char **argv) {
  if (argc !=2) {
    std::cerr << "Usage: " << argv[0] << " <cosmological parameter file>\n";
    exit (EXIT_FAILURE);
  }
  Cospar cos (argv [1]);
  double msc = phyconst::megaparsec, mpc3 = msc * msc * msc;
  std::cout << "H_0 in SI: " << cos.H_0 () << "\nH_0 in km/s/Mpc: "
	    << cos.H_0 () * phyconst::megaparsec / 1000.0 
	    << "\nh: " << cos.h_par () << "\nCritical density (SI): "
	    << cos.rhocrit_0 () << "\nOmega_radiation: " << cos.omega_rad ()
	    << "\nOmega_dust: " << cos.omega_mat () <<  "\nOmega_curvature: "
	    << cos.omega_curv () << "\nOmega_0: " << 1.0 - cos.omega_curv ()
	    << "\nOmega_vacuum: " << cos.omega_vac ()
	    << "\nMatter density at t_0 (SI): "
	    << cos.rhocdm () << "\nMatter density (Msun / Mpc^3): "
	    << cos.rhocdm () * mpc3 / phyconst::Msun << "\nMass_8 (SI): "
	    << cos.mass8 () << "\nMass_8 (Msun): " 
	    << cos.mass8 () / phyconst::Msun << std::endl;
  return EXIT_SUCCESS;
}
