// Provide physical and astronomical constants.
// This version reads the fundamental constants from a file.

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "phyconst.h"


namespace phyconst {

  // For reading values from the file of constants
  class Labelval {
  private:
    std::string label;
    double value;
    
  public:
    friend std::istream& operator>> (std::istream& sce, Labelval& res);
    friend std::ostream& operator<< (std::ostream& dest, Labelval& x);
    bool match (const char* s) {return label.find (s) == 0;}
    double val () {return value;}
  };

  std::istream& operator>> (std::istream& sce, Labelval& res) {
    // No format checks
    getline (sce, res.label, '=');
    sce >> res.value;
    // Ready for next entry
    sce.ignore (256, '\n');
    return sce;
  }

  std::ostream& operator<< (std::ostream& dest, Labelval& x) {
    dest << x.label << "= " << x.value;
    return dest;
  }


  // Physical constants
  class PAconst {
  private:
    std::ifstream *dfp; // Temporary used by constructor
    
    std::ifstream *opendata ();
    double getnext (const char *lab);

  public:
    const double c_light;
    const double hbar;
    const double G_newton;
    const double k_boltzmann;
    const double alpha_fs;
    const double q_electron;
    const double m_electron;
    const double m_proton;
    const double m_alpha;
    const double year_days;
    const double Astr_unit;
    const double GMsun;

    const double f_He, f_H;
    const double mu_0, epsilon_0;

    double mufun () const {return 1.0 / (2.0 * f_H + 3.0 * f_He * m_H / m_He);}
    double ntotfun () const {return (2.0 * f_H + 3.0 * f_He * m_H / m_He)
	/ (f_H + 2.0 * f_He * m_H / m_He);}
    double nhfun () const {return f_H / (f_H + 2.0 * f_He * m_H / m_He);}
    double rhofun () const {return m_H / (f_H + 2.0 * f_He * m_H / m_He);}
    double SBfun () const {double ttt = k_boltzmann / (c_light * hbar);
      return M_PI * M_PI / 15.0 * k_boltzmann * ttt * ttt * ttt;}
    double elecfun () const {return q_electron* q_electron
	/ (4.0 * M_PI * epsilon_0 * m_electron * c_light * c_light);}

    PAconst ();
  };

  std::ifstream* PAconst::opendata () {
    static const char fname [] = "Basic_Constants";
    std::ifstream *f = new std::ifstream (fname);
    if (!f->is_open ()) {
      std::cerr << "Failed to open " << fname << std::endl;
      exit (1);
    }
    return f;
  }

  double PAconst::getnext (const char *lab) {
    Labelval x;
    *dfp >> x;
    if (!x.match (lab)) {
      std::cerr << "Expecting " << lab << ", got " << x << std::endl;
      exit (1);
    }
    return x.val ();
  }

  PAconst::PAconst () 
    : dfp (opendata ()), 
      c_light (getnext ("Speed of light")),
      hbar (getnext ("Planck constant")),
      G_newton (getnext ("Newton constant")),
      k_boltzmann (getnext ("Boltzmann constant")),
      alpha_fs (getnext ("Fine structure constant")),
      q_electron (getnext ("Electron charge")),
      m_electron (getnext ("Electron mass")),
      m_proton (getnext ("Proton mass")),
      m_alpha (getnext ("Alpha particle mass")),
      year_days (getnext ("Tropical year (days)")),
      Astr_unit (getnext ("Astronomical unit")),
      GMsun (getnext ("GM for Sun")),
      f_He (0.25), f_H (1.0 - f_He),
      mu_0 (4e-7 * M_PI), epsilon_0 (1.0 / (mu_0 * c_light * c_light))
  {
    delete dfp;
  }

  // Read the basic constants
  static const PAconst pc;

  // Speed of light
  const double c_light = pc.c_light;

  // Planck's constant
  const double hbar = pc.hbar;

  // Newton's constant (recommended and recent)
  const double G_newton = pc.G_newton;

  // Boltmann constant
  const double k_boltzmann = pc.k_boltzmann;

  // Fine structure constant
  const double alpha_fs = pc.alpha_fs;

  // Electron
  const double q_electron = pc.q_electron;
  const double m_electron = pc.m_electron;

  // Proton and alpha particle
  const double m_proton = pc.m_proton;
  const double m_alpha = pc.m_alpha;

  // Tropical year
  const double year_days = pc.year_days;
  const double year_secs = year_days * 86400.0;

  // AU and parsec
  const double Astr_unit = pc.Astr_unit;
  const double parsec = Astr_unit * 180.0 * 3600.0 / M_PI;
  const double kiloparsec = 1000.0 * parsec;
  const double megaparsec = 1.0e6 * parsec;

  // Solar mass as GM
  const double GMsun = pc.GMsun;
  const double Msun = GMsun / G_newton;

  // Hydrogen and Helium (atomic mass deficits negligible)
  const double m_H = m_proton + m_electron;
  const double m_He = m_alpha + 2.0 * m_electron;

  // Energy conversion
  const double keV = 1000.0 * q_electron;

  // Mass fractions of helium and hydrogen
  const double f_He = pc.f_He;
  const double f_H = pc.f_H;

  // Mean mass per particle (only allows for H and He)
  const double mu = pc.mufun ();
  const double ntot_ne = pc.ntotfun ();
  const double nh_ne = pc.nhfun ();
  const double rho_ne = pc.rhofun ();

  // Electromagnetism
  const double mu_0 = pc.mu_0;
  const double epsilon_0 = pc.epsilon_0;

  // Stefan Boltzmann constants (sigma) and a */
  const double SB_density = pc.SBfun ();
  const double SB_sigma = SB_density * c_light / 4.0;

  // Classical electron radius and Thomson cross section
  const double r_electron = pc.elecfun ();
  const double sigma_thomson = 8.0 * M_PI / 3.0 * r_electron * r_electron;

}

int phyconst_main (int argc, char **argv) {
//int main (int argc, char **argv) {
  std::cout << std::setprecision (10) << "Speed of light " 
	    << phyconst::c_light << std::endl;
  std::cout << "Planck's constant, hbar " <<  phyconst::hbar
	    << "\nNewton's constant " << phyconst::G_newton
	    << "\nBoltzmann's constant " << phyconst::k_boltzmann
	    << "\nFine structure constant " << phyconst::alpha_fs
	    << "\nElectron charge " << phyconst::q_electron
	    << "\nElectron mass " << phyconst::m_electron
	    << "\nProton mass " << phyconst::m_proton
	    << "\nAlpha particle mass " << phyconst::m_alpha
	    << "\nTropical year (days) " << phyconst::year_days
	    << "\nTropical year (sec) " << phyconst::year_secs
	    << "\nAstronomical unit " << phyconst::Astr_unit
	    << "\nparsec " << phyconst::parsec
	    << "\nkiloparsec " << phyconst::kiloparsec
	    << "\nmegaparsec " << phyconst::megaparsec
	    << "\nGM for sun " << phyconst::GMsun
	    << "\nMsun " << phyconst::Msun
	    << "\nHydrogen mass " << phyconst::m_H
	    << "\nHelium mass " << phyconst::m_He
	    << "\nkeV " << phyconst::keV << "\nf_He " << phyconst::f_He
	    << "\nf_H " << phyconst::f_H << "\nmu " << phyconst::mu
	    << "\nn_tot / n_e " << phyconst::ntot_ne
	    << "\nn_H / n_e " << phyconst::nh_ne
	    << "\nrho / n_e " << phyconst::rho_ne
	    << "\nmu_0 " << phyconst::mu_0 
	    << "\nepsilon_0 " << phyconst::epsilon_0
	    << "\nStefan Boltzmann a " << phyconst::SB_density
	    << "\nStefan Boltzmann sigma " << phyconst::SB_sigma
	    << "\nClassical electron radius " << phyconst::r_electron
	    << "\nThomson cross section " << phyconst::sigma_thomson 
	    << std::endl;
  return 0;
}
