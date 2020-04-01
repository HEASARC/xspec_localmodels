// Provide cosmological parameters for LCDM

#ifndef COSPAR_H
#define COSPAR_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include "phyconst.h"

class Cospar {
private:
  std::ifstream *dfp; // Temporary used by constructor

  const double h;  // H_0 / 100 km/s/Mpc
  const double T_rad;  // Radiation temperature (K)
  // Contributions to Omega_0
  const double om_rad;  // Radiation and neutrinos
  const double om_mat;  // Dust (baryons plus dark matter)
  const double om_curv;  // 1 - Omega_0
  const double om_vac;  // Vacuum energy

  std::ifstream *opendata (const char* fname);
  double getnext (const char *lab);
  double oradfn ();

public:
  Cospar (const char* fname);
  double H_0 () const {return h * 1e5 / phyconst::megaparsec;}
  // Critical density at t_0
  double rhocrit_0 () const {
    double H0 = H_0 ();
    return 3.0 * H0 * H0 / (8.0 * M_PI * phyconst::G_newton);
  }
  // Density at t_0 of CDM (dark matter and baryons)
  double rhocdm () const {return om_mat * rhocrit_0 ();}
  // Mass corresponding to sigma_8
  double mass8 () const {
    double r8 = 8.0 * phyconst::megaparsec / h, r3 = r8 * r8 * r8;
    return (4.0 * M_PI / 3.0) * rhocdm () * r3;
  }
  const double& h_par () const {return h;}
  const double& omega_rad () const {return om_rad;}
  const double& omega_mat () const {return om_mat;}
  const double& omega_curv () const {return om_curv;}
  const double& omega_vac () const {return om_vac;}
};

#endif
