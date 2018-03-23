// Derived cosmological quantities in SI units


#ifndef COSBITS_H
#define COSBITS_H 1

#include "cospar.h"

class Cosbits {
private:
  static const double arcsec;
  const double z;  // redshift
  struct Dbits {
    const Cospar *cp;
    double time;  // time
    double tbp;  // time before present
    double dcm;  // comoving distance
    double dang;  // angular diameter distance
    double dlum;  // luminosity distance
    double H;   // Hubble constant at z
    double rhocrit;  // Critical density at z
    
    Dbits (const double& zz, const char *parfile);
    ~Dbits () {delete cp;}
  } const d;

  static double angdd (const double& ocom, const double& ell, const double& a);
  static double hvar (const Cospar& pcos, const double& z);
  static void costd (const Cospar& pcos, const double& z, double& tm, 
		     double& tb, double& dcom, double& angdd, double& lumd,
		     double& H, double& rhocrit);
  //  Dbits maked (const double& zz, const char *parfile);

public:
  Cosbits (const double& zz, const char *parfile = "stdcosmo.pars")
    : z(zz), d (zz, parfile) {}
  double H0 () const {return d.cp->H_0 ();}
  double H () const {return d.H;}
  double rhocrit0 () const {return d.cp->rhocrit_0 ();}
  double rhocrit () const {return d.rhocrit;}
  double metresperarcsec () const {return d.dang * arcsec;}
  double cosmictime () const {return d.time;}
  double timebefore () const {return d.tbp;}
  double comovingdist () const {return d.dcm;}
  double angulardiamd () const {return d.dang;}
  double luminosityd () const {return d.dlum;}
};

#endif
