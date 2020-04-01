// Fourth order integrator

#ifndef INTEG_H
#define INTEG_H 1


class Integrable {
private:
  static const int npmin = 16;
  static const int nhmax = 36;
  static const double fleft;
  class wobbly {};

  double ecrec (const double& a, const double& b, const double& ea, 
		const double& er, const double& res, int nhlf);

public:
  virtual double integrand (const double& x) const = 0;
  double basic (const double& a, const double& b, int n);
  double ecint (const double& a, const double& b, const double& e, 
		const double& er);
  virtual ~Integrable () {}
};

#endif
