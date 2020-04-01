// Fourth order integrator

#include <iostream>
#include <cmath>
#include <sstream>
#include "integ.h"


// Instanced for safety (should not be necessary)
const int Integrable::npmin;
const int Integrable::nhmax;
const double Integrable::fleft = 0.5;

// Basic fourth order integrator
double Integrable::basic (const double& a, const double& b, int n) {
  // n must be even
  if (n <= 0 || n % 2 == 1) {
    std::cerr << "Integrable::basic: bad n " << n << std::endl;
    throw wobbly ();
  }
  double last = integrand (a);
  double step = (b - a) / n;
  double sum = 0.0;
  for (int i = 1; i < n; i += 2) {
    sum += last + 4.0 * integrand (a + i * step);
    last = integrand (a + (i + 1) * step);
    sum += last;
  }
  return (1.0 / 3.0) * step * sum;
}

// Recursive integrator used by ecint
double Integrable::ecrec (const double& a, const double& b, const double& ea, 
			  const double& er, const double& res, int nhlf) {
  double mid = 0.5 * (a + b);
  double resl = basic (a, mid, npmin);
  double resr = basic (mid, b, npmin);
  // ntot += 2 * npmin;
  if (fabs (resl + resr - res) <= ea + er * fabs (res)) {
    return resl + resr;
  }
  if (nhlf >= nhmax) {
    std::cerr << "Integrable::ecrec: too many halvings" << std::endl;
    throw wobbly ();
  }
  ++nhlf;
  // Apportion absolute error (not relative error)
  double eal = fleft * ea;
  double ear = ea - eal;
  return ecrec (a, mid, eal, er, resl, nhlf) 
    + ecrec (mid, b, ear, er, resr, nhlf);
}

// Integrator with error control
double Integrable::ecint (const double& a, const double& b, const double& e, 
			  const double& er) {
  double res = basic (a, b, npmin);
  // ntot = npmin;
  int nhlf = 0;
  res = ecrec (a, b, e, er, res, nhlf);
  // std::cout << "Used " << ntot << " calls" << std::endl;
  return res;
}


// Test example

class Intex : public Integrable {
public:
  double integrand (const double& x) const {return sin (x);}
};

// int main (int argc, char **argv) {
int integ_main (int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] << " <a> <b> <e> <er>" << std::endl;
    exit (1);
  }
  double a;
  {
    std::istringstream ast (argv[1]);
    ast >> a;
  }
  a *= M_PI;
  double b;
  {
    std::istringstream ast (argv[2]);
    ast >> b;
  }
  b *= M_PI;
  double e;
  {
    std::istringstream ast (argv[3]);
    ast >> e;
  }
  double er;
  {
    std::istringstream ast (argv[4]);
    ast >> er;
  }
  Intex f;
  double res = f.ecint (a, b, e, er);
  double perf = cos (a) - cos (b);
  std::cout << res << ", " << perf << ", " << res - perf << std::endl;
  return 0;
}
