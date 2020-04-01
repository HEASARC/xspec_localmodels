// Newton's method as a virtual base class

#ifndef NEWTON_H
#define NEWTON_H 1

#include <iostream>

class Newton {
private:
  static const int itmax = 30; // Deliberately excessive

  const bool noisy;
  double current, previous, f, df;

  void refine () { 
    f = funct (current); 
    df = dfunct (current); 
    previous = current;
    current -= f / df;
  }

  class Failed {};

public:
  Newton (bool noise) : noisy (noise) {}
  virtual ~Newton () {}
  virtual double funct (const double& x) const = 0;
  virtual double dfunct (const double& x) const = 0;
  virtual bool converged (const double& p, const double& c, const double& fp,
			  const double& dfp) const = 0;
  double root (const double& xinit);
  friend std::ostream& operator<< (std::ostream& dest, const Newton& newt);
};

#endif
