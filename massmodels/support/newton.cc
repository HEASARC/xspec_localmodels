// Newton's method


#include "newton.h"

// Find root
double Newton::root (const double& xinit) {
  current = xinit;
  for (int i = 0; i < itmax; ++i) {
    refine ();
    if (noisy)
      std::cerr << *this << std::endl;
    if (converged (current, previous, f, df))
      return current;
  }
  throw Failed ();
}  

std::ostream& operator<< (std::ostream& dest, const Newton& newt) {
  return dest << newt.previous << ", " << newt.f << ", " << newt.df << ", " 
	      << newt.current;
}
