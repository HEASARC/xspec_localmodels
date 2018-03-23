// Angular diameter distance and critical density

#include <iostream>
#include <cstdlib>
#include "cosbits.h"


void usage (char **argv) {
  std::cerr << "Usage: " << argv[0] << " <cosmological par file> <redshift>\n";
  exit (EXIT_FAILURE);
}


int main (int argc, char **argv) {
  if (argc != 3)
    usage (argv);

  char *p;
  double z = strtod (argv[2], &p);
  if (p == argv[2])
    usage (argv);

  Cosbits cb (z, argv[1]);
  std::cout.precision (12);
  std::cout << "Angular diameter distance: " << cb.angulardiamd ()
	    << " m\nCritical density: " << cb.rhocrit () 
	    << " kg m^{-3}" << std::endl;

  return EXIT_SUCCESS;
}
