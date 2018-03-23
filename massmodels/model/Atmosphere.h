// Handles most of the calculation for a hydrostatic, spherical atmosphere.
// Gravitational potentials are implemented in derived classes.


#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H 1

#include "ClusterMassModel.h"


namespace Hydrostatic {

  // Computes mixing weights for a model atmosphere
  class Atmosphere {
  private:
    const ClusterMassModel& parent;
    RealArray kT;
    RealArray ng2;
    std::vector<RealArray> emissionMeasure;

    void shellPars (const std::vector<Real>& r);
    void mixWeight (const ClusterMassModel& mc);
    Real EMint (const std::vector<Real>& r, size_t i, size_t ny, 
		RealArray& res); 
    Real tpinteg (const std::vector<Real>& r, size_t i, size_t ny, 
		  const RealArray& yint, const RealArray& fint, 
		  RealArray& res);

  protected:
    static const int chatterInfo;
    static const int chatterDebug;
    // Fixed number of parameters before the temperatures
    static const size_t kTBase = 4;
    const std::vector<Real>& radii () const {return parent.rbounds;}
    size_t numShells () const {return parent.m_numberOfShells;}
    size_t shellOfRadius (size_t i) const {return parent.radialToShell [i];}

  public:
    // Finishes making mixing weights after derived class with the 
    // gravitational potential has been built
    void atmosProps ();
    Atmosphere (const std::vector<Real>& params, const ClusterMassModel& mc);
    virtual ~Atmosphere ();
    void fillMatrix (ClusterMassModel& mc);
    // \Delta \phi - to be defined in derived classes
    virtual Real dphi (size_t i, const Real& ri, const Real& r) const = 0;
    virtual void pdump () const = 0;
  };

}

#endif
