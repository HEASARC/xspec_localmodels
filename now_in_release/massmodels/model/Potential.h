// Potentials for use with the mass models


#ifndef POTENTIAL_H
#define POTENTIAL_H 1


#include "Atmosphere.h"


namespace Hydrostatic {

  // Implements gravitational potential for "model-independent" mass
  // determinations.  The potential model has constant mass density in each
  // shell.
  class MIAtmosphere : public Atmosphere {
  protected:
    RealArray rho; // Gravitating mass densities

  private:
    RealArray mass; // Cumulative masses
    RealArray pca;
    RealArray pcb;

  public:
    MIAtmosphere (const std::vector<Real>& params, const ClusterMassModel& mc,
		  bool partial = false);
    void shellPars (const std::vector<Real>& r);
    Real dphi (size_t i, const Real& ri, const Real& r) const;
    void pdump () const;
  };


  // Same as MIAtmosphere, but using gravitating mass density differences
  // as parameters.  Makes the monotonic constraint a lot easier.
  class MonoAtmosphere : public MIAtmosphere {
  public:
    MonoAtmosphere (const std::vector<Real>& params, 
		    const ClusterMassModel& mc);
  };


  // Implements NFW gravitational potential
  class NFWAtmosphere : public Atmosphere {
  private:
    Real NFWScale;
    Real NFWNorm;

  public:
    NFWAtmosphere (const std::vector<Real>& params, const ClusterMassModel& mc);
    Real dphi (size_t i, const Real& ri, const Real& r) const;
    Real dpaux (const Real& x) const;
    void pdump () const;
  };

}


// The following declarations derive concrete classes from ClusterMassModel
// to use the potentials declared above

// The "model-independent" mass model
class MIClusterMassModel : public ClusterMassModel {
 public:
  MIClusterMassModel (const string& name) : ClusterMassModel (name) {}
  Hydrostatic::Atmosphere *makeAtmos (const std::vector<Real>& params) {
    return new Hydrostatic::MIAtmosphere (params, *this);
  }
};


// The monotonic "model-independent" mass model
class MonoClusterMassModel : public ClusterMassModel {
 public:
  MonoClusterMassModel (const string& name) : ClusterMassModel (name) {}
  Hydrostatic::Atmosphere *makeAtmos (const std::vector<Real>& params) {
    return new Hydrostatic::MonoAtmosphere (params, *this);
  }
};


// NFW mass model
class NFWClusterMassModel : public ClusterMassModel {
 public:
  NFWClusterMassModel (const string& name) : ClusterMassModel (name) {}
  Hydrostatic::Atmosphere *makeAtmos (const std::vector<Real>& params) {
    return new Hydrostatic::NFWAtmosphere (params, *this);
  }
};

#endif
