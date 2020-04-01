//C++

#ifndef HYDROSTATICFUNC_H
#define HYDROSTATICFUNC_H

#include "xsTypes.h"
#include <XSFunctions/Utilities/XSModelFunction.h>

class MixUtility;

class HydrostaticFunc
{
   public:
      static void modFunction(const EnergyPointer& energyArray,
                                  const std::vector<Real>& parameterValues,
                                  GroupFluxContainer& flux,
                                  GroupFluxContainer& fluxError,
                                  MixUtility* mixGenerator,
                                  const std::string& modelName);
      static MixUtility* createUtility();
};

template <>
void XSCall<HydrostaticFunc>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const;

template <>
MixUtility* XSCall<HydrostaticFunc>::getUtilityObject() const;

#endif
