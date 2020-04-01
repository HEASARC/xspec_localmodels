//C++

#ifndef MONOHYDROSTATIC_H
#define MONOHYDROSTATIC_H

#include "xsTypes.h"
#include <XSFunctions/Utilities/XSModelFunction.h>

class MixUtility;

class Monohydrostatic
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
void XSCall<Monohydrostatic>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const;

template <>
MixUtility* XSCall<Monohydrostatic>::getUtilityObject() const;

#endif
