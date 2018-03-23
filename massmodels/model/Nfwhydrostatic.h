//C++

#ifndef NFWHYDROSTATIC_H
#define NFWHYDROSTATIC_H

#include "xsTypes.h"
#include <XSFunctions/Utilities/XSModelFunction.h>

class MixUtility;

class Nfwhydrostatic
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
void XSCall<Nfwhydrostatic>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const;

template <>
MixUtility* XSCall<Nfwhydrostatic>::getUtilityObject() const;

#endif
