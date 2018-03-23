// XSPEC interface to the clmass model

#include "Monohydrostatic.h"
#include <XSModel/GlobalContainer/ModelContainer.h>
#include <XSFunctions/Utilities/funcType.h>

#include "Potential.h"


void Monohydrostatic::modFunction(const EnergyPointer& energyArray,
                                  const std::vector<Real>& parameterValues,
                                  GroupFluxContainer& flux,
                                  GroupFluxContainer& fluxError,
                                  MixUtility* mixUtility,
                                  const std::string& modelName)
{
   mixUtility->perform(energyArray,parameterValues,flux,fluxError);
}                                  

MixUtility* Monohydrostatic::createUtility()
{
   // ASSUME this can't throw.
   return new MonoClusterMassModel("MonotonicMass");
}

template <>
void XSCall<Monohydrostatic>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const
{
   Monohydrostatic::modFunction(energyArray,parameterValues,flux,fluxError, mixGenerator,modelName);
}

template <>
MixUtility* XSCall<Monohydrostatic>::getUtilityObject() const
{
   return Monohydrostatic::createUtility();
}


