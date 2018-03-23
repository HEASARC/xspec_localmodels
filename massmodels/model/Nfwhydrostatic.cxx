// XSPEC interface to the clmass model

#include "Nfwhydrostatic.h"
#include <XSModel/GlobalContainer/ModelContainer.h>
#include <XSFunctions/Utilities/funcType.h>

#include "Potential.h"

void Nfwhydrostatic::modFunction(const EnergyPointer& energyArray,
                                  const std::vector<Real>& parameterValues,
                                  GroupFluxContainer& flux,
                                  GroupFluxContainer& fluxError,
                                  MixUtility* mixUtility,
                                  const std::string& modelName)
{
   mixUtility->perform(energyArray,parameterValues,flux,fluxError);
}                                  

MixUtility* Nfwhydrostatic::createUtility()
{
   // ASSUME this can't throw.
   return new NFWClusterMassModel("NFWPotential");
}

template <>
void XSCall<Nfwhydrostatic>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const
{
   Nfwhydrostatic::modFunction(energyArray,parameterValues,flux,fluxError, mixGenerator,modelName);
}

template <>
MixUtility* XSCall<Nfwhydrostatic>::getUtilityObject() const
{
   return Nfwhydrostatic::createUtility();
}


