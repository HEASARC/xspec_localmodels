// XSPEC interface to the clmass model

#include "HydrostaticFunc.h"
#include <XSModel/GlobalContainer/ModelContainer.h>
#include <XSFunctions/Utilities/funcType.h>

#include "Potential.h"


void HydrostaticFunc::modFunction(const EnergyPointer& energyArray,
                                  const std::vector<Real>& parameterValues,
                                  GroupFluxContainer& flux,
                                  GroupFluxContainer& fluxError,
                                  MixUtility* mixUtility,
                                  const std::string& modelName)
{
   mixUtility->perform(energyArray,parameterValues,flux,fluxError);
}                                  

MixUtility* HydrostaticFunc::createUtility()
{
   // ASSUME this can't throw.
   return new MIClusterMassModel("ModelIndependentMass");
}

template <>
void XSCall<HydrostaticFunc>::operator() (const EnergyPointer& energyArray, 
                        const std::vector<Real>& parameterValues, GroupFluxContainer& flux,
                        GroupFluxContainer& fluxError, MixUtility* mixGenerator, 
                        const string& modelName) const
{
   HydrostaticFunc::modFunction(energyArray,parameterValues,flux,fluxError, mixGenerator,modelName);
}

template <>
MixUtility* XSCall<HydrostaticFunc>::getUtilityObject() const
{
   return HydrostaticFunc::createUtility();
}

