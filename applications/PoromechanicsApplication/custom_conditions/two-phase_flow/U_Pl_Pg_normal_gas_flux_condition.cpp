//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//                   Danilo Cavalcanti
//


// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_normal_gas_flux_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgNormalGasFluxCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgNormalGasFluxCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgNormalGasFluxCondition<TDim,TNumNodes>::GetNormalFluidFluxVector(array_1d<double,TNumNodes>& rNormalFluidFluxVector, const GeometryType& Geom)
{
    for(unsigned int i=0; i<TNumNodes; i++)
    {
        // Multiplied by -1.0 to obtain an inflow from a positive value of NORMAL_GAS_FLUX
        rNormalFluidFluxVector[i] = -1.0*Geom[i].FastGetSolutionStepValue(NORMAL_GAS_FLUX);
    }
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgNormalGasFluxCondition<TDim,TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector, NormalFluidFluxVariables& rVariables)
{
    noalias(rVariables.PVector) = -rVariables.NormalFluidFlux * rVariables.Np * rVariables.IntegrationCoefficient;
    
    PoroElementUtilities::AssemblePgBlockVector< array_1d<double,TNumNodes> >(rRightHandSideVector,rVariables.PVector,TDim,TNumNodes);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgNormalGasFluxCondition<2,2>;
template class UPlPgNormalGasFluxCondition<3,3>;
template class UPlPgNormalGasFluxCondition<3,4>;

} // Namespace Kratos.
