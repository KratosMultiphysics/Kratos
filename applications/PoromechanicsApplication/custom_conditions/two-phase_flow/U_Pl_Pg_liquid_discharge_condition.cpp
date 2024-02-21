//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana and Danilo Cavalcanti
//


// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_liquid_discharge_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgLiquidDischargeCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgLiquidDischargeCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgLiquidDischargeCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    double DischargeScalar = this->GetGeometry()[0].FastGetSolutionStepValue( LIQUID_DISCHARGE );
    rRightHandSideVector[TDim] = DischargeScalar;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgLiquidDischargeCondition<2,1>;
template class UPlPgLiquidDischargeCondition<3,1>;

} // Namespace Kratos.
