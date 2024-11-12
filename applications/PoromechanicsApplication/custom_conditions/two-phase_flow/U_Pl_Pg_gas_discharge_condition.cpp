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
//


// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_gas_discharge_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgGasDischargeCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgGasDischargeCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgGasDischargeCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    double DischargeScalar = this->GetGeometry()[0].FastGetSolutionStepValue( GAS_DISCHARGE );
    rRightHandSideVector[TDim+1] = DischargeScalar;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgGasDischargeCondition<2,1>;
template class UPlPgGasDischargeCondition<3,1>;

} // Namespace Kratos.
