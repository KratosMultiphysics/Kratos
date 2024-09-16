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
#include "custom_conditions/U_Pw_discharge_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwDischargeCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwDischargeCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwDischargeCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    double DischargeScalar = this->GetGeometry()[0].FastGetSolutionStepValue( DISCHARGE );
    rRightHandSideVector[TDim] = DischargeScalar;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwDischargeCondition<2,1>;
template class UPwDischargeCondition<3,1>;

} // Namespace Kratos.
