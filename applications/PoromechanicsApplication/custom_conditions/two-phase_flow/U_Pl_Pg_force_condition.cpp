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
#include "custom_conditions/two-phase_flow/U_Pl_Pg_force_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgForceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgForceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgForceCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{    
    array_1d<double,3> ForceVector = this->GetGeometry()[0].FastGetSolutionStepValue( FORCE );
        
    for(unsigned int i = 0; i < TDim; i++)
    {
        rRightHandSideVector[i] = ForceVector[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgForceCondition<2,1>;
template class UPlPgForceCondition<3,1>;

} // Namespace Kratos.
