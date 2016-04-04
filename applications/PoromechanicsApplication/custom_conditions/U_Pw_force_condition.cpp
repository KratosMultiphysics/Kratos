//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

// Application includes
#include "custom_conditions/U_Pw_force_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwForceCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwForceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< >
void UPwForceCondition<2,1>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    array_1d<double,3> ForceVector = this->GetGeometry()[0].FastGetSolutionStepValue( FORCE );
    const double& Thickness = this->GetProperties()[THICKNESS];
    
    rRightHandSideVector[0] = ForceVector[0] * Thickness;
    rRightHandSideVector[1] = ForceVector[1] * Thickness;
}

//----------------------------------------------------------------------------------------

template< >
void UPwForceCondition<3,1>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    array_1d<double,3> ForceVector = this->GetGeometry()[0].FastGetSolutionStepValue( FORCE );
    
    rRightHandSideVector[0] = ForceVector[0];
    rRightHandSideVector[1] = ForceVector[1];
    rRightHandSideVector[2] = ForceVector[2];
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwForceCondition<2,1>;
template class UPwForceCondition<3,1>;

} // Namespace Kratos.
