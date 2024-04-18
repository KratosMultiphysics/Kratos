// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//


// Application includes
#include "custom_conditions/U_Pw_condition.hpp"
#include "custom_utilities/dof_utilities.h"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwCondition<TDim,TNumNodes>::
    Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes * (TDim + 1);

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != conditionSize )
        rLeftHandSideMatrix.resize( conditionSize, conditionSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( conditionSize, conditionSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != conditionSize )
        rRightHandSideVector.resize( conditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( conditionSize );

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_ERROR << "UPwCondition::CalculateLeftHandSide is not implemented" << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim, TNumNodes>::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes * (TDim + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != conditionSize )
        rRightHandSideVector.resize( conditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( conditionSize );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void UPwCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    this->CalculateRHS(rRightHandSideVector, CurrentProcessInfo);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{    
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::DofsVectorType UPwCondition<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractUPwDofsFromNodes(GetGeometry(), TDim);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template class UPwCondition<2,1>;
template class UPwCondition<2,2>;
template class UPwCondition<3,1>;
template class UPwCondition<3,3>;
template class UPwCondition<3,4>;

template class UPwCondition<2,3>;
template class UPwCondition<2,4>;
template class UPwCondition<2,5>;

} // Namespace Kratos.
