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
//  Main authors:    Aron Noordam
//


// Application includes
#include "custom_conditions/Pw_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer PwCondition<TDim,TNumNodes>::
    Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PwCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim,TNumNodes>::
    GetDofList(DofsVectorType& rConditionDofList,
               const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int nDof = 1;
    const unsigned int conditionSize = TNumNodes * nDof;
    if (rConditionDofList.size() != conditionSize)
        rConditionDofList.resize( conditionSize );

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rConditionDofList[i] = rGeom[i].pGetDof(WATER_PRESSURE);
    }
   

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim,TNumNodes>::
    CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;

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
void PwCondition<TDim,TNumNodes>::
    CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix,
                           const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_ERROR << "PwCondition::CalculateLeftHandSide is not implemented" << std::endl;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim, TNumNodes>::
    CalculateRightHandSide( VectorType& rRightHandSideVector,
                            const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes ;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != conditionSize )
        rRightHandSideVector.resize( conditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( conditionSize );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim,TNumNodes>::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    unsigned int nDof = 1;
    unsigned int conditionSize = TNumNodes * nDof;

    if (rResult.size() != conditionSize)
        rResult.resize( conditionSize,false );

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }
   
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim,TNumNodes>::
    CalculateAll( MatrixType& rLeftHandSideMatrix,
                  VectorType& rRightHandSideVector,
                  const ProcessInfo& CurrentProcessInfo )
{
    this->CalculateRHS(rRightHandSideVector, CurrentProcessInfo);
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void PwCondition<TDim,TNumNodes>::
    CalculateRHS(VectorType& rRightHandSideVector,
                 const ProcessInfo& CurrentProcessInfo)
{    
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
template class PwCondition<2,1>;
template class PwCondition<2,2>;
template class PwCondition<3,1>;
template class PwCondition<3,3>;
template class PwCondition<3,4>;

template class PwCondition<2,3>;

} // Namespace Kratos.
