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
//  Main authors:    Mohamed Nabi
//


// Application includes
#include "custom_conditions/T_condition.hpp"

namespace Kratos
{

// ============================================================================================
// ============================================================================================
// Default constructor
template<unsigned int TDim, unsigned int TNumNodes>
TCondition<TDim, TNumNodes>::TCondition() : Condition() {}

// Constructor 1
template<unsigned int TDim, unsigned int TNumNodes>
TCondition<TDim, TNumNodes>::TCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry) {}

// Constructor 2
template<unsigned int TDim, unsigned int TNumNodes>
TCondition<TDim, TNumNodes>::TCondition(IndexType NewId, GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties) {}

// Destructor
template<unsigned int TDim, unsigned int TNumNodes>
TCondition<TDim, TNumNodes>::~TCondition() = default;

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer TCondition<TDim,TNumNodes>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new TCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;
    const GeometryType& rGeom = GetGeometry();
    if (rConditionDofList.size() != conditionSize)
        rConditionDofList.resize( conditionSize );

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rConditionDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH( "" )
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
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

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;
    const GeometryType& rGeom = GetGeometry();

    if (rResult.size() != conditionSize)
        rResult.resize( conditionSize,false );

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
    }
   
    KRATOS_CATCH( "" )
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& CurrentProcessInfo)
{
    this->CalculateRHS(rRightHandSideVector, CurrentProcessInfo);
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::CalculateRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& CurrentProcessInfo)
{    
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

// ============================================================================================
// ============================================================================================
template class TCondition<2,1>;
template class TCondition<2,2>;
template class TCondition<3,1>;
template class TCondition<3,3>;
template class TCondition<3,4>;

template class TCondition<2,3>;

} // Namespace Kratos.
