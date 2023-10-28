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
//                   John van Esch
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

    if (rConditionDofList.size() != TNumNodes) {
        rConditionDofList.resize(TNumNodes);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rConditionDofList[i] = rGeom[i].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH( "" )
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim, TNumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //Resetting the LHS
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }
    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    //Resetting the RHS
    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }
    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim, TNumNodes>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != TNumNodes) {
        rResult.resize(TNumNodes, false);
    }

    const GeometryType& rGeom = GetGeometry();
    for (unsigned int i = 0; i < TNumNodes; ++i) {
        rResult[i] = rGeom[i].GetDof(TEMPERATURE).EquationId();
    }

    KRATOS_CATCH("")
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
}

// ============================================================================================
// ============================================================================================
template<unsigned int TDim, unsigned int TNumNodes>
void TCondition<TDim,TNumNodes>::CalculateRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{    
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular condition ... illegal operation!!" << std::endl;

    KRATOS_CATCH( "" )
}

// ============================================================================================
// ============================================================================================
template class TCondition<2,1>;
template class TCondition<2,2>;
template class TCondition<2,3>;
template class TCondition<2,4>;
template class TCondition<2,5>;
template class TCondition<3,1>;
template class TCondition<3,3>;
template class TCondition<3,4>;
template class TCondition<3,6>;
template class TCondition<3,8>;
template class TCondition<3,9>;

} // Namespace Kratos.
