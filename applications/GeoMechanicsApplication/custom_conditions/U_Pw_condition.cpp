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

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwCondition<TDim,TNumNodes>::
    Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    GetDofList(DofsVectorType& rConditionDofList,
               const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int nDof = TDim + 1;
    const unsigned int conditionSize = TNumNodes * nDof;
    if (rConditionDofList.size() != conditionSize)
        rConditionDofList.resize( conditionSize );

    if constexpr (TDim == 2) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
            rConditionDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
        }
    } else if constexpr (TDim == 3) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
            rConditionDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
        }
    } else {
        KRATOS_ERROR << "undefined dimension in U_Pw_condition!!" << std::endl;
    }

    KRATOS_CATCH( "" )
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
template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::
    EquationIdVector(EquationIdVectorType& rResult,
                     const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();

    unsigned int nDof = TDim + 1;
    unsigned int conditionSize = TNumNodes * nDof;

    if (rResult.size() != conditionSize)
        rResult.resize( conditionSize );

    if constexpr (TDim == 2) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
        }
    } else if constexpr (TDim == 3) {
        unsigned int index = 0;
        for (unsigned int i = 0; i < TNumNodes; ++i) {
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
            rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
        }
    } else {
        KRATOS_ERROR << "Undefined dimension in U_Pw_condition!!" << std::endl;
    }

    KRATOS_CATCH( "" )
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
