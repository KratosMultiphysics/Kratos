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
//                   Danilo Cavalcanti
//


// Application includes
#include "custom_conditions/two-phase_flow/U_Pl_Pg_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlPgCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlPgCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int condition_size = TNumNodes * (TDim + 2);
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    for (unsigned int i = 0; i < TNumNodes; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        if constexpr (TDim>2)
            rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(LIQUID_PRESSURE);
        rConditionDofList[index++] = rGeom[i].pGetDof(GAS_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes * (TDim + 2);

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != condition_size )
        rLeftHandSideMatrix.resize( condition_size, condition_size, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( condition_size, condition_size );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_THROW_ERROR(std::logic_error,"UPlPgCondition::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes * (TDim + 2);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgCondition<2,1>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 + 2;
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(LIQUID_PRESSURE).EquationId();
    rResult[index++] = rGeom[0].GetDof(GAS_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 * (2 + 2);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 2; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
        rResult[index++] = rGeom[i].GetDof(GAS_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgCondition<3,1>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 + 2;
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[index++] = rGeom[0].GetDof(LIQUID_PRESSURE).EquationId();
    rResult[index++] = rGeom[0].GetDof(GAS_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 2);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
        rResult[index++] = rGeom[i].GetDof(GAS_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlPgCondition<3,4>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 2);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
        rResult[index++] = rGeom[i].GetDof(GAS_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    this->CalculateRHS(rRightHandSideVector,rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlPgCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlPgCondition<2,1>;
template class UPlPgCondition<2,2>;
template class UPlPgCondition<3,1>;
template class UPlPgCondition<3,3>;
template class UPlPgCondition<3,4>;

} // Namespace Kratos.
