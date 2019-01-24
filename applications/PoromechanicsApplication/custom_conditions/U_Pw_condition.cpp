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
#include "custom_conditions/U_Pw_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPwCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<2,1>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 2 + 1;
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[index++] = rGeom[0].pGetDof(WATER_PRESSURE);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<2,2>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 2; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,1>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 3 + 1;
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Z);
    rConditionDofList[index++] = rGeom[0].pGetDof(WATER_PRESSURE);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,3>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 3; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,4>::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < 4; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int condition_size = TNumNodes * (TDim + 1);
    
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
void UPwCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;
    
    KRATOS_THROW_ERROR(std::logic_error,"UPwCondition::CalculateLeftHandSide not implemented","");
    
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    unsigned int condition_size = TNumNodes * (TDim + 1);
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != condition_size )
        rRightHandSideVector.resize( condition_size, false );
    noalias( rRightHandSideVector ) = ZeroVector( condition_size );
    
    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<2,1>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 + 1;
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(WATER_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 2; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,1>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 + 1;
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[index++] = rGeom[0].GetDof(WATER_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPwCondition<3,4>::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    this->CalculateRHS(rRightHandSideVector,CurrentProcessInfo);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPwCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo )
{    
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPwCondition<2,1>;
template class UPwCondition<2,2>;
template class UPwCondition<3,1>;
template class UPwCondition<3,3>;
template class UPwCondition<3,4>;

} // Namespace Kratos.
