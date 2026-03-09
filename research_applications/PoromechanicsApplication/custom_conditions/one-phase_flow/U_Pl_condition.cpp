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
#include "custom_conditions/one-phase_flow/U_Pl_condition.hpp"

namespace Kratos
{

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer UPlCondition<TDim,TNumNodes>::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPlCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<2,1>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 2 + 1;
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[index++] = rGeom[0].pGetDof(LIQUID_PRESSURE);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<2,2>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    for (unsigned int i = 0; i < 2; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(LIQUID_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,1>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const unsigned int condition_size = 3 + 1;
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_X);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Y);
    rConditionDofList[index++] = rGeom[0].pGetDof(DISPLACEMENT_Z);
    rConditionDofList[index++] = rGeom[0].pGetDof(LIQUID_PRESSURE);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,3>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    for (unsigned int i = 0; i < 3; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(LIQUID_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,4>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;

    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );

    for (unsigned int i = 0; i < 4; i++)
    {
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Y);
        rConditionDofList[index++] = rGeom[i].pGetDof(DISPLACEMENT_Z);
        rConditionDofList[index++] = rGeom[i].pGetDof(LIQUID_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlCondition<TDim,TNumNodes>::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
void UPlCondition<TDim,TNumNodes>::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    KRATOS_THROW_ERROR(std::logic_error,"UPlCondition::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlCondition<TDim,TNumNodes>::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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
void UPlCondition<2,1>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 + 1;
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(LIQUID_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 2 * (2 + 1);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 2; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,1>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 + 1;
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Y).EquationId();
    rResult[index++] = rGeom[0].GetDof(DISPLACEMENT_Z).EquationId();
    rResult[index++] = rGeom[0].GetDof(LIQUID_PRESSURE).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 3 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 3; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

template< >
void UPlCondition<3,4>::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    unsigned int condition_size = 4 * (3 + 1);
    unsigned int index = 0;

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < 4; i++)
    {
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = rGeom[i].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[index++] = rGeom[i].GetDof(LIQUID_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlCondition<TDim,TNumNodes>::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    this->CalculateRHS(rRightHandSideVector,rCurrentProcessInfo);
}

//----------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlCondition<TDim,TNumNodes>::CalculateRHS( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateRHS method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template< unsigned int TDim, unsigned int TNumNodes >
void UPlCondition<TDim,TNumNodes>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    const Variable<array_1d<double,3> >& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    GeometryType& rGeom = GetGeometry();

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL ) {
        // CD, VV

        for(SizeType i=0; i< TNumNodes; ++i) {

            SizeType index = (TDim + 1) * i;
            array_1d<double, 3 >& r_external_force = rGeom[i].FastGetSolutionStepValue(EXTERNAL_FORCE);

            for(SizeType j=0; j<TDim; ++j) {
                #pragma omp atomic
                r_external_force[j] += rRHSVector[index + j];
            }
        }
    } else if(rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == REACTION ) {
        // Residual/Reactions

        for(SizeType i=0; i< TNumNodes; ++i) {
            SizeType index = (TDim + 1) * i;
            array_1d<double, 3 >& r_force_residual = rGeom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            double& r_flux_residual = rGeom[i].FastGetSolutionStepValue(FLUX_RESIDUAL);

            for(SizeType j=0; j<TDim; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j];
            }

            #pragma omp atomic
            r_flux_residual += rRHSVector[index + TDim];
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template class UPlCondition<2,1>;
template class UPlCondition<2,2>;
template class UPlCondition<3,1>;
template class UPlCondition<3,3>;
template class UPlCondition<3,4>;

} // Namespace Kratos.
