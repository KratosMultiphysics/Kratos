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

#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"

// Project includes
#include "custom_conditions/one-phase_flow/general_U_Pl_diff_order_condition.hpp"

namespace Kratos
{

// Default Constructor
GeneralUPlDiffOrderCondition::GeneralUPlDiffOrderCondition() : Condition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
GeneralUPlDiffOrderCondition::GeneralUPlDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
GeneralUPlDiffOrderCondition::GeneralUPlDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = this->GetIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
GeneralUPlDiffOrderCondition::~GeneralUPlDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer GeneralUPlDiffOrderCondition::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new GeneralUPlDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();

    switch(NumUNodes)
    {
        case 3: //2D L3P2
            mpPressureGeometry = GeometryType::Pointer( new Line2D2< Node >(rGeom(0), rGeom(1)) );
            break;
        case 6: //3D T6P3
            mpPressureGeometry = GeometryType::Pointer( new Triangle3D3< Node >(rGeom(0), rGeom(1), rGeom(2)) );
            break;
        case 8: //3D Q8P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral3D4< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 9: //3D Q9P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral3D4< Node >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    if(rConditionDofList.size() != ConditionSize)
        rConditionDofList.resize(ConditionSize);

    SizeType Index = 0;
/*
    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( LIQUID_PRESSURE );
    }

    for(SizeType i=NumPNodes; i<NumUNodes; i++)
    {
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }
*/

    for(SizeType i = 0; i < NumUNodes; i++)
    {
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }

    for(SizeType i=0; i<NumPNodes; i++)
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( LIQUID_PRESSURE );

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ConditionSize )
        rLeftHandSideMatrix.resize( ConditionSize, ConditionSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ConditionSize, ConditionSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_THROW_ERROR(std::logic_error,"GeneralUPlDiffOrderCondition::CalculateLeftHandSide not implemented","");

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::EquationIdVector(EquationIdVectorType& rResult,const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const GeometryType& rGeom = GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType ConditionSize = NumUNodes * Dim + NumPNodes;

    if ( rResult.size() != ConditionSize )
        rResult.resize( ConditionSize, false );

    SizeType Index = 0;

    for ( SizeType i = 0; i < NumUNodes; i++ )
    {
        rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if(Dim > 2)
            rResult[Index++] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    for ( SizeType i = 0; i < NumPNodes; i++ )
        rResult[Index++] = GetGeometry()[i].GetDof( LIQUID_PRESSURE ).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo,
                                        bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    //Definition of variables
    ConditionVariables Variables;
    this->InitializeConditionVariables(Variables,rCurrentProcessInfo);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics (Np)
        this->CalculateKinematics(Variables,PointNumber);

        //Compute Condition Vector
        this->CalculateConditionVector(Variables,PointNumber);

        //Calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( Variables, PointNumber, integration_points[PointNumber].Weight() );

        //Contributions to the left hand side
        if ( CalculateLHSMatrixFlag )
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if ( CalculateResidualVectorFlag )
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGPoints = rGeom.IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType WorkingDim = rGeom.WorkingSpaceDimension();
    const SizeType LocalDim = rGeom.LocalSpaceDimension();

    (rVariables.NuContainer).resize(NumGPoints,NumUNodes,false);
    rVariables.NuContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );

    (rVariables.NpContainer).resize(NumGPoints,NumPNodes,false);
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

    (rVariables.Nu).resize(NumUNodes,false);
    (rVariables.Np).resize(NumPNodes,false);

    (rVariables.JContainer).resize(NumGPoints,false);
    for(SizeType i = 0; i<NumGPoints; i++)
        ((rVariables.JContainer)[i]).resize(WorkingDim,LocalDim,false);
    rGeom.Jacobian( rVariables.JContainer, mThisIntegrationMethod );
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateKinematics(ConditionVariables& rVariables,unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the shape function vector
    noalias(rVariables.Nu) = row( rVariables.NuContainer, PointNumber);
    noalias(rVariables.Np) = row( rVariables.NpContainer, PointNumber);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateConditionVector method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateIntegrationCoefficient method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{

}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPlDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAndAddConditionForce method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
