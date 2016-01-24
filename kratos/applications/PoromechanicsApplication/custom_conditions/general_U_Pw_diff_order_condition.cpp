//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes
#include <math.h>

// Project includes
#include "custom_conditions/general_U_Pw_diff_order_condition.hpp"

#include "poromechanics_application.h"

#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"

namespace Kratos
{

// Default Constructor
GeneralUPwDiffOrderCondition::GeneralUPwDiffOrderCondition() : Condition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
GeneralUPwDiffOrderCondition::GeneralUPwDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
GeneralUPwDiffOrderCondition::GeneralUPwDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
GeneralUPwDiffOrderCondition::~GeneralUPwDiffOrderCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer GeneralUPwDiffOrderCondition::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new GeneralUPwDiffOrderCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

int  GeneralUPwDiffOrderCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if (this->Id() < 1)
    {
        KRATOS_THROW_ERROR(std::logic_error, "Condition found with Id 0 or negative","")
    }
    /* //there are some issues when computing the area of some geometries
    if (this->GetGeometry().Area() < 0)
    {
        std::cout << "error on condition -> " << this->Id() << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, "Area cannot be less than 0","")
    }
    */
    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::Initialize()
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const SizeType NumUNodes = rGeom.PointsNumber();
    
    switch(NumUNodes)
    {
        case 3: //2D L3P2
            mpPressureGeometry = GeometryType::Pointer( new Line2D2< Node<3> >(rGeom(0), rGeom(1)) );
            break;
        case 6: //3D T6P3
            mpPressureGeometry = GeometryType::Pointer( new Triangle3D3< Node<3> >(rGeom(0), rGeom(1), rGeom(2)) );
            break;
        case 8: //3D Q8P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral3D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        case 9: //3D Q9P4
            mpPressureGeometry = GeometryType::Pointer( new Quadrilateral3D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
            break;
        default:
            KRATOS_THROW_ERROR(std::logic_error,"Unexpected geometry type for different order interpolation element","");
    }
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
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
    
    for(SizeType i = 0; i < NumPNodes; i++)
    {
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( WATER_PRESSURE );
    }
    
    for(SizeType i=NumPNodes; i<NumUNodes; i++)
    {
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_X );
        rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Y );
        if(Dim > 2)
            rConditionDofList[Index++] = GetGeometry()[i].pGetDof( DISPLACEMENT_Z );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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

void GeneralUPwDiffOrderCondition::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
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
        rResult[Index++] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
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

void GeneralUPwDiffOrderCondition::InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    rVariables.NuContainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    rVariables.NpContainer = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );
    GetGeometry().Jacobian( rVariables.JContainer, mThisIntegrationMethod );
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateKinematics(ConditionVariables& rVariables,unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the shape function vector
    rVariables.Nu = row( rVariables.NuContainer, PointNumber);
    rVariables.Np = row( rVariables.NpContainer, PointNumber);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateConditionVector method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateIntegrationCoefficient method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{    

}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwDiffOrderCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAndAddConditionForce method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
