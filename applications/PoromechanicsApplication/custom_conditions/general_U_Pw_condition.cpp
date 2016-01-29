//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes
#include <math.h>

// Project includes
#include "custom_conditions/general_U_Pw_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
GeneralUPwCondition::GeneralUPwCondition() : Condition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
GeneralUPwCondition::GeneralUPwCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
GeneralUPwCondition::GeneralUPwCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
GeneralUPwCondition::~GeneralUPwCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer GeneralUPwCondition::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new GeneralUPwCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * (dimension + 1);
    unsigned int index = 0;
    
    if (rConditionDofList.size() != condition_size)
      rConditionDofList.resize( condition_size );
    
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rConditionDofList[index++] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        rConditionDofList[index++] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        if( dimension > 2 )
            rConditionDofList[index++] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
        
        rConditionDofList[index++] = GetGeometry()[i].pGetDof(WATER_PRESSURE);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dimension + 1);
    
    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != MatSize )
        rLeftHandSideMatrix.resize( MatSize, MatSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );
    
    //Resetting the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( MatSize );
    
    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.PointsNumber();
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * (dimension + 1);
    unsigned int index = 0;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        rResult[index++] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        if( dimension > 2)
            rResult[index++] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

        rResult[index++] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const GeometryType& rGeom = GetGeometry();
    const unsigned int number_of_nodes = rGeom.size();
    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dimension + 1);
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( MatSize );
    
    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
    
    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
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

void GeneralUPwCondition::InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& rGeom = GetGeometry();
    rVariables.NContainer = rGeom.ShapeFunctionsValues( mThisIntegrationMethod );
    rGeom.Jacobian( rVariables.JContainer, mThisIntegrationMethod );
    
    //FIC condition variables
    double BulkModulus = GetProperties()[YOUNG_MODULUS]/(3.0*(1.0-2.0*GetProperties()[POISSON_RATIO]));
    double BulkModulusSolid = GetProperties()[BULK_MODULUS_SOLID];
    double BiotCoefficient = 1.0-BulkModulus/BulkModulusSolid;
    double Porosity = GetProperties()[POROSITY];
    rVariables.BiotModulusInverse = (BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/GetProperties()[BULK_MODULUS_FLUID];

    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    rVariables.NewmarkCoefficient = 1.0/(rCurrentProcessInfo[THETA_NEWMARK]*DeltaTime);

    const unsigned int dimension = rGeom.WorkingSpaceDimension();
    if(dimension==2)
        rVariables.ElementLength = rGeom.Length();
    else
        rVariables.ElementLength = sqrt(4.0*rGeom.Area()/M_PI);
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateKinematics(ConditionVariables& rVariables,unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the shape function vector
    rVariables.Np = row( rVariables.NContainer, PointNumber);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateConditionVector method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateIntegrationCoefficient method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{    

}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralUPwCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAndAddConditionForce method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
