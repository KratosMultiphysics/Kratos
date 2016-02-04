//
//   Project Name:         $
//   Last modified by:    $Author:             $
//   Date:                $Date:                $
//   Revision:            $Revision:        $
//

// System includes
#include <math.h>

// Project includes
#include "custom_conditions/general_condition.hpp"

#include "dam_application.h"

namespace Kratos
{

// Default Constructor
GeneralCondition::GeneralCondition() : Condition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
GeneralCondition::GeneralCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
GeneralCondition::GeneralCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
GeneralCondition::~GeneralCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer GeneralCondition::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new GeneralCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

void GeneralCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < GetGeometry().size(); i++)
    {
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        if( dimension == 3 )
            rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void GeneralCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dimension;
    
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

void GeneralCondition::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * dimension;
    int index;
    
    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * dimension;
        
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

void GeneralCondition::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
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
        //compute element kinematics (N)
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

void GeneralCondition::InitializeConditionVariables (ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    rVariables.NContainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    GetGeometry().Jacobian( rVariables.JContainer, mThisIntegrationMethod );
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateKinematics(ConditionVariables& rVariables,unsigned int PointNumber)
{
    KRATOS_TRY

    //Setting the shape function vector
    rVariables.N = row( rVariables.NContainer, PointNumber);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateConditionVector(ConditionVariables& rVariables, unsigned int PointNumber)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateConditionVector method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateIntegrationCoefficient(ConditionVariables& rVariables, unsigned int PointNumber, double weight)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateIntegrationCoefficient method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables)
{    

}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateAndAddRHS(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    this->CalculateAndAddConditionForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void GeneralCondition::CalculateAndAddConditionForce(VectorType& rRightHandSideVector, ConditionVariables& rVariables)
{
    KRATOS_TRY

    KRATOS_THROW_ERROR( std::logic_error, "calling the default CalculateAndAddConditionForce method for a particular condition ... illegal operation!!", "" )

    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
