//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/point_load_3D_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{


//***********************************************************************************
//***********************************************************************************
PointLoad3DCondition::PointLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ForceLoadCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
PointLoad3DCondition::PointLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ForceLoadCondition(NewId, pGeometry, pProperties)
{

    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();

    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointLoad3DCondition::PointLoad3DCondition( PointLoad3DCondition const& rOther )
    : ForceLoadCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer PointLoad3DCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer PointLoad3DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}


//***********************************************************************************
//***********************************************************************************
PointLoad3DCondition::~PointLoad3DCondition()
{
}


//************************************************************************************
//************************************************************************************

void PointLoad3DCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rVariables.DomainSize = 1;
    rVariables.N = ZeroVector(1);
    
    //Only one node:
    rVariables.N[0] = 1;

    KRATOS_CATCH( "" )

}


//***********************************************************************************
//***********************************************************************************

Vector& PointLoad3DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    //PRESSURE CONDITION:
    rVectorForce = ZeroVector(3);
   
    //FORCE CONDITION:
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( POINT_LOAD ) ) //temporary, will be checked once at the beginning only
	  rVectorForce += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( POINT_LOAD );
      }

    return rVectorForce;

    KRATOS_CATCH( "" )
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void PointLoad3DCondition::CalculateKinematics(GeneralVariables& rVariables,
					       const double& rPointNumber)
{
    KRATOS_TRY

      rVariables.Jacobian = 1.0;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void PointLoad3DCondition::CalculateConditionSystem(LocalSystemComponents& rLocalSystem,
						    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create and initialize condition variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //reading integration points

    //force terms
    Vector VectorForce;

    for ( unsigned int PointNumber = 0; PointNumber < 1; PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = 1;
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

        if ( rLocalSystem.CalculationFlags.Is(PointLoad3DCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
            //contributions to stiffness matrix calculated on the reference config
	    this->CalculateAndAddLHS ( rLocalSystem, Variables, IntegrationWeight );
        }

        if ( rLocalSystem.CalculationFlags.Is(PointLoad3DCondition::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
            //contribution to external forces
   	    VectorForce  = this->CalculateVectorForce( VectorForce, Variables );
	 
	    this->CalculateAndAddRHS ( rLocalSystem, Variables, VectorForce, IntegrationWeight );
        }

    }

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

//***********************************************************************************
//************************************************************************************

double& PointLoad3DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{

    return rIntegrationWeight;
}


//***********************************************************************************
//***********************************************************************************


int PointLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void PointLoad3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ForceLoadCondition )
}

void PointLoad3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ForceLoadCondition )
}



} // Namespace Kratos.
