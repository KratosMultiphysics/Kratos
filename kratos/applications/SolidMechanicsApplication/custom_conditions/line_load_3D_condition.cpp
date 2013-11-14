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
#include "custom_conditions/line_load_3D_condition.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::LineLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ForceLoadCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::LineLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ForceLoadCondition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoad3DCondition::LineLoad3DCondition( LineLoad3DCondition const& rOther )
    : ForceLoadCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineLoad3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************
Condition::Pointer LineLoad3DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}


//***********************************************************************************
//***********************************************************************************
LineLoad3DCondition::~LineLoad3DCondition()
{
}

//************* GETTING METHODS

//************************************************************************************
//************************************************************************************

void LineLoad3DCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

  ForceLoadCondition::InitializeGeneralVariables(rVariables, rCurrentProcessInfo);

  //Initialize Line Variables
  rVariables.Tangent1 = ZeroVector(3);
  rVariables.Normal   = ZeroVector(3);
    

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LineLoad3DCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Set Shape Functions Values for this integration point
    rVariables.N =row( Ncontainer, rPointNumber);

    //Get domain size
    rVariables.DomainSize = GetGeometry().Length();

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

Vector& LineLoad3DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //PRESSURE CONDITION:
    rVectorForce = ZeroVector(dimension);
   
    //FORCE CONDITION:
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( FACE_LOAD ) ) //temporary, will be checked once at the beginning only
	  rVectorForce += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( FACE_LOAD );
      }

    return rVectorForce;

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

double& LineLoad3DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    rIntegrationWeight *= ( GetGeometry().Length() * 0.5 );

    return rIntegrationWeight;
}


//***********************************************************************************
//***********************************************************************************


int LineLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
