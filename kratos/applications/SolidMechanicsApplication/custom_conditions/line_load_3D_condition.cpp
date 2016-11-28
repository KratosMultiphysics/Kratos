//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/line_load_3D_condition.hpp"

#include "solid_mechanics_application_variables.h"

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
  LineLoad3DCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());

  //-----------//      
  return Condition::Pointer( new LineLoad3DCondition(NewCondition) );
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

  //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
  rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );
  
  //Calculate Delta Position
  rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

  ///calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& LineLoad3DCondition::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    rDeltaPosition.resize(number_of_nodes , dimension, false);
    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LineLoad3DCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent1[1] = rVariables.J[rPointNumber](1, 0); // x_2,e
    rVariables.Tangent1[2] = rVariables.J[rPointNumber](2, 0); // x_3,e

    //normal in the x-y plane (must be generalized)
    rVariables.Normal[0] = -rVariables.J[rPointNumber](1, 0); //-x_2,e
    rVariables.Normal[1] =  rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Normal[2] =  rVariables.J[rPointNumber](2, 0); // x_3,e

    //Jacobian to the last known configuration
    rVariables.Jacobian = norm_2(rVariables.Tangent1);

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
    
    if( rVectorForce.size() != dimension )
      rVectorForce.resize(dimension,false);

    noalias(rVectorForce) = ZeroVector(dimension);
    
    //PRESSURE CONDITION:
    rVectorForce = rVariables.Normal;
    rVariables.Pressure = 0.0;

    //defined on condition
    if( this->Has( NEGATIVE_FACE_PRESSURE ) ){
      double& NegativeFacePressure = this->GetValue( NEGATIVE_FACE_PRESSURE );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	rVariables.Pressure += rVariables.N[i] * NegativeFacePressure;
    }

    if( this->Has( POSITIVE_FACE_PRESSURE ) ){
      double& PositiveFacePressure = this->GetValue( POSITIVE_FACE_PRESSURE );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	rVariables.Pressure -= rVariables.N[i] * PositiveFacePressure;
    }

    //defined on condition nodes
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) 
	rVariables.Pressure += rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) );
      if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) 
	rVariables.Pressure -= rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE ) );     
    }
    
    rVectorForce *= rVariables.Pressure;
   
    //FORCE CONDITION:
    
    //defined on condition
    if( this->Has( LINE_LOAD ) ){
      array_1d<double, 3 > & LineLoad = this->GetValue( LINE_LOAD );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * LineLoad[k];
	}
    }

    //defined on condition nodes
    if( this->Has( LINE_LOADS_VECTOR ) ){
      Vector& LineLoads = this->GetValue( LINE_LOADS_VECTOR );
      unsigned int counter = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( unsigned int k = 0; k < dimension; k++ )
	    {
	      rVectorForce[k] += rVariables.N[i] * LineLoads[counter+k];
	    }
	  
	}
    }
    
    //defined on geometry nodes
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( LINE_LOAD ) ){
	  array_1d<double, 3 > & LineLoad = GetGeometry()[i].FastGetSolutionStepValue( LINE_LOAD );
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * LineLoad[k];
	}
      }

    //KRATOS_WATCH( rVectorForce )

    return rVectorForce;

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

int LineLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void LineLoad3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ForceLoadCondition )
}

void LineLoad3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ForceLoadCondition )
}



} // Namespace Kratos.
