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
  PointLoad3DCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  NewCondition.SetData(this->GetData());
  NewCondition.SetFlags(this->GetFlags());

  //-----------//      
  return Condition::Pointer( new PointLoad3DCondition(NewCondition) );
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
      
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int local_dimension = GetGeometry().LocalSpaceDimension();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVariables.Initialize(dimension, local_dimension, number_of_nodes);
   
    //Only one node:
    rVariables.N[0] = 1.0;

    KRATOS_CATCH( "" )

}


//***********************************************************************************
//***********************************************************************************

Vector& PointLoad3DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVectorForce.size() != dimension )
      rVectorForce.resize(dimension,false);

    rVectorForce = ZeroVector(dimension);
   
    //FORCE CONDITION:
    //defined on condition
    if( this->Has( POINT_LOAD ) ){
      array_1d<double, 3 > & PointLoad = this->GetValue( POINT_LOAD );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * PointLoad[k];
	}
    }
    
    //defined on condition nodes
    if( this->Has( POINT_LOADS_VECTOR ) ){
      Vector& PointLoads = this->GetValue( POINT_LOADS_VECTOR );
      unsigned int counter = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( unsigned int k = 0; k < dimension; k++ )
	    {
	      rVectorForce[k] += rVariables.N[i] * PointLoads[counter+k];
	    }
	  
	}
    }
    
    //defined on condition nodes      
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( POINT_LOAD ) ){
	  array_1d<double, 3 > & PointLoad = GetGeometry()[i].FastGetSolutionStepValue( POINT_LOAD );
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVectorForce[k] += rVariables.N[i] * PointLoad[k];

// 	  std::cout<<" PointLoad "<<PointLoad<<" N "<<rVariables.N[0]<<std::endl;
 
	}
      }

//     std::cout<<" rVerctorForce "<<rVectorForce<<std::endl;

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
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Vector VectorForce = ZeroVector(dimension);

    for ( unsigned int PointNumber = 0; PointNumber < 1; PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = 1;

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
