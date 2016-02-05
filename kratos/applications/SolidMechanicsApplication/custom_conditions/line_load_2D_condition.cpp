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
#include "custom_conditions/line_load_2D_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
LineLoad2DCondition::LineLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : LineLoad3DCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
LineLoad2DCondition::LineLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LineLoad3DCondition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoad2DCondition::LineLoad2DCondition( LineLoad2DCondition const& rOther )
    : LineLoad3DCondition(rOther)     
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineLoad2DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new LineLoad2DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************
Condition::Pointer LineLoad2DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}


//***********************************************************************************
//***********************************************************************************
LineLoad2DCondition::~LineLoad2DCondition()
{
}

//************* GETTING METHODS

//************************************************************************************
//************************************************************************************

void LineLoad2DCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
  
  LineLoad3DCondition::InitializeGeneralVariables(rVariables, rCurrentProcessInfo);

}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void LineLoad2DCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.j[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent1[1] = rVariables.j[rPointNumber](1, 0); // x_2,e

    rVariables.Normal[0] = -rVariables.j[rPointNumber](1, 0); //-x_2,e
    rVariables.Normal[1] =  rVariables.j[rPointNumber](0, 0); // x_1,e

    //Jacobian to the deformed configuration
    rVariables.Jacobian = norm_2(rVariables.Normal);

    //std::cout<< " jacobian "<<rVariables.Jacobian<<std::endl;

    //Compute the unit normal and weighted tangents
    if(rVariables.Jacobian>0){
      rVariables.Normal   /= rVariables.Jacobian;
      rVariables.Tangent1 /= rVariables.Jacobian;
    }


    //Jacobian to the last known configuration
    rVariables.Tangent2[0] = rVariables.J[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent2[1] = rVariables.J[rPointNumber](1, 0); // x_2,e

    rVariables.Jacobian = norm_2(rVariables.Tangent2);

    //std::cout<< " Jacobian "<<rVariables.Jacobian<<std::endl;

    //Set Shape Functions Values for this integration point
    rVariables.N =row( Ncontainer, rPointNumber);

    //Get domain size
    rVariables.DomainSize = GetGeometry().Length();
    
    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

Vector& LineLoad2DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    //const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //PRESSURE CONDITION:
    rVectorForce = rVariables.Normal;
    rVariables.Pressure = 0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      if( GetGeometry()[j].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) && GetGeometry()[j].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) //temporary, will be checked once at the beginning only
	rVariables.Pressure += rVariables.N[j] * ( GetGeometry()[j].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) - GetGeometry()[j].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE ) );
      
    }
    
    rVectorForce *= rVariables.Pressure;
   
    //FORCE CONDITION:
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( LINE_LOAD ) ) //temporary, will be checked once at the beginning only
	  rVectorForce += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( LINE_LOAD );
      }

    //KRATOS_WATCH( rVectorForce )

    return rVectorForce;

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

double& LineLoad2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    if( GetProperties()[THICKNESS] > 0 )
      rIntegrationWeight *=  GetProperties()[THICKNESS];

  return rIntegrationWeight;
}

//***********************************************************************************
//***********************************************************************************

void LineLoad2DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
					      GeneralVariables& rVariables,
					      double& rIntegrationWeight)

{
    KRATOS_TRY

    if( rVariables.Pressure == 0 )
      {
	rLeftHandSideMatrix = ZeroMatrix( 4, 4 );
      }
    else
      {
	const unsigned int number_of_nodes = GetGeometry().PointsNumber();
	
	Matrix Kij     ( 2, 2 );
	Matrix SkewSymmMatrix( 2, 2 );
	
	//Compute the K sub matrix
	SkewSymmMatrix( 0, 0 ) =  0.0;
	SkewSymmMatrix( 0, 1 ) = -1.0;
	SkewSymmMatrix( 1, 0 ) = +1.0;
	SkewSymmMatrix( 1, 1 ) =  0.0;

	double DiscretePressure=0;
	unsigned int RowIndex = 0;
	unsigned int ColIndex = 0;

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    RowIndex = i * 2;
	    
	    for ( unsigned int j = 0; j < number_of_nodes; j++ )
	      {
		ColIndex = j * 2;
		
		DiscretePressure = rVariables.Pressure * rVariables.N[i] * rVariables.DN_De( j, 0 ) * rIntegrationWeight;
		Kij = DiscretePressure * SkewSymmMatrix;
		
		MathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
	      }
	  }
      }

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************


int LineLoad2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void LineLoad2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoad3DCondition )
}

void LineLoad2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoad3DCondition )
}


} // Namespace Kratos.
