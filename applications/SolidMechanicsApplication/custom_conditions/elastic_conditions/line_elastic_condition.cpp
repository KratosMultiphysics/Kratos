//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/elastic_conditions/line_elastic_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  LineElasticCondition::LineElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ElasticCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  LineElasticCondition::LineElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ElasticCondition(NewId, pGeometry, pProperties)
  {
  }

  //************************************************************************************
  //************************************************************************************
  LineElasticCondition::LineElasticCondition( LineElasticCondition const& rOther )
    : ElasticCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer LineElasticCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LineElasticCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************
  Condition::Pointer LineElasticCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    LineElasticCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    return Kratos::make_shared<LineElasticCondition>(NewCondition);

  }


  //***********************************************************************************
  //***********************************************************************************
  LineElasticCondition::~LineElasticCondition()
  {
  }



  //************************************************************************************
  //************************************************************************************

  void LineElasticCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    ElasticCondition::InitializeConditionVariables(rVariables, rCurrentProcessInfo);

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    //rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    //rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    //Calculate Total Delta Position
    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_0/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )

  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void LineElasticCondition::CalculateKinematics(ConditionVariables& rVariables,
						const double& rPointNumber)
  {
    KRATOS_TRY

    const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.j[rPointNumber](0, 0); // x_1,e
    rVariables.Tangent1[1] = rVariables.j[rPointNumber](1, 0); // x_2,e

    rVariables.Normal[0] = -rVariables.j[rPointNumber](1, 0); //-x_2,e
    rVariables.Normal[1] =  rVariables.j[rPointNumber](0, 0); // x_1,e

    if(dimension==3){
      rVariables.Tangent1[2] = rVariables.J[rPointNumber](2, 0); // x_3,e
      rVariables.Normal[2]   = rVariables.J[rPointNumber](2, 0); // x_3,e
    }

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
    if(dimension==3){
      rVariables.Tangent2[2] = rVariables.J[rPointNumber](2, 0); // x_3,e
    }

    rVariables.Jacobian = norm_2(rVariables.Tangent2);

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Set Shape Functions Derivatives [dN/d£] for this integration point
    rVariables.DN_De = DN_De[rPointNumber];

    //Get geometry size
    rVariables.GeometrySize = GetGeometry().Length();

    //Get external load
    this->CalculateExternalStiffness(rVariables);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void LineElasticCondition::CalculateExternalStiffness(ConditionVariables& rVariables)
  {

    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    //noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);

    //PRESSURE CONDITION:
    rVariables.ExternalVectorValue = rVariables.Normal;
    rVariables.ExternalScalarValue = 0.0;

    //defined on condition
    if( this->Has( BALLAST_COEFFICIENT ) ){
      double& BallastCoefficient = this->GetValue( BALLAST_COEFFICIENT );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	rVariables.ExternalScalarValue += rVariables.N[i] * fabs(BallastCoefficient);
    }


    //defined on condition nodes
    if( this->Has( BALLAST_COEFFICIENT_VECTOR ) ){
      Vector& BallastCoefficient = this->GetValue( BALLAST_COEFFICIENT_VECTOR );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  rVariables.ExternalScalarValue += rVariables.N[i] * fabs(BallastCoefficient[i]);
	}
    }

    //defined on geometry nodes
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas( BALLAST_COEFFICIENT ) )
	  rVariables.ExternalScalarValue += rVariables.N[i] * fabs( GetGeometry()[i].FastGetSolutionStepValue( BALLAST_COEFFICIENT ) );
      }

    rVariables.ExternalVectorValue *= rVariables.ExternalScalarValue;


    //STIFFNESS CONDITION:

    //defined on condition
    if( this->Has( ELASTIC_LOAD ) ){
      array_1d<double, 3 > & LineStiffness = this->GetValue( ELASTIC_LOAD );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(LineStiffness[k]);
	}
    }

    //defined on condition nodes
    if( this->Has( ELASTIC_LOAD_VECTOR ) ){
      Vector& LineStiffness = this->GetValue( ELASTIC_LOAD_VECTOR );
      unsigned int counter = 0;
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( SizeType k = 0; k < dimension; k++ )
	    {
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(LineStiffness[counter+k]);
	    }

	}
    }

    //defined on geometry nodes
    for (SizeType i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( ELASTIC_LOAD ) ){
	  array_1d<double, 3 > & LineStiffness = GetGeometry()[i].FastGetSolutionStepValue( ELASTIC_LOAD );
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(LineStiffness[k]);
	}
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void LineElasticCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						 ConditionVariables& rVariables,
						 double& rIntegrationWeight)

  {
    KRATOS_TRY

      const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.ExternalScalarValue == 0 )
	{
	  ElasticCondition::CalculateAndAddKuug(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
	}
      else
	{
	  if( dimension == 2 ){

	    ElasticCondition::CalculateAndAddKuug(rLeftHandSideMatrix, rVariables, rIntegrationWeight);

	    const SizeType number_of_nodes = GetGeometry().PointsNumber();

	    BoundedMatrix<double, 2, 2 > Kij;
	    BoundedMatrix<double, 2, 2 > SkewSymmMatrix;

	    //Compute the K sub matrix
	    SkewSymmMatrix( 0, 0 ) =  0.0;
	    SkewSymmMatrix( 0, 1 ) = -1.0;
	    SkewSymmMatrix( 1, 0 ) = +1.0;
	    SkewSymmMatrix( 1, 1 ) =  0.0;

	    double coeff;
	    unsigned int RowIndex = 0;
	    unsigned int ColIndex = 0;

	    for ( SizeType i = 0; i < number_of_nodes; i++ )
	      {
		RowIndex = i * 2;

		for ( SizeType j = 0; j < number_of_nodes; j++ )
		  {
		    ColIndex = j * 2;

		    coeff = rVariables.ExternalScalarValue * rVariables.N[i] * rVariables.DN_De( j, 0 ) * rIntegrationWeight;
		    Kij = coeff * SkewSymmMatrix;

		    BeamMathUtils<double>::AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
		  }
	      }

	  }
	  else{ //3D line pressure not considered here
            unsigned int MatSize = this->GetDofsSize();
            if(rLeftHandSideMatrix.size1() != MatSize ){
              rLeftHandSideMatrix.resize(MatSize,MatSize,false);
              noalias(rLeftHandSideMatrix) = ZeroMatrix( MatSize, MatSize );
            }
	  }

	}

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************


  int LineElasticCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = ElasticCondition::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(BALLAST_COEFFICIENT);
    KRATOS_CHECK_VARIABLE_KEY(BALLAST_COEFFICIENT_VECTOR);


    KRATOS_CHECK_VARIABLE_KEY(ELASTIC_LOAD);
    KRATOS_CHECK_VARIABLE_KEY(ELASTIC_LOAD_VECTOR);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void LineElasticCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticCondition )
  }

  void LineElasticCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticCondition )
  }


} // Namespace Kratos.
