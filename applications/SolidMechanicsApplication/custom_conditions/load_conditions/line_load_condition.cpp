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
#include "custom_conditions/load_conditions/line_load_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  LineLoadCondition::LineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : LoadCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!
  }

  //***********************************************************************************
  //***********************************************************************************
  LineLoadCondition::LineLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : LoadCondition(NewId, pGeometry, pProperties)
  {
    //DO NOT ADD DOFS HERE!!!
  }

  //************************************************************************************
  //************************************************************************************
  LineLoadCondition::LineLoadCondition( LineLoadCondition const& rOther )
    : LoadCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer LineLoadCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LineLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************
  Condition::Pointer LineLoadCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    LineLoadCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    return Kratos::make_shared<LineLoadCondition>(NewCondition);

  }


  //***********************************************************************************
  //***********************************************************************************
  LineLoadCondition::~LineLoadCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void LineLoadCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    LoadCondition::InitializeConditionVariables(rVariables, rCurrentProcessInfo);

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

  void LineLoadCondition::CalculateKinematics(ConditionVariables& rVariables,
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
      rVariables.Tangent1[2] = rVariables.j[rPointNumber](2, 0); // x_3,e
      rVariables.Normal[2]   = rVariables.j[rPointNumber](2, 0); // x_3,e
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
    this->CalculateExternalLoad(rVariables);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void LineLoadCondition::CalculateExternalLoad(ConditionVariables& rVariables)
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
    if( this->Has( NEGATIVE_FACE_PRESSURE ) ){
      double& NegativeFacePressure = this->GetValue( NEGATIVE_FACE_PRESSURE );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	rVariables.ExternalScalarValue += rVariables.N[i] * NegativeFacePressure;
    }

    if( this->Has( POSITIVE_FACE_PRESSURE ) ){
      double& PositiveFacePressure = this->GetValue( POSITIVE_FACE_PRESSURE );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	rVariables.ExternalScalarValue -= rVariables.N[i] * PositiveFacePressure;
    }


    //defined on condition nodes
    if( this->Has( NEGATIVE_FACE_PRESSURE_VECTOR ) ){
      Vector& Pressures = this->GetValue( NEGATIVE_FACE_PRESSURE_VECTOR );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  rVariables.ExternalScalarValue += rVariables.N[i] * Pressures[i];
	}
    }

    if( this->Has( POSITIVE_FACE_PRESSURE_VECTOR ) ){
      Vector& Pressures = this->GetValue( POSITIVE_FACE_PRESSURE_VECTOR );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  rVariables.ExternalScalarValue -= rVariables.N[i] * Pressures[i];
	}
    }

    //defined on geometry nodes
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
	  rVariables.ExternalScalarValue += rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) );
	if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
	  rVariables.ExternalScalarValue -= rVariables.N[i] * ( GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE ) );
      }

    rVariables.ExternalVectorValue *= rVariables.ExternalScalarValue;

    //FORCE CONDITION:

    //defined on condition
    if( this->Has( FORCE_LOAD ) ){
      array_1d<double, 3 > & LineLoad = this->GetValue( FORCE_LOAD );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoad[k];
	}
    }

    //defined on condition nodes
    if( this->Has( FORCE_LOAD_VECTOR ) ){
      Vector& LineLoads = this->GetValue( FORCE_LOAD_VECTOR );
      unsigned int counter = 0;
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
     counter = 3*i;
	  for( SizeType k = 0; k < dimension; k++ )
	    {
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoads[counter];
	      counter++;
	    }

	}
    }

    //defined on geometry nodes
    for (SizeType i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( FORCE_LOAD ) ){
	  array_1d<double, 3 > & LineLoad = GetGeometry()[i].FastGetSolutionStepValue( FORCE_LOAD );
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoad[k];
	}
      }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void LineLoadCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						ConditionVariables& rVariables,
						double& rIntegrationWeight)

  {
    KRATOS_TRY

      const SizeType& dimension = GetGeometry().WorkingSpaceDimension();

      if( rVariables.ExternalScalarValue == 0 )
	{
          unsigned int MatSize = this->GetDofsSize();
          if(rLeftHandSideMatrix.size1() != MatSize ){
            rLeftHandSideMatrix.resize(MatSize,MatSize,false);
            noalias(rLeftHandSideMatrix) = ZeroMatrix( MatSize, MatSize );
          }
	}
      else
	{
	  if( dimension == 2 ){

	    const SizeType number_of_nodes = GetGeometry().PointsNumber();

	    BoundedMatrix<double, 2, 2 > Kij;
	    BoundedMatrix<double, 2, 2 > SkewSymmMatrix;

	    //Compute the K sub matrix
	    SkewSymmMatrix( 0, 0 ) =  0.0;
	    SkewSymmMatrix( 0, 1 ) = -1.0;
	    SkewSymmMatrix( 1, 0 ) = +1.0;
	    SkewSymmMatrix( 1, 1 ) =  0.0;

	    double DiscretePressure;
	    unsigned int RowIndex = 0;
	    unsigned int ColIndex = 0;

	    for ( SizeType i = 0; i < number_of_nodes; i++ )
	      {
		RowIndex = i * 2;

		for ( SizeType j = 0; j < number_of_nodes; j++ )
		  {
		    ColIndex = j * 2;

		    DiscretePressure = rVariables.ExternalScalarValue * rVariables.N[i] * rVariables.DN_De( j, 0 ) * rIntegrationWeight;
		    Kij = DiscretePressure * SkewSymmMatrix;

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


  int LineLoadCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = LoadCondition::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(NEGATIVE_FACE_PRESSURE);
    KRATOS_CHECK_VARIABLE_KEY(NEGATIVE_FACE_PRESSURE_VECTOR);

    KRATOS_CHECK_VARIABLE_KEY(POSITIVE_FACE_PRESSURE);
    KRATOS_CHECK_VARIABLE_KEY(POSITIVE_FACE_PRESSURE_VECTOR);

    KRATOS_CHECK_VARIABLE_KEY(FORCE_LOAD);
    KRATOS_CHECK_VARIABLE_KEY(FORCE_LOAD_VECTOR);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void LineLoadCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LoadCondition )
  }

  void LineLoadCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LoadCondition )
  }


} // Namespace Kratos.
