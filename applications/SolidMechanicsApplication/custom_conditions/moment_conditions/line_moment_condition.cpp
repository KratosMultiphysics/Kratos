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
#include "custom_conditions/moment_conditions/line_moment_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  LineMomentCondition::LineMomentCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : MomentCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  LineMomentCondition::LineMomentCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : MomentCondition(NewId, pGeometry, pProperties)
  {
  }

  //************************************************************************************
  //************************************************************************************
  LineMomentCondition::LineMomentCondition( LineMomentCondition const& rOther )
    : MomentCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer LineMomentCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<LineMomentCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************CLONE*******************************************
  //************************************************************************************
  Condition::Pointer LineMomentCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    LineMomentCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//
    return Kratos::make_shared<LineMomentCondition>(NewCondition);

  }


  //***********************************************************************************
  //***********************************************************************************
  LineMomentCondition::~LineMomentCondition()
  {
  }

  //************* GETTING METHODS

  //************************************************************************************
  //************************************************************************************

  void LineMomentCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    MomentCondition::InitializeConditionVariables(rVariables, rCurrentProcessInfo);

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

    //Calculate Delta Position
    //ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    //rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    //Calculate Total Delta Position
    ElementUtilities::CalculateTotalDeltaPosition(rVariables.DeltaPosition,this->GetGeometry());

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_0/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

    KRATOS_CATCH( "" )

  }


  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void LineMomentCondition::CalculateKinematics(ConditionVariables& rVariables,
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
    this->CalculateExternalMoment(rVariables);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void LineMomentCondition::CalculateExternalMoment(ConditionVariables& rVariables)
  {

    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType& dimension       = GetGeometry().WorkingSpaceDimension();

    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);

    //PRESSURE CONDITION:
    rVariables.ExternalVectorValue = rVariables.Normal;
    rVariables.ExternalScalarValue = 0.0;

    //MOMENT CONDITION:

    //defined on condition
    if( this->Has( MOMENT_LOAD ) ){
      array_1d<double, 3 > & LineLoad = this->GetValue( MOMENT_LOAD );
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoad[k];
	}
    }

    //defined on condition nodes
    if( this->Has( MOMENT_LOAD_VECTOR ) ){
      Vector& LineLoads = this->GetValue( MOMENT_LOAD_VECTOR );
      unsigned int counter = 0;
      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  counter = i*3;
	  for( SizeType k = 0; k < dimension; k++ )
	    {
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoads[counter+k];
	    }

	}
    }

    //defined on geometry nodes
    for (SizeType i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( MOMENT_LOAD ) ){
	  array_1d<double, 3 > & LineLoad = GetGeometry()[i].FastGetSolutionStepValue( MOMENT_LOAD );
	  for( SizeType k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * LineLoad[k];
	}
      }

    KRATOS_CATCH( "" )
  }


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************

  void LineMomentCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						ConditionVariables& rVariables,
						double& rIntegrationWeight)

  {
    KRATOS_TRY

    unsigned int MatSize = this->GetDofsSize();
    if(rLeftHandSideMatrix.size1() != MatSize ){
      rLeftHandSideMatrix.resize(MatSize,MatSize,false);
      noalias(rLeftHandSideMatrix) = ZeroMatrix( MatSize, MatSize );
    }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************


  int LineMomentCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = MomentCondition::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(MOMENT_LOAD);
    KRATOS_CHECK_VARIABLE_KEY(MOMENT_LOAD_VECTOR);

    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void LineMomentCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MomentCondition )
  }

  void LineMomentCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MomentCondition )
  }


} // Namespace Kratos.
