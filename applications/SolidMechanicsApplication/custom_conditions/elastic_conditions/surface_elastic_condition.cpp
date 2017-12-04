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
#include "custom_conditions/elastic_conditions/surface_elastic_condition.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

  //***********************************************************************************
  //***********************************************************************************
  SurfaceElasticCondition::SurfaceElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ElasticCondition(NewId, pGeometry)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  SurfaceElasticCondition::SurfaceElasticCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ElasticCondition(NewId, pGeometry, pProperties)
  {
  }

  //************************************************************************************
  //************************************************************************************
  SurfaceElasticCondition::SurfaceElasticCondition( SurfaceElasticCondition const& rOther )
    : ElasticCondition(rOther)
  {
  }

  //***********************************************************************************
  //***********************************************************************************
  Condition::Pointer SurfaceElasticCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new SurfaceElasticCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************CLONE*******************************************
  //************************************************************************************
  Condition::Pointer SurfaceElasticCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    SurfaceElasticCondition NewCondition( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewCondition.SetData(this->GetData());
    NewCondition.SetFlags(this->GetFlags());

    //-----------//      
    return Condition::Pointer( new SurfaceElasticCondition(NewCondition) );
  }


  //***********************************************************************************
  //***********************************************************************************
  SurfaceElasticCondition::~SurfaceElasticCondition()
  {
  }


  //************************************************************************************
  //************************************************************************************

  void SurfaceElasticCondition::InitializeConditionVariables(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
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

  void SurfaceElasticCondition::CalculateKinematics(ConditionVariables& rVariables,
						 const double& rPointNumber)
  {
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    
    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.J[rPointNumber](0, 0);
    rVariables.Tangent1[1] = rVariables.J[rPointNumber](1, 0);
    rVariables.Tangent1[2] = rVariables.J[rPointNumber](2, 0);

    //get second vector of the plane
    rVariables.Tangent2[0] = rVariables.J[rPointNumber](0, 1);
    rVariables.Tangent2[1] = rVariables.J[rPointNumber](1, 1);
    rVariables.Tangent2[2] = rVariables.J[rPointNumber](2, 1);

    //Compute the  normal
    CrossProduct( rVariables.Normal, rVariables.Tangent1, rVariables.Tangent2);

    //Jacobian to the last known configuration
    double Jacobian =  norm_2(rVariables.Normal);

    //auxiliar computation

    //get first vector of the plane
    rVariables.Tangent1[0] = rVariables.j[rPointNumber](0, 0);
    rVariables.Tangent1[1] = rVariables.j[rPointNumber](1, 0);
    rVariables.Tangent1[2] = rVariables.j[rPointNumber](2, 0);

    //get second vector of the plane
    rVariables.Tangent2[0] = rVariables.j[rPointNumber](0, 1);
    rVariables.Tangent2[1] = rVariables.j[rPointNumber](1, 1);
    rVariables.Tangent2[2] = rVariables.j[rPointNumber](2, 1);

    //Compute the  normal
    CrossProduct( rVariables.Normal, rVariables.Tangent1, rVariables.Tangent2);

    //Jacobian to the deformed configuration
    rVariables.Jacobian = norm_2(rVariables.Normal);

    //Compute the unit normal and weighted tangents
    if(rVariables.Jacobian>0){
      rVariables.Normal   /= rVariables.Jacobian;
      rVariables.Tangent1 /= rVariables.Jacobian;
      rVariables.Tangent2 /= rVariables.Jacobian;
    }

    //Jacobian to the last known configuration
    rVariables.Jacobian =  Jacobian;

    //Set Shape Functions Values for this integration point
    rVariables.N =row( Ncontainer, rPointNumber);

    //Set Shape Functions Derivatives [dN/d£] for this integration point
    rVariables.DN_De = DN_De[rPointNumber];

    //Get geometry size
    rVariables.GeometrySize = GetGeometry().Area();

    //Get external stiffness
    this->CalculateExternalStiffness(rVariables);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::CalculateExternalStiffness(ConditionVariables& rVariables)
  {

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    
    if( rVariables.ExternalVectorValue.size() != dimension )
      rVariables.ExternalVectorValue.resize(dimension,false);

    noalias(rVariables.ExternalVectorValue) = ZeroVector(dimension);
    
    rVariables.ExternalScalarValue = 0;
   
    //STIFFNESS CONDITION:
    
    //defined on condition
    if( this->Has( SURFACE_STIFFNESS ) ){
      array_1d<double, 3 > & SurfaceStiffness = this->GetValue( SURFACE_STIFFNESS );
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(SurfaceStiffness[k]);
	}
    }

    //defined on condition nodes
    if( this->Has( SURFACE_STIFFNESS_VECTOR ) ){
      Vector& SurfaceStiffnesss = this->GetValue( SURFACE_STIFFNESS_VECTOR );
      unsigned int counter = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
	{
	  for( unsigned int k = 0; k < dimension; k++ )
	    {
	      rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(SurfaceStiffnesss[counter]);
	      counter++;
	    }
	}
    }
        
    //defined on geometry nodes
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( SURFACE_STIFFNESS ) ){
	  array_1d<double, 3 > & SurfaceStiffness = GetGeometry()[i].FastGetSolutionStepValue( SURFACE_STIFFNESS );
	  for( unsigned int k = 0; k < dimension; k++ )
	    rVariables.ExternalVectorValue[k] += rVariables.N[i] * fabs(SurfaceStiffness[k]);
	}
      }

    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************

  void SurfaceElasticCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						 ConditionVariables& rVariables,
						 double& rIntegrationWeight)

  {
    KRATOS_TRY

      if( rVariables.ExternalScalarValue == 0 )
	{
	  ElasticCondition::CalculateAndAddKuug(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
	}
      else
	{
	  boost::numeric::ublas::bounded_matrix<double, 3, 3 > Kij;
	  boost::numeric::ublas::bounded_matrix<double, 3, 3 > Cross_ge;
	  boost::numeric::ublas::bounded_matrix<double, 3, 3 > Cross_gn;

	  double coeff;
	  const unsigned int number_of_nodes = GetGeometry().PointsNumber();

	  this->MakeCrossMatrix(Cross_ge, rVariables.Tangent1);
	  this->MakeCrossMatrix(Cross_gn, rVariables.Tangent2);

	  unsigned int RowIndex = 0;
	  unsigned int ColIndex = 0;

	  for (unsigned int i = 0; i < number_of_nodes; i++)
	    {
	      RowIndex = i * 3;
	      for (unsigned int j = 0; j < number_of_nodes; j++)
		{
		  ColIndex = j * 3;

		  coeff = rVariables.ExternalScalarValue * rVariables.N[i] * rVariables.DN_De(j, 1) * rIntegrationWeight;
		  noalias(Kij) = coeff * Cross_ge;

		  coeff = rVariables.ExternalScalarValue * rVariables.N[i] * rVariables.DN_De(j, 0) * rIntegrationWeight;

		  noalias(Kij) -= coeff * Cross_gn;

		  this->AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
		}
	    }
	}

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::MakeCrossMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M, 
					     Vector& U)
  {
    M(0, 0) = 0.00;
    M(0, 1) = U[2];
    M(0, 2) = -U[1];
    M(1, 0) = -U[2];
    M(1, 1) = 0.00;
    M(1, 2) = U[0];
    M(2, 0) = U[1];
    M(2, 1) = -U[0];
    M(2, 2) = 0.00;
  }

  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::CrossProduct(Vector & cross,
					  Vector & a,
					  Vector & b)
  {
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
  }

  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::ExpandReducedMatrix(Matrix& Destination,
						 Matrix& ReducedMatrix)
  {
    KRATOS_TRY

    unsigned int size = ReducedMatrix.size2();
    unsigned int rowindex = 0;
    unsigned int colindex = 0;

    for (unsigned int i = 0; i < size; i++)
      {
        rowindex = i * 3;
        for (unsigned int j = 0; j < size; j++)
	  {
            colindex = j * 3;
            for (unsigned int ii = 0; ii < 3; ii++)
	      Destination(rowindex + ii, colindex + ii) += ReducedMatrix(i, j);
	  }
      }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::AddMatrix(MatrixType& Destination,
				       boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InputMatrix,
				       int InitialRow,
				       int InitialCol)
  {
    KRATOS_TRY

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
	Destination(InitialRow + i, InitialCol + j) += InputMatrix(i, j);
    
    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::SubtractMatrix(MatrixType& Destination,
					    boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InputMatrix,
					    int InitialRow,
					    int InitialCol)
  {
    KRATOS_TRY

    for (unsigned int i = 0; i < 3; i++)
      for (unsigned int j = 0; j < 3; j++)
	Destination(InitialRow + i, InitialCol + j) -= InputMatrix(i, j);

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************


  int SurfaceElasticCondition::Check( const ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

    // Perform base condition checks
    int ErrorCode = 0;
    ErrorCode = ElasticCondition::Check(rCurrentProcessInfo);
      
    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(SURFACE_STIFFNESS);
    KRATOS_CHECK_VARIABLE_KEY(SURFACE_STIFFNESS_VECTOR);
        
    return ErrorCode;

    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void SurfaceElasticCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticCondition )
  }

  void SurfaceElasticCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticCondition )
  }


} // Namespace Kratos.
