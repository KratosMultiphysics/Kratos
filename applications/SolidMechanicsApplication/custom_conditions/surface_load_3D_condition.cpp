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
#include "custom_conditions/surface_load_3D_condition.hpp"

#include "solid_mechanics_application.h"


namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : ForceLoadCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : ForceLoadCondition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
SurfaceLoad3DCondition::SurfaceLoad3DCondition( SurfaceLoad3DCondition const& rOther )
    : ForceLoadCondition(rOther)
{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer SurfaceLoad3DCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new SurfaceLoad3DCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}


//************************************CLONE*******************************************
//************************************************************************************
Condition::Pointer SurfaceLoad3DCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return (this->Create( NewId, rThisNodes, pGetProperties() ) );
}


//***********************************************************************************
//***********************************************************************************
SurfaceLoad3DCondition::~SurfaceLoad3DCondition()
{
}

//************* GETTING METHODS

//************************************************************************************
//************************************************************************************

void SurfaceLoad3DCondition::InitializeGeneralVariables(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

  ForceLoadCondition::InitializeGeneralVariables(rVariables, rCurrentProcessInfo);
  
  //Initialize Surface Variables
  rVariables.Tangent1 = ZeroVector(3);
  rVariables.Tangent2 = ZeroVector(3);
  rVariables.Normal   = ZeroVector(3);

  //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
  rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );


  //Calculate Delta Position
  rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

  ///calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}

//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& SurfaceLoad3DCondition::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

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

void SurfaceLoad3DCondition::CalculateKinematics(GeneralVariables& rVariables,
					     const double& rPointNumber)
{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //auxiliar computation

    //std::cout<<" PointNumber "<<rPointNumber<<" rVarables.J "<<rVariables.J[rPointNumber]<<std::endl;

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

    //Get domain size
    rVariables.DomainSize = GetGeometry().Area();
    

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

Vector& SurfaceLoad3DCondition::CalculateVectorForce(Vector& rVectorForce, GeneralVariables& rVariables)
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
      else if( GetGeometry()[j].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) 
	rVariables.Pressure += rVariables.N[j] * ( GetGeometry()[j].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) );
      else if( GetGeometry()[j].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) 
	rVariables.Pressure -= rVariables.N[j] * ( GetGeometry()[j].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE ) );
      
    }
    
    rVectorForce *= rVariables.Pressure;

    //FORCE CONDITION:
    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
	if( GetGeometry()[i].SolutionStepsDataHas( SURFACE_LOAD ) ) //temporary, will be checked once at the beginning only
	  rVectorForce += rVariables.N[i] * GetGeometry()[i].FastGetSolutionStepValue( SURFACE_LOAD );
      }


    //KRATOS_WATCH( rVectorForce )

      
    return rVectorForce;

    KRATOS_CATCH( "" )
}


//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************

double& SurfaceLoad3DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{ 
    return rIntegrationWeight;
}


//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
						 GeneralVariables& rVariables,
						 double& rIntegrationWeight)

{
    KRATOS_TRY

    if( rVariables.Pressure == 0 )
	{
		//rLeftHandSideMatrix = ZeroMatrix( 9, 9 ); // This assumes a triangle!!!
		// no need to do anything. LHS is already set to zero from force_load_condition.InitializeSystemMatrices
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

				coeff = rVariables.Pressure * rVariables.N[i] * rVariables.DN_De(j, 1) * rIntegrationWeight;
				noalias(Kij) = coeff * Cross_ge;

				coeff = rVariables.Pressure * rVariables.N[i] * rVariables.DN_De(j, 0) * rIntegrationWeight;

				noalias(Kij) -= coeff * Cross_gn;

				this->AddMatrix( rLeftHandSideMatrix, Kij, RowIndex, ColIndex );
			}
		}
	}

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::MakeCrossMatrix(boost::numeric::ublas::bounded_matrix<double, 3, 3 > & M, 
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

void SurfaceLoad3DCondition::CrossProduct(Vector & cross,
				      Vector & a,
				      Vector & b)
{
    cross[0] = a[1] * b[2] - a[2] * b[1];
    cross[1] = a[2] * b[0] - a[0] * b[2];
    cross[2] = a[0] * b[1] - a[1] * b[0];
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::ExpandReducedMatrix(Matrix& Destination,
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

void SurfaceLoad3DCondition::AddMatrix(MatrixType& Destination,
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

void SurfaceLoad3DCondition::SubtractMatrix(MatrixType& Destination,
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


int SurfaceLoad3DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************

void SurfaceLoad3DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ForceLoadCondition )
}

void SurfaceLoad3DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ForceLoadCondition )
}


} // Namespace Kratos.
