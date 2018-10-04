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
#include "custom_conditions/thermal_conditions/line_heat_flux_condition.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
LineHeatFluxCondition::LineHeatFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  this->Set(THERMAL);
}

//************************************************************************************
//************************************************************************************
LineHeatFluxCondition::LineHeatFluxCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
  this->Set(THERMAL);
}

Condition::Pointer LineHeatFluxCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
  return Condition::Pointer( new LineHeatFluxCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

LineHeatFluxCondition::~LineHeatFluxCondition()
{
}


//************************************************************************************
//************************************************************************************
void LineHeatFluxCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  //calculation flags
  bool CalculateStiffnessMatrixFlag = false;
  bool CalculateResidualVectorFlag = true;
  MatrixType temp = Matrix();

  CalculateElementalSystem( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void LineHeatFluxCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
  //calculation flags
  bool CalculateStiffnessMatrixFlag = true;
  bool CalculateResidualVectorFlag = true;

  CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

}

//************************************************************************************
//************************************************************************************
void LineHeatFluxCondition::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                                    ProcessInfo& rCurrentProcessInfo,
                                                    bool CalculateStiffnessMatrixFlag,
                                                    bool CalculateResidualVectorFlag )
{
  KRATOS_TRY

  unsigned int number_of_nodes = GetGeometry().size();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  //resizing as needed the LHS
  unsigned int MatSize = number_of_nodes;

  if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
  {
    if ( rLeftHandSideMatrix.size1() != MatSize )
      rLeftHandSideMatrix.resize( MatSize, MatSize, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
  }

  //resizing as needed the RHS
  if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
  {
    if ( rRightHandSideVector.size( ) != MatSize )
      rRightHandSideVector.resize( MatSize, false );

    noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
  }

  //reading integration points and local gradients
  const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

  const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();

  const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

  //calculating actual jacobian
  GeometryType::JacobiansType J;

  J = GetGeometry().Jacobian( J );

  //FLUX CONDITION:
  Vector FluxOnNodes = ZeroVector( number_of_nodes );
  for ( unsigned int i = 0; i < number_of_nodes; i++ )
  {
    FluxOnNodes[i]=GetGeometry()[i].FastGetSolutionStepValue( FACE_HEAT_FLUX );
  }

  Vector NormalVector  = ZeroVector( 2 ); //normal direction (not normalized)
  Vector TangentVector = ZeroVector( 2 ); //normal direction (not normalized)

  for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
  {
    //calculating tangent
    TangentVector[0] =  J[PointNumber]( 0, 0 ); // x_1,e
    TangentVector[1] =  J[PointNumber]( 1, 0 ); // x_2,e

    double IntegrationWeight = integration_points[PointNumber].Weight() * norm_2(TangentVector);

    if ( dimension == 2 && GetProperties()[THICKNESS]>0 )
      IntegrationWeight *= GetProperties()[THICKNESS];

    //calculating normal
    NormalVector[0] = -J[PointNumber]( 1, 0 ); //-x_2,e
    NormalVector[1] =  J[PointNumber]( 0, 0 ); // x_1,e

    //calculating the flux on the gauss point
    double gauss_flux = 0.00;
    Vector Flux=ZeroVector( 3 );
    for ( unsigned int ii = 0; ii < number_of_nodes; ii++ )
    {
      gauss_flux += Ncontainer( PointNumber, ii ) * FluxOnNodes[ii];

    }

    if ( CalculateStiffnessMatrixFlag == true )
    {
      if ( gauss_flux != 0.00 )
      {
        CalculateAndSubKheatflux( rLeftHandSideMatrix, DN_De[PointNumber], row( Ncontainer, PointNumber ), gauss_flux, IntegrationWeight );
      }
    }

    //adding contributions to the residual vector
    if ( CalculateResidualVectorFlag == true )
    {
      if ( gauss_flux != 0.00 )
        CalculateAndAddFaceHeatFlux( rRightHandSideVector, row( Ncontainer, PointNumber ), NormalVector, gauss_flux, IntegrationWeight );

    }

  }

  KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
void LineHeatFluxCondition::CalculateAndSubKheatflux( Matrix& rK,
                                                    const Matrix& rDN_De,
                                                    const Vector& rN,
                                                    double rFlux,
                                                    double rIntegrationWeight)
{
  KRATOS_TRY

  const unsigned int number_of_nodes = GetGeometry().PointsNumber();

  Matrix Kij     ( 2, 2 );
  Matrix SkewSymmMatrix( 2, 2 );

  //Compute the K sub matrix
  SkewSymmMatrix( 0, 0 ) =  0.0;
  SkewSymmMatrix( 0, 1 ) = -1.0;
  SkewSymmMatrix( 1, 0 ) = -1.0;
  SkewSymmMatrix( 1, 1 ) =  0.0;

  double DiscreteFlux=0;

  for ( unsigned int i = 0; i < number_of_nodes; i++ )
  {
    unsigned int RowIndex = i;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
      unsigned int ColIndex = j;

      DiscreteFlux = rFlux * rN[i] * rDN_De( j, 0 ) * rIntegrationWeight;
      Kij = - DiscreteFlux * SkewSymmMatrix;

      //TAKE CARE: the load correction matrix should be SUBSTRACTED not added
      MathUtils<double>::SubtractMatrix( rK, Kij, RowIndex, ColIndex );
    }
  }


  KRATOS_CATCH( "" )
}


//***********************************************************************
//***********************************************************************
void LineHeatFluxCondition::CalculateAndAddFaceHeatFlux(Vector& rF,
                                                      const Vector& rN,
                                                      Vector& rNormal,
                                                      double rFlux,
                                                      double rIntegrationWeight )
{
  unsigned int number_of_nodes = GetGeometry().size();

  for ( unsigned int i = 0; i < number_of_nodes; i++ )
  {
    double DiscreteFlux = rFlux * rN[i] * rIntegrationWeight;
    rF[i]  += DiscreteFlux;
  }

}

//************************************************************************************
//************************************************************************************
void LineHeatFluxCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
  int number_of_nodes     = GetGeometry().PointsNumber(); //=Geometry().size();

  unsigned int conditiondimension = number_of_nodes;

  if ( rResult.size() != conditiondimension )
    rResult.resize( conditiondimension, false );

  for ( int i = 0; i < number_of_nodes; i++ )
  {
    rResult[i]   =  GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
  }

}

//************************************************************************************
//************************************************************************************
void LineHeatFluxCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rConditionalDofList.resize( 0 );
  for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
  {
    rConditionalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
  }
}

//************************************************************************************
//************************************************************************************
void  LineHeatFluxCondition::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

      if ( rMassMatrix.size1() != 0 )
        rMassMatrix.resize( 0, 0, false );

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void  LineHeatFluxCondition::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

      if ( rDampingMatrix.size1() != 0 )
        rDampingMatrix.resize( 0, 0, false );

  KRATOS_CATCH( "" )
}

int LineHeatFluxCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

} // Namespace Kratos


