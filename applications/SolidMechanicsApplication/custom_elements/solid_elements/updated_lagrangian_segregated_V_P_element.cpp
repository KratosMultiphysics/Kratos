//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:           May 2018 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/updated_lagrangian_segregated_V_P_element.hpp"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement::UpdatedLagrangianSegregatedVPElement()
    : LargeDisplacementSegregatedVPElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement::UpdatedLagrangianSegregatedVPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementSegregatedVPElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement::UpdatedLagrangianSegregatedVPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementSegregatedVPElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement::UpdatedLagrangianSegregatedVPElement( UpdatedLagrangianSegregatedVPElement const& rOther)
    :LargeDisplacementSegregatedVPElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement&  UpdatedLagrangianSegregatedVPElement::operator=(UpdatedLagrangianSegregatedVPElement const& rOther)
{
  LargeDisplacementSegregatedVPElement::operator=(rOther);

  return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianSegregatedVPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< UpdatedLagrangianSegregatedVPElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianSegregatedVPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

  UpdatedLagrangianSegregatedVPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

  //-----------//

  NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


  if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
  {
    NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
      KRATOS_ERROR << "constitutive law not has the correct size " << NewElement.mConstitutiveLawVector.size() << std::endl;
  }

  for(unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
  {
    NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
  }

  //-----------//

  if ( NewElement.mDeformationGradientF0.size() != mDeformationGradientF0.size() )
    NewElement.mDeformationGradientF0.resize(mDeformationGradientF0.size());

  for(unsigned int i=0; i<mDeformationGradientF0.size(); i++)
  {
    NewElement.mDeformationGradientF0[i] = mDeformationGradientF0[i];
  }

  NewElement.mDeterminantF0 = mDeterminantF0;

  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Kratos::make_shared< UpdatedLagrangianSegregatedVPElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedVPElement::~UpdatedLagrangianSegregatedVPElement()
{
}


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                                        std::vector<double>& rValues,
                                                                        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();


    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
      {
	mDeterminantF0[PointNumber] = rValues[PointNumber];

	mConstitutiveLawVector[PointNumber]->SetValue(rVariable, rValues[PointNumber], rCurrentProcessInfo);
      }

  }
  else{
    LargeDisplacementSegregatedVPElement::SetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
  }


}

//**********************************GET DOUBLE VALUE**********************************
//************************************************************************************


void UpdatedLagrangianSegregatedVPElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                                        std::vector<double>& rValues,
                                                                        const ProcessInfo& rCurrentProcessInfo )
{

  if (rVariable == DETERMINANT_F){

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
      rValues.resize( integration_points_number );

    for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
    {
      rValues[PointNumber] = mDeterminantF0[PointNumber];
    }

  }
  else{

    LargeDisplacementSegregatedVPElement::GetValueOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

  }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::Initialize()
{
    KRATOS_TRY

    LargeDisplacementSegregatedVPElement::Initialize();

    SizeType integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    //Resize historical variables
    if ( mDeformationGradientF0.size() != integration_points_number )
        mDeformationGradientF0.resize( integration_points_number );

    if ( mDeterminantF0.size() != integration_points_number )
        mDeterminantF0.resize( integration_points_number, false );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
    {
        mDeterminantF0[PointNumber] = 1;
        mDeformationGradientF0[PointNumber] = identity_matrix<double> (dimension);
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::FinalizeStepVariables( ElementDataType & rVariables, const double& rPointNumber )
{
    //update internal (historical) variables
    mDeterminantF0[rPointNumber] = rVariables.detF * rVariables.detF0;
    noalias(mDeformationGradientF0[rPointNumber]) = prod(rVariables.F, rVariables.F0);
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SolidElement::InitializeElementData(rVariables,rCurrentProcessInfo);

        break;
      }
    case PRESSURE_STEP:
      {

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        //initialize element variables
        rVariables.N.resize(number_of_nodes,false);
        rVariables.H.resize(dimension,dimension,false);
        rVariables.F.resize(dimension,dimension,false);
        rVariables.DN_DX.resize(number_of_nodes, dimension,false);
        rVariables.DeltaPosition.resize(number_of_nodes, dimension,false);

        //reading shape functions
        rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

        //reading shape functions local gradients
        rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

        //set process info
        rVariables.SetProcessInfo(rCurrentProcessInfo);

        //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
        rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  //Calculate Delta Position
  rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

  //set variables including all integration points values

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateMaterialResponse(ElementDataType& rVariables,
                                                                     ConstitutiveLaw::Parameters& rValues,
                                                                     const int & rPointNumber)
{
    KRATOS_TRY

    switch( mStepVariable )
    {
      case VELOCITY_STEP:
        {
          SolidElement::CalculateMaterialResponse(rVariables,rValues,rPointNumber);
          break;
        }
      case PRESSURE_STEP:
        {
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

    KRATOS_CATCH( "" )
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);


    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Deformation Gradient F [dx_n+1/dx_n] to be updated
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ ); //overwrites detJ

    //Compute cartesian derivatives [dN/dx_n+1]
    noalias(rVariables.DN_DX) = prod( DN_De[rPointNumber], Invj ); //overwrites DX now is the current position dx

    //Determinant of the Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {

          //Parent to reference configuration
          rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

          //Deformation Gradient F0
          rVariables.F0    = mDeformationGradientF0[rPointNumber];

          //Compute the deformation matrix B
          const GeometryType& rGeometry = GetGeometry();
          ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B,rGeometry,rVariables.DN_DX);

          break;
        }
      case PRESSURE_STEP:
        {

          GeometryType&  rGeometry = GetGeometry();
          //Calculate velocity gradient matrix
          ElementUtilities::CalculateVelocityGradient( rVariables.H, rGeometry, rVariables.DN_DX );

          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }

    KRATOS_CATCH( "" )
}



//*********************************COMPUTE KINETICS***********************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateKinetics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //TotalDeltaPosition must not be used in this element as mDeterminantF0 and mDeformationGradientF0 are stored for reduced order
    //however then the storage of variables in the full integration order quadrature must be considered

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    rVariables.DeltaPosition = this->CalculateTotalDeltaPosition(rVariables.DeltaPosition);

    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod, rVariables.DeltaPosition );

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], InvJ, rVariables.detJ);

    //Calculating the cartesian derivatives [dN/dx_n] = [dN/d£][d£/dx_0]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

    //Deformation Gradient F [dx_n+1/dx_0] = [dx_n+1/d£] [d£/dx_0]
    noalias( rVariables.F ) = prod( rVariables.j[rPointNumber], InvJ );

    //Determinant of the deformation gradient F
    rVariables.detF  = MathUtils<double>::Det(rVariables.F);

    //Determinant of the Deformation Gradient F0
    // (in this element F = F0, then F0 is set to the identity for coherence in the constitutive law)
    rVariables.detF0 = 1;
    rVariables.F0    = identity_matrix<double> ( dimension );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                                                              ElementDataType& rVariables)

{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = rGeometry.PointsNumber();

  // operation performed: calculate stabilization factor
  this->CalculateStabilizationTau(rVariables);

  // Get Free surface Faces
  std::vector<std::vector<SizeType> > Faces;
  this->GetFreeSurfaceFaces(Faces);

  // Add Boundary Matrix
  if( Faces.size() != 0 ){ //if there are free surfaces

    double SideWeight = 0;
    double side_normal_size = 0;
    double BoundFactor = 0;

    for( SizeType i=0; i<Faces.size(); ++i ){

      GetFaceWeight(Faces[i], rVariables, SideWeight, side_normal_size);

      BoundFactor = rVariables.Tau * 2.0 / side_normal_size;

      //(lumped)
      // for( SizeType j=0; j<Faces[i].size(); ++j ){
      //   rLeftHandSideMatrix(Faces[i][j],Faces[i][j]) += rVariables.N[Faces[i][j]] * BoundFactor * SideWeight;
      // }
      //(reduced integration)
      for( SizeType j=0; j<Faces[i].size(); ++j ){
        for( SizeType k=0; k<Faces[i].size(); ++k ){
          rLeftHandSideMatrix(Faces[i][j],Faces[i][k]) += BoundFactor * SideWeight * rVariables.N[Faces[i][j]] * rVariables.N[Faces[i][k]];
        }
      }
    }

  }

  // Add Stabilized Laplacian Matrix
  double StabilizationFactor = rVariables.Tau * rVariables.IntegrationWeight;
  for( SizeType i=0; i<number_of_nodes; ++i )
  {
    for( SizeType j=0; j<number_of_nodes; ++j )
    {
      for ( SizeType k = 0; k<dimension; ++k )
      {
        rLeftHandSideMatrix(i,j) += StabilizationFactor * rVariables.DN_DX(i,k) * rVariables.DN_DX(j,k);
      }
    }
  }


  // Add Bulk Matrix
  const double& YoungModulus = GetProperties()[YOUNG_MODULUS];
  const double& Poisson = GetProperties()[POISSON_RATIO];
  double LameMu = YoungModulus/(2.0*(1.0+Poisson));
  double BulkModulus = (YoungModulus * Poisson)/((1.0+Poisson)*(1.0-2.0*Poisson)) + (2.0/3.0) * LameMu;

  if( GetProperties().Has(BULK_MODULUS) )
    BulkModulus = GetProperties()[BULK_MODULUS];

  const double& Density = GetProperties()[DENSITY];
  const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];

  double MassFactor = rVariables.IntegrationWeight / (BulkModulus * TimeStep);
  double BulkFactor = MassFactor * Density * rVariables.Tau / TimeStep;

  // (LUMPED)
  // double coefficient = rGeometry.IntegrationPointsNumber() * (1 + dimension); //integration points independent
  // MassFactor /= coefficient;
  // BulkFactor /= coefficient;
  // for( SizeType i=0; i<number_of_nodes; ++i)
  // {
  //   rLeftHandSideMatrix(i,i) += MassFactor + BulkFactor;
  // }

  // (REDUCED INTEGRATION)
  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    for( SizeType j=0; j<number_of_nodes; ++j)
    {
      rLeftHandSideMatrix(i,j) += (MassFactor + BulkFactor) * rVariables.N[i] * rVariables.N[j];
    }
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                                                         ElementDataType & rVariables)
{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = rGeometry.PointsNumber();

  // operation performed: calculate stabilization factor
  this->CalculateStabilizationTau(rVariables);

  // Get Free surface Faces
  std::vector<std::vector<SizeType> > Faces;
  this->GetFreeSurfaceFaces(Faces);

  // Add Boundary Vector
  if( Faces.size() != 0 ){ //if there are free surfaces

    Vector Normal(dimension);
    noalias(Normal) = ZeroVector(dimension);
    double ProjectionVelocityGradient = 0;
    double BoundFactor  = 0;
    double BoundFactorA = 0;
    double BoundFactorB = 0;
    double SideWeight = 0;
    double side_normal_size = 0;

    const double& YoungModulus = GetProperties()[YOUNG_MODULUS];
    const double& Poisson = GetProperties()[POISSON_RATIO];

    // Get element properties
    const double& Density   = GetProperties()[DENSITY];
    double LameMu = YoungModulus/(2.0*(1.0+Poisson));
    const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];

    //h_n (normal h)
    Matrix D(dimension,dimension);
    noalias(D) = 0.5 * (trans(rVariables.H)+rVariables.H);


    for( SizeType i=0; i<Faces.size(); ++i ){

      GetFaceNormal(Faces[i], rVariables, Normal);

      GetFaceWeight(Faces[i], rVariables, SideWeight, side_normal_size);

      Vector Acceleration (dimension);
      noalias(Acceleration) = ZeroVector(dimension);

      ProjectionVelocityGradient = inner_prod(Normal, prod(D,Normal));

      BoundFactor = rVariables.Tau * 2.0 / side_normal_size;

      BoundFactorA = rVariables.Tau * Density;
      BoundFactorB = rVariables.Tau * 4.0 * ProjectionVelocityGradient * (LameMu * TimeStep) / side_normal_size;

      // Vector NodeNormal (dimension);
      // noalias(NodeNormal) = ZeroVector(dimension);

      for( SizeType j=0; j<Faces[i].size(); ++j )
      {

        for ( SizeType k = 0; k<dimension; ++k )
        {
          Acceleration[k] = rGeometry[Faces[i][j]].FastGetSolutionStepValue(ACCELERATION)[k];
          //NodeNormal[k] = rGeometry[Faces[i][j]].FastGetSolutionStepValue(NORMAL)[k];
        }

        //rRightHandSideVector[Faces[i][j]] += SideWeight * rVariables.N[Faces[i][j]] * (BoundFactorA * inner_prod(Acceleration,NodeNormal) - BoundFactorB);
        rRightHandSideVector[Faces[i][j]] += SideWeight * rVariables.N[Faces[i][j]] * (BoundFactorA * inner_prod(Acceleration,Normal) - BoundFactorB);


        // Add LHS to RHS: boundary terms (incremental pressure formulation)
        //(lumped)
        //rRightHandSideVector[Faces[i][j]] -=  SideWeight * BoundFactor * rVariables.N[Faces[i][j]] * rGeometry[Faces[i][j]].FastGetSolutionStepValue(PRESSURE);

        //(reduced integration)
        for( SizeType k=0; k<Faces[i].size(); ++k ){
          rRightHandSideVector[Faces[i][j]] -= SideWeight * BoundFactor * rVariables.N[Faces[i][j]] * rVariables.N[Faces[i][k]] * rGeometry[Faces[i][k]].FastGetSolutionStepValue(PRESSURE);
        }
      }

    }

  }

  // Add Divergence and volume acceleration vector
  double TraceVelocityGradient = 0;
  for (SizeType i = 0; i < dimension; i++)
  {
    TraceVelocityGradient += rVariables.H(i,i);
  }

  Vector VolumeForce(dimension);
  noalias(VolumeForce) = ZeroVector(dimension);
  VolumeForce = this->CalculateVolumeForce( VolumeForce, rVariables );

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    // Velocity divergence
    rRightHandSideVector[i] += rVariables.IntegrationWeight * rVariables.N[i] * TraceVelocityGradient;

    // Volume forces
    for (SizeType j=0; j<dimension; ++j)
    {
      rRightHandSideVector[i] -= rVariables.Tau * rVariables.IntegrationWeight * rVariables.DN_DX(i,j) * VolumeForce[j];
    }
  }

  // Add Dynamic Bulk Vector
  const double& YoungModulus = GetProperties()[YOUNG_MODULUS];
  const double& Poisson = GetProperties()[POISSON_RATIO];
  double LameMu = YoungModulus/(2.0*(1.0+Poisson));
  double BulkModulus = (YoungModulus * Poisson)/((1.0+Poisson)*(1.0-2.0*Poisson)) + (2.0/3.0) * LameMu;
  if( GetProperties().Has(BULK_MODULUS) )
    BulkModulus = GetProperties()[BULK_MODULUS];

  const double& Density = GetProperties()[DENSITY];

  double MassFactor = rVariables.IntegrationWeight / BulkModulus;
  double BulkFactor = MassFactor * Density * rVariables.Tau;

  // (LUMPED)
  // double coefficient = rGeometry.IntegrationPointsNumber() * (1 + dimension); //integration points independent
  // MassFactor /= coefficient;
  // BulkFactor /= coefficient;
  // const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];
  // for( SizeType i=0; i<number_of_nodes; ++i)
  // {
  //   rRightHandSideVector[i] -= MassFactor * (1.0/TimeStep) * rGeometry[j].FastGetSolutionStepValue(PRESSURE_VELOCITY);
  //   rRightHandSideVector[i] -= BulkFactor * (1.0/(TimeStep*TimeStep)) * rGeometry[j].FastGetSolutionStepValue(PRESSURE_ACCELERATION);
  // }

  // (REDUCED INTEGRATION)
  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    for( SizeType j=0; j<number_of_nodes; ++j)
    {
      rRightHandSideVector[i] -= MassFactor * rVariables.N[i] * rVariables.N[j] * rGeometry[j].FastGetSolutionStepValue(PRESSURE_VELOCITY);
      rRightHandSideVector[i] -= BulkFactor * rVariables.N[i] * rVariables.N[j] * rGeometry[j].FastGetSolutionStepValue(PRESSURE_ACCELERATION);
    }
  }

  // Add LHS to RHS: stabilization terms (incremental pressure formulation)

  // Add Stabilized Laplacian Matrix to RHS
  double StabilizationFactor = rVariables.Tau * rVariables.IntegrationWeight;

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    for( SizeType j=0; j<number_of_nodes; ++j)
    {
      for ( SizeType k = 0; k < dimension; ++k )
      {
        rRightHandSideVector[i] -= (StabilizationFactor * rVariables.DN_DX(i,k) * rVariables.DN_DX(j,k)) * rGeometry[j].FastGetSolutionStepValue(PRESSURE);
      }
    }
  }


  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::CalculateStabilizationTau(ElementDataType& rVariables)

{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  SizeType number_of_nodes = rGeometry.PointsNumber();

  // Get mean velocity norm
  array_1d<double,3> MeanVelocity;
  noalias(MeanVelocity) = ZeroVector(3);
  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    MeanVelocity += rGeometry[i].FastGetSolutionStepValue(VELOCITY);
  }
  double mean_velocity = norm_2(MeanVelocity)/double(number_of_nodes);

  // Calculate FIC stabilization coefficient
  rVariables.Tau = 0;
  if( mean_velocity != 0 ){

    const double& YoungModulus = GetProperties()[YOUNG_MODULUS];
    const double& Poisson = GetProperties()[POISSON_RATIO];

    // Get element properties
    const double& Density   = GetProperties()[DENSITY];
    double LameMu = YoungModulus/(2.0*(1.0+Poisson));
    const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];

    // Get element size
    double element_size = rGeometry.AverageEdgeLength();

    rVariables.Tau = (element_size * element_size * TimeStep) / ( Density * mean_velocity * TimeStep * element_size + Density * element_size * element_size +  8.0 * LameMu * TimeStep * TimeStep );

  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::GetFreeSurfaceFaces(std::vector<std::vector<SizeType> >& Faces)
{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  DenseMatrix<unsigned int> NodesInFaces;
  rGeometry.NodesInFaces(NodesInFaces);

  //based in existance of neighbour elements (proper detection for triangles and tetrahedra)
  WeakPointerVector<Element>& neighb_elems = this->GetValue(NEIGHBOUR_ELEMENTS);
  unsigned int counter=0;
  for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
  {
    if (ne->Id() == this->Id())  // If there is no shared element in face nf (the Id coincides)
    {
      std::vector<SizeType> Nodes;
      unsigned int FixedNodes  = 0;

      for(unsigned int i = 1; i < NodesInFaces.size1(); i++)
      {
        Nodes.push_back(NodesInFaces(i,counter));  //set boundary nodes
        if(rGeometry[NodesInFaces(i,counter)].IsFixed(VELOCITY_X) || rGeometry[NodesInFaces(i,counter)].IsFixed(VELOCITY_Y) ){
          ++FixedNodes;
        }
        else if(dimension == 3){
          if(rGeometry[NodesInFaces(i,counter)].IsFixed(VELOCITY_Z)) {
            ++FixedNodes;
          }
        }

      }
      if( FixedNodes < Nodes.size() )
        Faces.push_back(Nodes);
    }

    counter++;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::GetFaceNormal(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, Vector& rNormal)
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  // for triangles and tetrahedra
  if( rNormal.size() != dimension )
    rNormal.resize(dimension,false);

  noalias(rNormal) = ZeroVector(dimension);
  for( SizeType j=0; j<rFace.size(); ++j )
  {
    for(unsigned int d=0; d<dimension; d++)
    {
      rNormal[d] += rVariables.DN_DX(rFace[j],d);
    }
  }

  double norm = norm_2(rNormal);
  if(norm!=0)
    rNormal /= norm;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::GetFaceWeight(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, double& rWeight, double& rNormalSize)
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  // for triangles and tetrahedra
  Vector An(dimension);
  noalias(An) = ZeroVector(dimension);
  for( SizeType j=0; j<rFace.size(); ++j )
  {
    for(unsigned int d=0; d<dimension; d++)
    {
      An[d] += rVariables.DN_DX(rFace[j],d);
    }
  }

  double norm = norm_2(An);
  rNormalSize = 1.0/norm;
  rWeight = dimension * rVariables.IntegrationWeight * norm;

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::GetFaceNormal(const std::vector<SizeType>& rFace, Vector& rNormal)
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  bool computed = false;
  if( dimension == 2 ){

    if( rNormal.size() != 2 )
      rNormal.resize(2,false);

    if( rFace.size() == 2 ) {
      rNormal[0] =    GetGeometry()[rFace[1]].Y() - GetGeometry()[rFace[0]].Y();
      rNormal[1] = -( GetGeometry()[rFace[1]].X() - GetGeometry()[rFace[0]].X());

      double norm = norm_2(rNormal);
      if( norm != 0 )
        rNormal /= norm_2(rNormal);

      computed = true;
    }

  }
  else if( dimension == 3 ){

    if( rNormal.size() != 3 )
      rNormal.resize(3,false);

    if( rFace.size() == 3 ) {

      Vector v1(3);
      Vector v2(3);
      v1[0] = GetGeometry()[rFace[1]].X() - GetGeometry()[rFace[0]].X();
      v1[1] = GetGeometry()[rFace[1]].Y() - GetGeometry()[rFace[0]].Y();
      v1[2] = GetGeometry()[rFace[1]].Z() - GetGeometry()[rFace[0]].Z();

      v2[0] = GetGeometry()[rFace[2]].X() - GetGeometry()[rFace[0]].X();
      v2[1] = GetGeometry()[rFace[2]].Y() - GetGeometry()[rFace[0]].Y();
      v2[2] = GetGeometry()[rFace[2]].Z() - GetGeometry()[rFace[0]].Z();

      MathUtils<double>::CrossProduct(rNormal,v1,v2);
      double norm = norm_2(rNormal);
      if( norm != 0 )
        rNormal /= norm_2(rNormal);

      computed = true;
    }
  }

  if( !computed ){

     if( rNormal.size() != dimension )
       rNormal.resize(dimension,false);

    noalias(rNormal) = ZeroVector(dimension);

    double coefficient = 1.0 / double(rFace.size());
    for( SizeType i = 0; i<rFace.size(); ++i )
    {
      for ( SizeType k = 0; k<dimension; ++k )
      {
        rNormal[k] += coefficient * GetGeometry()[rFace[i]].FastGetSolutionStepValue(NORMAL)[k];
        //here the normal of the boundary can be calculated (more precise) if normals not updated
      }
    }

    double norm = norm_2(rNormal);
    if( norm != 0 )
      rNormal /= norm_2(rNormal);
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedVPElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber )
{
    LargeDisplacementElement::GetHistoricalVariables(rVariables,rPointNumber);

    //Deformation Gradient F0
    rVariables.detF0 = mDeterminantF0[rPointNumber];
    rVariables.F0    = mDeformationGradientF0[rPointNumber];
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianSegregatedVPElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0);

    return rVolumeChange;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int  UpdatedLagrangianSegregatedVPElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // Perform base element checks
  int ErrorCode = 0;
  ErrorCode = LargeDisplacementSegregatedVPElement::Check(rCurrentProcessInfo);

  // Check compatibility with the constitutive law
  ConstitutiveLaw::Features LawFeatures;
  this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

  // Check that the constitutive law has the correct dimension
  SizeType dimension = this->GetGeometry().WorkingSpaceDimension();
  if( dimension == 2 )
  {
    if( LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRAIN_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::AXISYMMETRIC_LAW) )
      KRATOS_ERROR <<  "wrong constitutive law used. This is a 2D element. Expected plane state or axisymmetric :: element id = " << this->Id() << std::endl;
  }

  // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
  for(SizeType i=0; i<this->GetGeometry().size(); ++i)
  {
    // Nodal data
    Node<3> &rNode = this->GetGeometry()[i];
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE,rNode);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_VELOCITY,rNode);
    KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_ACCELERATION,rNode);
    //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUME_ACCELERATION,rNode);

    // Nodal dofs
    KRATOS_CHECK_DOF_IN_NODE(PRESSURE,rNode);
  }
  // Check compatibility with the constitutive law

  // Check that all required variables have been registered

  return ErrorCode;

  KRATOS_CATCH( "" );
}


//************************************************************************************
//************************************************************************************


void UpdatedLagrangianSegregatedVPElement::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementSegregatedVPElement )
  rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
  rSerializer.save("DeterminantF0",mDeterminantF0);
}

void UpdatedLagrangianSegregatedVPElement::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementSegregatedVPElement )
  rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
  rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos
