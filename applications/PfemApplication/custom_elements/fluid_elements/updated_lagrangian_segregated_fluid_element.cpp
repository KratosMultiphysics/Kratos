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
#include "custom_elements/fluid_elements/updated_lagrangian_segregated_fluid_element.hpp"
#include "pfem_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement()
    : FluidElement()
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : FluidElement( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : FluidElement( NewId, pGeometry, pProperties )
{
  mStepVariable = VELOCITY_STEP;
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement( UpdatedLagrangianSegregatedFluidElement const& rOther)
    :FluidElement(rOther)
    ,mStepVariable(rOther.mStepVariable)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement&  UpdatedLagrangianSegregatedFluidElement::operator=(UpdatedLagrangianSegregatedFluidElement const& rOther)
{
  FluidElement::operator=(rOther);

  mStepVariable = rOther.mStepVariable;

  return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianSegregatedFluidElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared< UpdatedLagrangianSegregatedFluidElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianSegregatedFluidElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

  UpdatedLagrangianSegregatedFluidElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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
  
  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());
  
  return Kratos::make_shared< UpdatedLagrangianSegregatedFluidElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::~UpdatedLagrangianSegregatedFluidElement()
{
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetHistoricalVariables( ElementDataType& rVariables, const double& rPointNumber)
{
  //Deformation Gradient F ( set to identity )
  SizeType size =  rVariables.F.size1();

  rVariables.Alpha = 1.0;
  rVariables.detF = 1;
  noalias(rVariables.F) = IdentityMatrix(size);

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rElementalDofList.resize( 0 );

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  switch( StepType(rCurrentProcessInfo[SEGREGATED_STEP]) )
  {
    case VELOCITY_STEP:
      {
        for ( SizeType i = 0; i < GetGeometry().size(); i++ )
        {
          rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
          rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

          if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
        }
        break;
      }
    case PRESSURE_STEP:
      {
        for ( SizeType i = 0; i < GetGeometry().size(); i++ )
        {
          rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
        }
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  unsigned int   dofs_size       = GetDofsSize();

  if ( rResult.size() != dofs_size )
    rResult.resize( dofs_size, false );
  
  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          int index = i * dimension;
          rResult[index]     = GetGeometry()[i].GetDof( VELOCITY_X ).EquationId();
          rResult[index + 1] = GetGeometry()[i].GetDof( VELOCITY_Y ).EquationId();

          if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof( VELOCITY_Z ).EquationId();
        }
        break;
      }
    case PRESSURE_STEP:
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          rResult[i] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  this->SetProcessInformation(rCurrentProcessInfo);
  
  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY


  KRATOS_CATCH( "" )
}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetValuesVector( Vector& rValues, int Step )
{
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  unsigned int   dofs_size       = GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SizeType index = 0;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          index = i * dimension;
          rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
          rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

          if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

        }
        break;
      }
    case PRESSURE_STEP:
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
        }
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

}


//************************************VELOCITY****************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  unsigned int   dofs_size       = GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SizeType index = 0;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          index = i * dimension;
          rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
          rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

          if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
        }
        break;
      }
    case PRESSURE_STEP:
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_VELOCITY, Step );
        }
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

}

//*********************************ACCELERATION***************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  unsigned int   dofs_size       = GetDofsSize();

  if ( rValues.size() != dofs_size )
    rValues.resize( dofs_size, false );

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        SizeType index = 0;
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          index = i * dimension;
          rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
          rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

          if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
        break;
      }
    case PRESSURE_STEP:
      {
        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
          rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_ACCELERATION, Step );
        }
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    FluidElement::InitializeElementData(rVariables,rCurrentProcessInfo);

    //Calculate Delta Position
    rVariables.DeltaPosition = this->CalculateDeltaPosition(rVariables.DeltaPosition);

    //set variables including all integration points values

    //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetStepAlpha(double& rAlpha)
{
    KRATOS_TRY

    switch( mStepVariable )
    {
      case VELOCITY_STEP:
        {
          rAlpha = 0.5;
          break;
        }
      case PRESSURE_STEP:
        {
          rAlpha = 1.0;
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }
    
    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //Get integration point Alpha parameter
    GetStepAlpha(rVariables.Alpha);

    //to take in account previous step for output print purposes
    if( this->Is(FluidElement::FINALIZED_STEP) ){
      this->GetHistoricalVariables(rVariables,rPointNumber);
    }
    
    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Parent to reference configuration
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives [dN/dx_n]
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber], InvJ );

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
    //
    //rVariables.detF0 = mDeterminantF0[rPointNumber];
    //rVariables.F0    = mDeformationGradientF0[rPointNumber];

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);
   
    //Calculate velocity gradient matrix
    CalculateVelocityGradient( rVariables.L, rVariables.DN_DX, rVariables.Alpha );

    
    // const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    // for (SizeType k = 0; k < number_of_nodes; k++)
    // {
    //   std::cout<<" Velocity["<<k<<"] "<<GetGeometry()[k].FastGetSolutionStepValue(VELOCITY)<<std::endl;
    // }
    
    KRATOS_CATCH( "" )
}

//*********************************COMPUTE KINETICS***********************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateKinetics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    //TotalDeltaPosition must not be used in this element as mDeterminantF0 and mDeformationGradientF0 are stored for reduced order
    //however then the storage of variables in the full integration order quadrature must be considered

    //const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

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
    //rVariables.detF0 = 1;
    //rVariables.F0    = identity_matrix<double> ( dimension );

    //Set Shape Functions Values for this integration point
    noalias(rVariables.N) = matrix_row<const Matrix>( Ncontainer, rPointNumber);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateDeformationMatrix(Matrix& rB,
                                                                         const Matrix& rF,
                                                                         const Matrix& rDN_DX)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    rB.clear(); //set all components to zero

    if( dimension == 2 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        for ( SizeType i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 3 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 2 ) = rDN_DX( i, 2 );

            rB( 3, index + 0 ) = rDN_DX( i, 1 );
            rB( 3, index + 1 ) = rDN_DX( i, 0 );

            rB( 4, index + 1 ) = rDN_DX( i, 2 );
            rB( 4, index + 2 ) = rDN_DX( i, 1 );

            rB( 5, index + 0 ) = rDN_DX( i, 2 );
            rB( 5, index + 2 ) = rDN_DX( i, 0 );

        }
    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing B " << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateVelocityGradient(Matrix& rL,
                                                                        const Matrix& rDN_DX,
                                                                        double Alpha)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    rL = ZeroMatrix(dimension,dimension);

    array_1d<double,3> Velocity; 
    if( dimension == 2 )
    {

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        array_1d<double,3>& rPreviousVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double,3>& rCurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
          
        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);
          
        rL ( 0 , 0 ) += Velocity[0]*rDN_DX ( i , 0 );
        rL ( 0 , 1 ) += Velocity[0]*rDN_DX ( i , 1 );
        rL ( 1 , 0 ) += Velocity[1]*rDN_DX ( i , 0 );
        rL ( 1 , 1 ) += Velocity[1]*rDN_DX ( i , 1 );
      }

    }
    else if( dimension == 3)
    {

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {

        array_1d<double,3>& rPreviousVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double,3>& rCurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
          
        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);

        rL ( 0 , 0 ) += Velocity[0]*rDN_DX ( i , 0 );
        rL ( 0 , 1 ) += Velocity[0]*rDN_DX ( i , 1 );
        rL ( 0 , 2 ) += Velocity[0]*rDN_DX ( i , 2 );
        rL ( 1 , 0 ) += Velocity[1]*rDN_DX ( i , 0 );
        rL ( 1 , 1 ) += Velocity[1]*rDN_DX ( i , 1 );
        rL ( 1 , 2 ) += Velocity[1]*rDN_DX ( i , 2 );
        rL ( 2 , 0 ) += Velocity[2]*rDN_DX ( i , 0 );
        rL ( 2 , 1 ) += Velocity[2]*rDN_DX ( i , 1 );
        rL ( 2 , 2 ) += Velocity[2]*rDN_DX ( i , 2 );
      }

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateVelocityGradientVector(Vector& rL,
                                                                              const Matrix& rDN_DX,
                                                                              double Alpha)
{
    KRATOS_TRY

    const SizeType number_of_nodes  = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    rL = ZeroVector(dimension * dimension );
    array_1d<double,3> Velocity; 
    if( dimension == 2 )
    {
      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        array_1d<double,3>& rPreviousVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double,3>& rCurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
          
        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);

        rL[0] += Velocity[0]*rDN_DX ( i , 0 );
        rL[1] += Velocity[1]*rDN_DX ( i , 1 );
        rL[2] += Velocity[0]*rDN_DX ( i , 1 );
        rL[3] += Velocity[1]*rDN_DX ( i , 0 );
      }
    }
    else if( dimension == 3)
    {

      for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
        array_1d<double,3>& rPreviousVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,1);
        array_1d<double,3>& rCurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
          
        Velocity = rCurrentVelocity * Alpha + rPreviousVelocity * (1.0-Alpha);

        rL[0] += Velocity[0]*rDN_DX ( i , 0 );
        rL[1] += Velocity[1]*rDN_DX ( i , 1 );
        rL[2] += Velocity[2]*rDN_DX ( i , 2 );

        rL[3] += Velocity[0]*rDN_DX ( i , 1 );
        rL[4] += Velocity[1]*rDN_DX ( i , 2 );
        rL[5] += Velocity[2]*rDN_DX ( i , 0 );

        rL[6] += Velocity[1]*rDN_DX ( i , 0 );
        rL[7] += Velocity[2]*rDN_DX ( i , 1 );
        rL[8] += Velocity[0]*rDN_DX ( i , 2 );
      }

    }
    else
    {
      KRATOS_ERROR << " something is wrong with the dimension when computing velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateSymmetricVelocityGradient(const Matrix& rL,
                                                                                 Vector& rStrainVector)
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = rL( 0, 0 );
        rStrainVector[1] = rL( 1, 1 );
        rStrainVector[2] = 0.5 * (rL( 0, 1 ) + rL( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = rL( 0, 0 );
        rStrainVector[1] = rL( 1, 1 );
        rStrainVector[2] = rL( 2, 2 );
        rStrainVector[3] = 0.5 * ( rL( 0, 1 ) + rL( 1, 0 ) ); // xy
        rStrainVector[4] = 0.5 * ( rL( 1, 2 ) + rL( 2, 1 ) ); // yz
        rStrainVector[5] = 0.5 * ( rL( 0, 2 ) + rL( 2, 0 ) ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension symmetric velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateSkewSymmetricVelocityGradient(const Matrix& rL,
                                                                                     Vector& rStrainVector)
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = 0.0;
        rStrainVector[1] = 0.0;
        rStrainVector[2] = 0.5 * (rL( 0, 1 ) - rL( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = 0.0;
        rStrainVector[1] = 0.0;
        rStrainVector[2] = 0.0;
        rStrainVector[3] = 0.5 * ( rL( 0, 1 ) - rL( 1, 0 ) ); // xy
        rStrainVector[4] = 0.5 * ( rL( 1, 2 ) - rL( 2, 1 ) ); // yz
        rStrainVector[5] = 0.5 * ( rL( 0, 2 ) - rL( 2, 0 ) ); // xz

    }
    else
    {
        KRATOS_ERROR << " something is wrong with the dimension symmetric velocity gradient " << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::SetElementData(ElementDataType& rVariables,
                                                             ConstitutiveLaw::Parameters& rValues,
                                                             const int & rPointNumber)
{

  if(rVariables.detF<0){

    std::cout<<" Element: "<<this->Id()<<std::endl;

    SizeType number_of_nodes = GetGeometry().PointsNumber();

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
      array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
      array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
      array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
      std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
      std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
    }

    KRATOS_ERROR << " UPDATED LAGRANGIAN SEGREGATED FLUID ELEMENT INVERTED: |F|<0  detF = " << rVariables.detF << std::endl;
  }

  //Compute strain rate measures if they are required by the constitutive law
  ConstitutiveLaw::Features LawFeatures;
  mConstitutiveLawVector[rPointNumber]->GetLawFeatures(LawFeatures);

  //Compute symmetric spatial velocity gradient [DN_DX = dN/dx_n*1] stored in a vector
  this->CalculateSymmetricVelocityGradient( rVariables.L, rVariables.StrainVector );

  //std::cout<<" SpatialDefRate "<<rVariables.StrainVector<<std::endl;
  
  Flags& ConstitutiveLawOptions = rValues.GetOptions();
  ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

  rValues.SetDeterminantF(rVariables.detH);
  rValues.SetDeformationGradientF(rVariables.F);
  rValues.SetStrainVector(rVariables.StrainVector);
  rValues.SetStressVector(rVariables.StressVector);
  rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
  rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
  rValues.SetShapeFunctionsValues(rVariables.N);

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables)
{
  KRATOS_TRY

  MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        // std::cout<<"V ELEMENT["<<this->Id()<<"]"<<std::endl;
        // operation performed: add stiffness term to the rLefsHandSideMatrix
        this->CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );

        //std::cout<<" L "<<rVariables.L<<" StressVector "<<rVariables.StressVector<<std::endl;
        
        // operation performed: add Kg to the rLefsHandSideMatrix
        // this->CalculateAndAddKvvg( rLeftHandSideMatrix, rVariables );

        //std::cout<< " velocity LHS "<<  rLeftHandSideMatrix << "(" << this->Id() << ")" <<std::endl;
        
        break;
      }
    case PRESSURE_STEP:
      {
        //std::cout<<"P ELEMENT["<<this->Id()<<"]"<<std::endl;

        // operation performed: calculate stabilization factor
        this->CalculateStabilizationTau(rVariables);

        // operation performed: add Kpp to the rLefsHandSideMatrix
        this->CalculateAndAddKpp( rLeftHandSideMatrix, rVariables );

        //std::cout<< " pressure LHS "<<  rLeftHandSideMatrix << "(" << this->Id() << ")" <<std::endl;

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  //KRATOS_WATCH( rLeftHandSideMatrix )

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables)
{
  KRATOS_TRY

  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        // operation performed: add InternalForces to the rRightHandSideVector
        this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables );

        // operation performed: add ExternalForces to the rRightHandSideVector
        this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables );

        //std::cout<< " velocity RHS "<<  rRightHandSideVector << "(" << this->Id() << ")" <<std::endl;

        break;
      }
    case PRESSURE_STEP:
      {
        // operation performed: calculate stabilization factor
        this->CalculateStabilizationTau(rVariables);

        // operation performed: add PressureForces to the rRightHandSideVector
        this->CalculateAndAddPressureForces( rRightHandSideVector, rVariables );

        //std::cout<< " pressure RHS "<<  rRightHandSideVector << "(" << this->Id() << ")" <<std::endl;

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  //KRATOS_WATCH( rRightHandSideVector )

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
                                                                            ElementDataType & rVariables)
{
    KRATOS_TRY

    //Add volumetric term to stress
    const SizeType number_of_nodes = this->GetGeometry().PointsNumber();
    double MeanPressure = 0;
    for( SizeType i=0; i<number_of_nodes; ++i)
    {
      // Pressure
      MeanPressure += rVariables.N[i] * this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
    }
    this->AddVolumetricPart(rVariables.StressVector, MeanPressure);

    //std::cout<<" MeanPressure "<<MeanPressure<<std::endl;
    
    FluidElement::CalculateAndAddInternalForces(rRightHandSideVector, rVariables);

    this->RemoveVolumetricPart(rVariables.StressVector, MeanPressure);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables)
{
  KRATOS_TRY

  if( StepType(rVariables.GetProcessInfo()[SEGREGATED_STEP]) == VELOCITY_STEP )
  {
    FluidElement::CalculateAndAddDynamicLHS(rLeftHandSideMatrix, rVariables);
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables)
{
  KRATOS_TRY

  if( StepType(rVariables.GetProcessInfo()[SEGREGATED_STEP]) == VELOCITY_STEP )
  {
    FluidElement::CalculateAndAddDynamicRHS(rRightHandSideVector, rVariables);
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // the internal step variable must be set because InitializeNonLinearIteration is not called before this method
  this->SetProcessInformation(rCurrentProcessInfo);
  
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        FluidElement::CalculateMassMatrix( rMassMatrix, rCurrentProcessInfo );
        break;
      }
    case PRESSURE_STEP:
      {
        const unsigned int MatSize = this->GetDofsSize();
        if ( rMassMatrix.size1() != MatSize )
          rMassMatrix.resize( MatSize, MatSize, false );

        noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  } 

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // the internal step variable must be set because InitializeNonLinearIteration is not called before this method
  this->SetProcessInformation(rCurrentProcessInfo);
  
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        FluidElement::CalculateDampingMatrix( rDampingMatrix, rCurrentProcessInfo );
        break;
      }
    case PRESSURE_STEP:
      {
        const unsigned int MatSize = this->GetDofsSize();
        if ( rDampingMatrix.size1() != MatSize )
          rDampingMatrix.resize( MatSize, MatSize, false );

        noalias(rDampingMatrix) = ZeroMatrix( MatSize, MatSize );

        // Add Bulk Damping Matrix
        // const double& BulkModulus = GetProperties()[BULK_MODULUS];
        // const double& Density     = GetProperties()[DENSITY];
        // const double& TimeStep    = rCurrentProcessInfo[DELTA_TIME];

        // double DampingFactor = GetGeometry().DomainSize() * Density * rVariables.Tau / (BulkModulus * TimeStep * TimeStep);

        // const SizeType number_of_nodes = GetGeometry().size();
        // for( SizeType i=0; i<number_of_nodes; ++i)
        // {
        //   rDampingMatrix(i,i) = DampingFactor;
        // }
        
        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
  } 

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

unsigned int UpdatedLagrangianSegregatedFluidElement::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = GetGeometry().PointsNumber();

  SizeType size = 0;
  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      size = number_of_nodes * dimension; //size for velocity
      break;
    case PRESSURE_STEP:
      size = number_of_nodes; //size for pressure
      break;
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;

  }
  return size;

  KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::SetProcessInformation(const ProcessInfo& rCurrentProcessInfo)
{
  KRATOS_TRY

  mStepVariable = StepType(rCurrentProcessInfo[SEGREGATED_STEP]);
      
  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKvvm(MatrixType& rLeftHandSideMatrix,
                                                                  ElementDataType& rVariables)
{
  KRATOS_TRY

  // Add volumetric part to the constitutive tensor to compute Kvvm = Kdev + Kvol(quasi incompressible)

  // add fluid bulk modulus for quasi-incompressibility (must be implemented in a ConstituiveModel for Quasi-Incompressible Newtonian Fluid.
  double BulkFactor = GetProperties()[BULK_MODULUS] * rVariables.TimeStep;

  // add fluid bulk modulus
  this->AddVolumetricPart(rVariables.ConstitutiveMatrix,BulkFactor);

  FluidElement::CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );

  rLeftHandSideMatrix *= rVariables.Alpha;

  //remove fluid bulk modulus
  this->RemoveVolumetricPart(rVariables.ConstitutiveMatrix,BulkFactor);

  // operation performed: set regularized stiffness term matrix to rLeftHandSideMatrix
  this->CalculateRegularizedKvvm( rLeftHandSideMatrix, rVariables );


  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateRegularizedKvvm(MatrixType& rLeftHandSideMatrix,
                                                                       ElementDataType& rVariables)
{
  KRATOS_TRY

  SizeType number_of_nodes = GetGeometry().PointsNumber();

  double StiffnessFactor = 0.0;
  this->CalculateDenseMatrixMeanValue(rLeftHandSideMatrix,StiffnessFactor);

  // Alternative calculation:
  // this->CalculateLumpedMatrixMeanValue(rLeftHandSideMatrix,StiffnessFactor);

  double MassFactor = GetGeometry().DomainSize() * GetProperties()[DENSITY] / double(number_of_nodes);

  double BulkFactor = 0.0;
  if(StiffnessFactor!=0 && MassFactor!=0)
    BulkFactor = GetProperties()[BULK_MODULUS] * rVariables.TimeStep * ( MassFactor * 2.0 / (rVariables.TimeStep * StiffnessFactor) );

  this->AddVolumetricPart(rVariables.ConstitutiveMatrix, BulkFactor);

  rLeftHandSideMatrix.clear();
  FluidElement::CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );
  rLeftHandSideMatrix *= rVariables.Alpha;

  this->RemoveVolumetricPart(rVariables.ConstitutiveMatrix, BulkFactor);

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKvvg(MatrixType& rLeftHandSideMatrix,
                                                                  ElementDataType& rVariables)

{
  KRATOS_TRY

  SizeType dimension = GetGeometry().WorkingSpaceDimension();
  Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
  Matrix ReducedKg = prod( rVariables.DN_DX, rVariables.IntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
  MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::AddVolumetricPart(Matrix& rConstitutiveMatrix, double& rBulkFactor)
{
  KRATOS_TRY

  //add fluid bulk modulus
  SizeType dimension = GetGeometry().WorkingSpaceDimension();
  for(SizeType i = 0; i < dimension; ++i)
  {
    for(SizeType j = 0; j < dimension; ++j)
    {
      rConstitutiveMatrix(i,j) += rBulkFactor;
    }
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::RemoveVolumetricPart(Matrix& rConstitutiveMatrix, double& rBulkFactor)
{
  KRATOS_TRY

  rBulkFactor *= (-1);
  this->AddVolumetricPart(rConstitutiveMatrix,rBulkFactor);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::AddVolumetricPart(Vector& rStressVector, double& rMeanPressure)
{
  KRATOS_TRY

  //add fluid bulk modulus
  SizeType dimension = GetGeometry().WorkingSpaceDimension();
  for(SizeType i = 0; i < dimension; ++i)
  {
    rStressVector[i] += rMeanPressure;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::RemoveVolumetricPart(Vector& rStressVector, double& rMeanPressure)
{
  KRATOS_TRY

  rMeanPressure *= (-1);
  this->AddVolumetricPart(rStressVector,rMeanPressure);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateDenseMatrixMeanValue(MatrixType& rMatrix, double& rMeanValue)
{
  KRATOS_TRY

  //Mean of all components of the matrix
  SizeType size1 = rMatrix.size1();
  SizeType size2 = rMatrix.size2();

  for(SizeType i = 0; i < size1; ++i)
  {
    for (SizeType j = 0; j < size2; ++j)
    {
      rMeanValue += fabs(rMatrix(i,j));
    }
  }

  rMeanValue /= double(size1*size2);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateLumpedMatrixMeanValue(MatrixType& rMatrix, double& rMeanValue)
{
  KRATOS_TRY

  //Mean of all components of the matrix
  SizeType size1 = rMatrix.size1();

  for(SizeType i = 0; i < size1; ++i)
  {
    rMeanValue += fabs(rMatrix(i,i));
  }

  rMeanValue /= double(size1);

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                                                                 ElementDataType& rVariables)

{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = rGeometry.PointsNumber();
  double element_size = rGeometry.AverageEdgeLength();

  // Get Free surface Faces
  std::vector<std::vector<SizeType> > Faces;
  this->GetFreeSurfaceFaces(Faces);

  // Add Boundary Matrix
  if( Faces.size() != 0 ){ //if there are free surfaces
    double BoundFactor = rVariables.Tau*4.0*rVariables.IntegrationWeight/(element_size*element_size);
    for( SizeType i=0; i<Faces.size(); ++i ){
      double coefficient = 1.0 / double(Faces[i].size()+1);
      for( SizeType j=0; j<Faces[i].size(); ++j ){
        rLeftHandSideMatrix(Faces[i][j],Faces[i][j]) += coefficient * BoundFactor;
      }
      //std::cout<<" BoundFactor: "<<BoundFactor<<" Kpp "<<rLeftHandSideMatrix<<std::endl;
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
  //std::cout<<" StabFactor: "<<StabilizationFactor<<" Kpp "<<rLeftHandSideMatrix<<std::endl;

  // Add Bulk Matrix
  const double& BulkModulus = GetProperties()[BULK_MODULUS];
  const double& Density     = GetProperties()[DENSITY];

  double coefficient = rGeometry.IntegrationPointsNumber(); //integration points independent

  double MassFactor = GetGeometry().DomainSize() / (BulkModulus*rVariables.TimeStep);
  MassFactor /= coefficient;

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    rLeftHandSideMatrix(i,i) += MassFactor;
  }

  //std::cout<<" MassFactor: "<<MassFactor<<" Kpp "<<rLeftHandSideMatrix<<std::endl;
  // Add Bulk Damping Matrix
  double BulkFactor = MassFactor * Density * rVariables.Tau / (rVariables.TimeStep);
  BulkFactor /= coefficient;

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    rLeftHandSideMatrix(i,i) += BulkFactor;
  }

  //std::cout<<" BulkFactor: "<<BulkFactor<<" Kpp "<<rLeftHandSideMatrix<<std::endl;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                                                            ElementDataType & rVariables)
{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = rGeometry.PointsNumber();
  double element_size = rGeometry.AverageEdgeLength();

  // Get Free surface Faces
  std::vector<std::vector<SizeType> > Faces;
  this->GetFreeSurfaceFaces(Faces);

  // Add Boundary Vector
  if( Faces.size() != 0 ){ //if there are free surfaces

    Vector Normal(dimension);
    noalias(Normal) = ZeroVector(dimension);
    double ProjectionVelocityGradient = 0;
    double BoundFactorA = 0; 
    double BoundFactorB = 0;
    const double& Viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double& Density   = GetProperties()[DENSITY];

    BoundFactorA = rVariables.Tau * Density * 2.0 * rVariables.IntegrationWeight / element_size;
    
    Matrix D = 0.5 * (trans(rVariables.L)+rVariables.L);
    
    for( SizeType i=0; i<Faces.size(); ++i ){

      //std::cout<<" Face "<<i<<std::endl;
      double coefficient = 1.0 / double(Faces[i].size());
      Normal.clear();
      for( SizeType j=0; j<Faces[i].size(); ++j ){
        for ( SizeType k = 0; k<dimension; ++k )
        {
          Normal[k] += coefficient * rGeometry[Faces[i][j]].FastGetSolutionStepValue(NORMAL)[k];
          //here the normal of the boundary can be calculated (more precise) if normals not updated
        }
      }
      
      // double norm = norm_2(Normal);
      // if( norm != 0 )
      //   Normal /= norm;
          
      //std::cout<<" Normal "<<Normal<<std::endl;
      //std::cout<<" DefRate "<<D<<std::endl;

      Vector Acceleration (dimension);
      noalias(Acceleration) = ZeroVector(dimension);

      Vector NodeNormal (dimension);
      noalias(NodeNormal) = ZeroVector(dimension);
      
      ProjectionVelocityGradient = inner_prod(Normal, prod(D,Normal));

      BoundFactorB = rVariables.Tau * 8.0 * ProjectionVelocityGradient * Viscosity * rVariables.IntegrationWeight / (element_size*element_size);

      double factor = 1.0/double(number_of_nodes);
      for( SizeType j=0; j<Faces[i].size(); ++j )
      {
        
        for ( SizeType k = 0; k<dimension; ++k )
        {
          //Acceleration[k] = rGeometry[Faces[i][j]].FastGetSolutionStepValue(ACCELERATION)[k];
          Acceleration[k] = 0.5 * ( rGeometry[Faces[i][j]].FastGetSolutionStepValue(VELOCITY)[k]-rGeometry[Faces[i][j]].FastGetSolutionStepValue(VELOCITY,1)[k] ) / rVariables.TimeStep - rGeometry[Faces[i][j]].FastGetSolutionStepValue(ACCELERATION,1)[k];
          NodeNormal[k] = rGeometry[Faces[i][j]].FastGetSolutionStepValue(NORMAL)[k];
        }

        // double norm = norm_2(Normal);
        // if( norm != 0 )
        //   Normal /= norm;
          
        //std::cout<<" Normal "<<Normal<<std::endl;
        //std::cout<<" Acceleration "<<Acceleration<<std::endl;
        rRightHandSideVector[Faces[i][j]] += factor * (BoundFactorA * inner_prod(Acceleration,NodeNormal) + BoundFactorB);
      }
    }
    //std::cout<<" NprojSpatial "<<ProjectionVelocityGradient<<" BoundRHSAcc "<<BoundFactorA<<std::endl;

    //std::cout<<" BoundRHSCoeffDev: "<<BoundFactorB<<" RHS "<<rRightHandSideVector<<std::endl;
  }

  //std::cout<<" RHS not clean "<<rRightHandSideVector<<std::endl;
  //rRightHandSideVector.clear();
  
  // Add Stabilized Laplacian Vector
  double StabilizationFactor = rVariables.Tau * rVariables.IntegrationWeight;

  double TraceVelocityGradient = 0;
  for (SizeType i = 0; i < dimension; i++)
  {
    TraceVelocityGradient += rVariables.L(i,i);
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
      rRightHandSideVector[i] -= rVariables.IntegrationWeight * rVariables.DN_DX(i,j) * rVariables.Tau * VolumeForce[j];
    }
  }

  //std::cout<<" Stab RHS "<<rRightHandSideVector<<" DefRate "<<TraceVelocityGradient<<std::endl;

  // Add Dynamic Bulk Vector
  const double& BulkModulus = GetProperties()[BULK_MODULUS];
  const double& Density     = GetProperties()[DENSITY];

  double coefficient = rGeometry.IntegrationPointsNumber(); //integration points independent

  double MassFactor = GetGeometry().DomainSize() / (BulkModulus*rVariables.TimeStep);
  MassFactor /= coefficient;

  double BulkFactor = MassFactor * Density * rVariables.Tau / (rVariables.TimeStep);
  BulkFactor /= coefficient;

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    rRightHandSideVector[i] -= MassFactor * (rGeometry[i].FastGetSolutionStepValue(PRESSURE,0) - rGeometry[i].FastGetSolutionStepValue(PRESSURE,1));
    rRightHandSideVector[i] -= BulkFactor * (rGeometry[i].FastGetSolutionStepValue(PRESSURE,0) - rGeometry[i].FastGetSolutionStepValue(PRESSURE,1) - rVariables.TimeStep * rGeometry[i].FastGetSolutionStepValue(PRESSURE_VELOCITY));
  }

  // Add LHS to RHS: stabilization and boundary terms

  // Add Boundary Matrix to RHS
  if( Faces.size() != 0 ){ //if there are free surfaces
    double BoundFactor = rVariables.Tau*4.0*rVariables.IntegrationWeight/(element_size*element_size);
    for( SizeType i=0; i<Faces.size(); ++i ){
      double coefficient = 1.0 / double(number_of_nodes); //integration points independent
      for( SizeType j=0; j<Faces[i].size(); ++j ){
        rRightHandSideVector[Faces[i][j]] -= coefficient * BoundFactor * rGeometry[Faces[i][j]].FastGetSolutionStepValue(PRESSURE,0);
      }
    }
  }

  // Add Stabilized Laplacian Matrix to RHS
  StabilizationFactor = rVariables.Tau * rVariables.IntegrationWeight;

  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    for( SizeType j=0; j<number_of_nodes; ++j)
    {
      for ( SizeType k = 0; k < dimension; ++k )
      {
        rRightHandSideVector[i] -= (StabilizationFactor * rVariables.DN_DX(i,k) * rVariables.DN_DX(j,k)) * rGeometry[j].FastGetSolutionStepValue(PRESSURE,0);
      }
    }
  }
  //std::cout<<" End RHS "<<rRightHandSideVector<<std::endl;
  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateStabilizationTau(ElementDataType& rVariables)

{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();
  SizeType number_of_nodes = rGeometry.PointsNumber();

  // Get mean velocity
  double mean_velocity = 0;
  for( SizeType i=0; i<number_of_nodes; ++i)
  {
    const array_1d<double,3>& Velocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
    mean_velocity += norm_2(Velocity);
  }
  mean_velocity /= double(number_of_nodes);

  // Calculate FIC stabilization coefficient
  rVariables.Tau = 0;
  if( mean_velocity != 0 ){

    // Get element properties
    const double& Density   = GetProperties()[DENSITY];
    const double& Viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    // Get element size
    double element_size = rGeometry.AverageEdgeLength();

    rVariables.Tau = (element_size * element_size * rVariables.TimeStep) / ( Density * mean_velocity * rVariables.TimeStep * element_size + Density * element_size * element_size +  8.0 * Viscosity * rVariables.TimeStep );

    //std::cout<<" Tau: "<<rVariables.Tau<<"(ElemSize:"<<element_size<<" Viscosity:"<<Viscosity<<")"<<std::endl;

  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetFreeSurfaceFaces(std::vector<std::vector<SizeType> >& Faces)
{
  KRATOS_TRY

  GeometryType& rGeometry = GetGeometry();

  DenseMatrix<unsigned int> NodesInFaces;
  rGeometry.NodesInFaces(NodesInFaces);

  //based on node flags (fail in edge elements)
  // for( SizeType i=0; i<NodesInFaces.size2(); ++i ){
  //   bool free_surface = true;
  //   for( SizeType j=1; j<NodesInFaces.size1(); ++j ){
  //     if( rGeometry[NodesInFaces(j,i)].IsNot(FREE_SURFACE) ){
  //       free_surface = false;
  //       break;
  //     }
  //   }
  //   if( free_surface ){
  //     std::vector<SizeType> Nodes;
  //     for( SizeType j=1; j<NodesInFaces.size1(); ++j ){
  //       if( rGeometry[NodesInFaces(j,i)].IsNot(INLET) ){
  //         Nodes.push_back(NodesInFaces(j,i));
  //       }
  //     }
  //     Faces.push_back(Nodes);
  //   }
  // }

  //based in existance of neighbour elements (proper detection for triangles and tetrahedra)
  WeakPointerVector<Element>& neighb_elems = this->GetValue(NEIGHBOUR_ELEMENTS);
  unsigned int counter=0;
  for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
  {
    if (ne->Id() == this->Id())  // If there is no shared element in face nf (the Id coincides)
    {
      std::vector<SizeType> Nodes;
      unsigned int WallNodes = 0;
      for(unsigned int i = 0; i < rGeometry.FacesNumber(); i++)
      {
        if(i!=counter){
          Nodes.push_back(NodesInFaces(i,0));  //set boundary nodes
          if(rGeometry[NodesInFaces(i,0)].Is(RIGID) || rGeometry[NodesInFaces(i,0)].Is(SOLID)){
            ++WallNodes;
          }
        }
      }
      if( WallNodes < Nodes.size() )
        Faces.push_back(Nodes);
    }

    counter++;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

int  UpdatedLagrangianSegregatedFluidElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
  KRATOS_TRY

  // Perform base element checks
  int ErrorCode = 0;
  ErrorCode = FluidElement::Check(rCurrentProcessInfo);

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
    //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_VELOCITY,rNode);
    //KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE_ACCELERATION,rNode);
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


void UpdatedLagrangianSegregatedFluidElement::save( Serializer& rSerializer ) const
{
  KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, FluidElement )
  int IntStepType = int(mStepVariable);
  rSerializer.save("StepVariable",IntStepType);
}

void UpdatedLagrangianSegregatedFluidElement::load( Serializer& rSerializer )
{
  KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidElement )
  int IntStepType;
  rSerializer.load("StepVariable",IntStepType);
  mStepVariable = StepType(IntStepType);
}


} // Namespace Kratos
