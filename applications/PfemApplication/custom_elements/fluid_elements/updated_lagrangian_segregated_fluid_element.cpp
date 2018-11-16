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
#include "custom_utilities/element_utilities.hpp"
#include "pfem_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement()
    : FluidElement()
{
  mStepVariable = VELOCITY_STEP;
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : FluidElement( NewId, pGeometry )
{
  mStepVariable = VELOCITY_STEP;
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

  NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


  if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
  {
    NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

    if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
      KRATOS_ERROR << "constitutive law not has the correct size " << NewElement.mConstitutiveLawVector.size() << std::endl;
  }

  for(unsigned int i=0; i<mConstitutiveLawVector.size(); ++i)
  {
    NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
  }

  //commented: it raises a segmentation fault in make_shared bad alloc
  NewElement.mStepVariable = mStepVariable;

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

void UpdatedLagrangianSegregatedFluidElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
  rElementalDofList.resize( 0 );

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
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
  this->SetProcessInformation(rCurrentProcessInfo);

  const SizeType number_of_nodes = GetGeometry().size();
  const SizeType dimension       = GetGeometry().WorkingSpaceDimension();
  unsigned int   dofs_size       = GetDofsSize();

  if ( rResult.size() != dofs_size )
    rResult.resize( dofs_size, false );

  switch(mStepVariable)
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

void UpdatedLagrangianSegregatedFluidElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    this->SetProcessInformation(rCurrentProcessInfo);

    FluidElement::InitializeExplicitContributions();

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {

          for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                                                               GetGeometry(),
                                                               row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                                                               rCurrentProcessInfo );
          break;
        }
      case PRESSURE_STEP:
        {
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }


    KRATOS_CATCH( "" )

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

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    this->SetProcessInformation(rCurrentProcessInfo);

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          FluidElement::FinalizeSolutionStep(rCurrentProcessInfo);
          break;
        }
      case PRESSURE_STEP:
        {
          //set as VELOCITY STEP for gauss point calculations:
          mStepVariable = VELOCITY_STEP;
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

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

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        FluidElement::InitializeElementData(rVariables,rCurrentProcessInfo);

        break;
      }
    case PRESSURE_STEP:
      {

        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        //initialize element variables
        rVariables.N.resize(number_of_nodes,false);
        rVariables.L.resize(dimension,dimension,false);
        rVariables.F.resize(dimension,dimension,false);
        rVariables.DN_DX.resize(number_of_nodes,dimension,false);
        rVariables.DeltaPosition.resize(number_of_nodes,dimension,false);

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

  const GeometryType&  rGeometry = GetGeometry();
  //Calculate Delta Position
  ElementUtilities::CalculateDeltaPosition(rVariables.DeltaPosition,rGeometry);

  //set variables including all integration points values

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateMaterialResponse(ElementDataType& rVariables,
                                                                        ConstitutiveLaw::Parameters& rValues,
                                                                        const int & rPointNumber)
{
    KRATOS_TRY

    switch( mStepVariable )
    {
      case VELOCITY_STEP:
        {
          FluidElement::CalculateMaterialResponse(rVariables,rValues,rPointNumber);
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

    switch(mStepVariable)
    {
      case VELOCITY_STEP:
        {
          //Parent to reference configuration
          rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

          const GeometryType&  rGeometry = GetGeometry();
          const SizeType dimension = GetGeometry().WorkingSpaceDimension();

          //Compute the deformation matrix B
          ElementUtilities::CalculateLinearDeformationMatrix(rVariables.B, rGeometry, rVariables.DN_DX);

          //Calculate velocity gradient matrix
          ElementUtilities::CalculateVelocityGradient( rVariables.L, rGeometry, rVariables.DN_DX, rVariables.Alpha );

          //Compute symmetric spatial velocity gradient [DN_DX = dN/dx_n*1] stored in a vector
          ElementUtilities::CalculateSymmetricVelocityGradientVector( rVariables.L, rVariables.StrainVector, dimension );

          break;
        }
      case PRESSURE_STEP:
        {

          const GeometryType&  rGeometry = GetGeometry();
          //Calculate velocity gradient matrix
          ElementUtilities::CalculateVelocityGradient( rVariables.L, rGeometry, rVariables.DN_DX, rVariables.Alpha );

          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::SetElementData(ElementDataType& rVariables,
                                                             ConstitutiveLaw::Parameters& rValues,
                                                             const int & rPointNumber)
{
  //to be accurate calculus must stop
  if(rVariables.detF<0){

    if(this->IsNot(SELECTED)){

      KRATOS_WARNING(" [Element Ignored]") << "ULSFluidElement["<<this->Id()<<"] (|F|=" << rVariables.detF <<")  (Iter:"<<rVariables.GetProcessInfo()[NL_ITERATION_NUMBER]<<")"<<std::endl;
      this->Set(SELECTED,true);

      // SizeType number_of_nodes  = GetGeometry().PointsNumber();

      // for ( SizeType i = 0; i < number_of_nodes; i++ )
      // {
      //   array_1d<double, 3> & CurrentPosition      = GetGeometry()[i].Coordinates();
      //   array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
      //   array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
      //   array_1d<double, 3> PreviousPosition       = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
      //   KRATOS_WARNING("")<<" Node["<<GetGeometry()[i].Id()<<"]: (Position: (pre)"<<PreviousPosition<<",(cur)"<<CurrentPosition<<")"<<std::endl;
      //   //KRATOS_WARNING("")<<" (Displacement: (pre)"<<CurrentDisplacement<<",(cur)"<<PreviousDisplacement<<")"<<std::endl;
      // }
      // for ( SizeType i = 0; i < number_of_nodes; i++ )
      // {
      //   if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
      //     array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
      //     array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
      //     KRATOS_WARNING("")<<" (Contact: (pre)"<<PreContactForce<<",(cur)"<<ContactForce<<")["<<GetGeometry()[i].Id()<<"]"<<std::endl;
      //   }

      // }

      // KRATOS_ERROR<<" [Element Failed] ["<<this->Id()<<"]"<<std::endl;

    }

    rVariables.detJ = 0;

  }
  else{

    if(this->Is(SELECTED) && this->Is(ACTIVE)){
      this->Set(SELECTED,false);
      KRATOS_WARNING("")<<" Undo SELECTED ULSFluidElement "<<this->Id()<<std::endl;
    }
  }



  Flags& ConstitutiveLawOptions = rValues.GetOptions();
  ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);

  rValues.SetDeterminantF(rVariables.detF);
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
        // operation performed: add stiffness term to the rLefsHandSideMatrix
        this->CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );

        // operation performed: add Kg to the rLefsHandSideMatrix
        // this->CalculateAndAddKvvg( rLeftHandSideMatrix, rVariables );

        break;
      }
    case PRESSURE_STEP:
      {
        // operation performed: add Kpp to the rLefsHandSideMatrix
        this->CalculateAndAddKpp( rLeftHandSideMatrix, rVariables );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

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

        break;
      }
    case PRESSURE_STEP:
      {
        // operation performed: add PressureForces to the rRightHandSideVector
        this->CalculateAndAddPressureForces( rRightHandSideVector, rVariables );

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables)
{
  KRATOS_TRY


  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        FluidElement::CalculateAndAddDynamicLHS(rLeftHandSideMatrix, rVariables);

        //KRATOS_WARNING("")<<" DynamicLHS "<<rLeftHandSideMatrix<<std::endl;

        break;
      }
    case PRESSURE_STEP:
      {

        const unsigned int MatSize = this->GetDofsSize();

        if(rLeftHandSideMatrix.size1() != MatSize)
          rLeftHandSideMatrix.resize (MatSize, MatSize, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix( MatSize, MatSize ); //resetting LHS

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables)
{
  KRATOS_TRY

  switch(mStepVariable)
  {
    case VELOCITY_STEP:
      {
        FluidElement::CalculateAndAddDynamicRHS(rRightHandSideVector, rVariables);

        //KRATOS_WARNING("")<<" DynamicRHS "<<rRightHandSideVector<<std::endl;

        break;
      }
    case PRESSURE_STEP:
      {
        const unsigned int MatSize = this->GetDofsSize();

        if ( rRightHandSideVector.size() != MatSize )
          rRightHandSideVector.resize( MatSize, false );

        noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

        break;
      }
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << mStepVariable << std::endl;
  }

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
      MeanPressure += rVariables.N[i] * ( this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) * rVariables.Alpha +  this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE,1) * (1.0 - rVariables.Alpha) );
    }

    this->AddVolumetricPart(rVariables.StressVector, MeanPressure);

    FluidElement::CalculateAndAddInternalForces(rRightHandSideVector, rVariables);

    this->RemoveVolumetricPart(rVariables.StressVector, MeanPressure);

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

  SizeType number_of_nodes = GetGeometry().PointsNumber();

  // 1.- Calculate Stiffness Matrix to get the bulk correction
  const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];

  double BulkFactor = GetProperties()[BULK_MODULUS] * TimeStep;

  double StiffnessFactor = 0.0;

  this->CalculateStiffnessFactor(rLeftHandSideMatrix, rVariables, BulkFactor, StiffnessFactor);


  // 2.- Estimate the bulk correction coefficient and the new tangent
  double MassFactor = GetGeometry().DomainSize() * GetProperties()[DENSITY] / double(number_of_nodes);

  BulkFactor = 0.0;
  if(StiffnessFactor!=0 && MassFactor!=0)
    BulkFactor = GetProperties()[BULK_MODULUS] * TimeStep * ( MassFactor * 2.0 / (TimeStep * StiffnessFactor) );

  // Add volumetric part to the constitutive tensor to compute:
  // Kvvm = Kdev + Kvol(quasi incompressible)
  // add fluid bulk modulus for quasi-incompressibility
  // (can be implemented in a ConstituiveModel for Quasi-Incompressible Newtonian Fluid).
  this->AddVolumetricPart(rVariables.ConstitutiveMatrix, BulkFactor);

  rLeftHandSideMatrix.clear();
  FluidElement::CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );

  rLeftHandSideMatrix *= rVariables.Alpha;

  this->RemoveVolumetricPart(rVariables.ConstitutiveMatrix, BulkFactor);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateStiffnessFactor(MatrixType& rLeftHandSideMatrix,
                                                                       ElementDataType& rVariables,
                                                                       const double& rBulkFactor,
                                                                       double& rStiffnessFactor)

{
  KRATOS_TRY

  // First approach:
  // {
  // // add fluid bulk modulus
  // this->AddVolumetricPart(rVariables.ConstitutiveMatrix,rBulkFactor);

  // FluidElement::CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables );

  // rLeftHandSideMatrix *= rVariables.Alpha;

  // //remove fluid bulk modulus
  // this->RemoveVolumetricPart(rVariables.ConstitutiveMatrix,rBulkFactor);

  // //a. Proposed calculation:
  // // this->CalculateDenseMatrixMeanValue(rLeftHandSideMatrix,rStiffnessFactor);

  // //b. Alternative calculation:
  // this->CalculateLumpedMatrixMeanValue(rLeftHandSideMatrix,rStiffnessFactor);
  // //rStiffnessFactor *= 0.75;
  // }

  //std::cout<<" StiffnessFactor A "<<rStiffnessFactor<<std::endl;
  // Seccond approach:
  // {
  rStiffnessFactor = 0;

  double diagonal = 0;
  for(SizeType i=0; i<rVariables.B.size2(); ++i)
  {
    diagonal = 0;
    for(SizeType j=0; j<rVariables.B.size1(); ++j)
    {
      for(SizeType k=0; k<rVariables.B.size1(); ++k)
      {
        diagonal += rVariables.B(j,i) * (rVariables.ConstitutiveMatrix(j,k) + rBulkFactor) * rVariables.B(k,i);
      }
    }

    rStiffnessFactor += fabs(diagonal);
  }

  double factor = 1.0/ double(GetGeometry().WorkingSpaceDimension() * rLeftHandSideMatrix.size1());

  rStiffnessFactor *= factor * rVariables.IntegrationWeight * rVariables.Alpha;
  // }

  //std::cout<<" StiffnessFactor B "<<rStiffnessFactor<<std::endl;

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

void UpdatedLagrangianSegregatedFluidElement::AddVolumetricPart(Matrix& rConstitutiveMatrix, const double& rBulkFactor)
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

void UpdatedLagrangianSegregatedFluidElement::RemoveVolumetricPart(Matrix& rConstitutiveMatrix, const double& rBulkFactor)
{
  KRATOS_TRY

  double BulkFactor = -rBulkFactor;
  this->AddVolumetricPart(rConstitutiveMatrix, BulkFactor);

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::AddVolumetricPart(Vector& rStressVector, const double& rMeanPressure)
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

void UpdatedLagrangianSegregatedFluidElement::RemoveVolumetricPart(Vector& rStressVector, const double& rMeanPressure)
{
  KRATOS_TRY

  double MeanPressure = -rMeanPressure;
  this->AddVolumetricPart(rStressVector,MeanPressure);

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
  const double& BulkModulus = GetProperties()[BULK_MODULUS];
  const double& Density     = GetProperties()[DENSITY];
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

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
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
    const double& Viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double& Density   = GetProperties()[DENSITY];

    //h_n (normal h)
    Matrix D(dimension,dimension);
    noalias(D) = 0.5 * (trans(rVariables.L)+rVariables.L);

    for( SizeType i=0; i<Faces.size(); ++i ){

      GetFaceNormal(Faces[i], rVariables, Normal);

      GetFaceWeight(Faces[i], rVariables, SideWeight, side_normal_size);

      Vector Acceleration (dimension);
      noalias(Acceleration) = ZeroVector(dimension);

      ProjectionVelocityGradient = inner_prod(Normal, prod(D,Normal));

      BoundFactor = rVariables.Tau * 2.0 / side_normal_size;

      BoundFactorA = rVariables.Tau * Density;
      BoundFactorB = rVariables.Tau * 4.0 * ProjectionVelocityGradient * Viscosity / side_normal_size;

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
        //rRightHandSideVector[Faces[i][j]] -=  SideWeight * BoundFactor * rVariables.N[Faces[i][j]] * rGeometry[Faces[i][j]].FastGetSolutionStepValue(PRESSURE,0);

        //(reduced integration)
        for( SizeType k=0; k<Faces[i].size(); ++k ){
          rRightHandSideVector[Faces[i][j]] -= SideWeight * BoundFactor * rVariables.N[Faces[i][j]] * rVariables.N[Faces[i][k]] * rGeometry[Faces[i][k]].FastGetSolutionStepValue(PRESSURE,0);
        }
      }

    }

  }

  // Add Divergence and volume acceleration vector
  double TraceVelocityGradient = 0;
  for (SizeType i = 0; i < dimension; ++i)
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
      rRightHandSideVector[i] -= rVariables.Tau * rVariables.IntegrationWeight * rVariables.DN_DX(i,j) * VolumeForce[j];
    }
  }

  // Add Dynamic Bulk Vector
  const double& BulkModulus = GetProperties()[BULK_MODULUS];
  const double& Density     = GetProperties()[DENSITY];

  double MassFactor = rVariables.IntegrationWeight / BulkModulus;
  double BulkFactor = MassFactor * Density * rVariables.Tau;

  // (LUMPED)
  // double coefficient = rGeometry.IntegrationPointsNumber() * (1 + dimension); //integration points independent
  // MassFactor /= coefficient;
  // BulkFactor /= coefficient;
  // const double& TimeStep = rVariables.GetProcessInfo()[DELTA_TIME];
  // for( SizeType i=0; i<number_of_nodes; ++i)
  // {
  //   rRightHandSideVector[i] -= MassFactor * (1.0/ TimeStep) * rGeometry[j].FastGetSolutionStepValue(PRESSURE_VELOCITY);
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
        rRightHandSideVector[i] -= (StabilizationFactor * rVariables.DN_DX(i,k) * rVariables.DN_DX(j,k)) * rGeometry[j].FastGetSolutionStepValue(PRESSURE,0);
      }
    }
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateStabilizationTau(ElementDataType& rVariables)

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

    // Get element properties
    const double& Density   = GetProperties()[DENSITY];
    const double& Viscosity = GetProperties()[DYNAMIC_VISCOSITY];
    const double& TimeStep  = rVariables.GetProcessInfo()[DELTA_TIME];

    // Get element size
    double element_size = rGeometry.AverageEdgeLength();

    rVariables.Tau = (element_size * element_size * TimeStep) / ( Density * mean_velocity * TimeStep * element_size + Density * element_size * element_size +  8.0 * Viscosity * TimeStep );
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
  unsigned int face=0;
  for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
  {
    if (ne->Id() == this->Id())  // If there is no shared element in face nf (the Id coincides)
    {
      std::vector<SizeType> Nodes;
      unsigned int WallNodes  = 0;
      unsigned int InletNodes = 0;
      unsigned int SliverNodes = 0;
      unsigned int FreeSurfaceNodes = 0;

      for(unsigned int i = 1; i < NodesInFaces.size1(); ++i)
      {
        Nodes.push_back(NodesInFaces(i,face));  //set boundary nodes
        if(rGeometry[NodesInFaces(i,face)].Is(RIGID) || rGeometry[NodesInFaces(i,face)].Is(SOLID)){
          ++WallNodes;
        }
        if(rGeometry[NodesInFaces(i,face)].Is(FREE_SURFACE)){
          ++FreeSurfaceNodes;
        }
        if(rGeometry[NodesInFaces(i,face)].Is(INLET)){
          ++InletNodes;
        }
        if(rGeometry[NodesInFaces(i,face)].Is(SELECTED)){
          ++SliverNodes;
        }
      }
      if( SliverNodes != Nodes.size() )
        if( WallNodes < Nodes.size() &&  InletNodes < Nodes.size() )
          if( FreeSurfaceNodes == Nodes.size() )
            Faces.push_back(Nodes);
    }

    face++;
  }

  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetFaceNormal(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, Vector& rNormal)
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  // for triangles and tetrahedra
  if( rNormal.size() != dimension )
    rNormal.resize(dimension,false);

  noalias(rNormal) = ZeroVector(dimension);
  for( SizeType j=0; j<rFace.size(); ++j )
  {
    for(unsigned int d=0; d<dimension; ++d)
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

void UpdatedLagrangianSegregatedFluidElement::GetFaceWeight(const std::vector<SizeType>& rFace, const ElementDataType & rVariables, double& rWeight, double& rNormalSize)
{
  KRATOS_TRY

  const SizeType dimension = GetGeometry().WorkingSpaceDimension();

  // for triangles and tetrahedra
  Vector An(dimension);
  noalias(An) = ZeroVector(dimension);
  for( SizeType j=0; j<rFace.size(); ++j )
  {
    for(unsigned int d=0; d<dimension; ++d)
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

void UpdatedLagrangianSegregatedFluidElement::GetFaceNormal(const std::vector<SizeType>& rFace, Vector& rNormal)
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

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& UpdatedLagrangianSegregatedFluidElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1.0 / (rVariables.detF);

    return rVolumeChange;

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
