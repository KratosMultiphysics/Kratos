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
#include "custom_elements/solid_elements/fluid_element.hpp"
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
  mFinalizedStep = true; // the creation is out of the time step, it must be true 
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::UpdatedLagrangianSegregatedFluidElement( UpdatedLagrangianSegregatedFluidElement const& rOther)
    :FluidElement(rOther)
    ,mFinalizedStep(rOther.mFinalizedStep)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement&  UpdatedLagrangianSegregatedFluidElement::operator=(UpdatedLagrangianSegregatedFluidElement const& rOther)
{
    FluidElement::operator=(rOther);
    
    mFinalizedStep = rOther.mFinalizedStep;

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianSegregatedFluidElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianSegregatedFluidElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
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
    
    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());
       
    return Element::Pointer( new UpdatedLagrangianSegregatedFluidElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianSegregatedFluidElement::~UpdatedLagrangianSegregatedFluidElement()
{
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetHistoricalVariables( ElementVariables& rVariables, const double& rPointNumber )
{
    //Deformation Gradient F ( set to identity )
    unsigned int size =  rVariables.F.size1();

    rVariables.detF = 1;
    noalias(rVariables.F) = IdentityMatrix(size);

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_X ) );
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Y ) );

            if( dimension == 3 )
              rElementalDofList.push_back( GetGeometry()[i].pGetDof( VELOCITY_Z ) );
          }
        }
      case PRESSURE_STEP:
        {
          for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
          {
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ) );
          }
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
          {
            rResult[i] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

}


//*********************************DISPLACEMENT***************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }

}


//************************************VELOCITY****************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_VELOCITY, Step );
          }
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }  
 
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = GetDofsSize();
    
    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          unsigned int index = 0;
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
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
          for ( unsigned int i = 0; i < number_of_nodes; i++ )
          {
            rValues[i]     = GetGeometry()[i].GetSolutionStepValue( PRESSURE_ACCELERATION, Step );
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

    FluidElement::InitializeSolutionStep(rCurrentProcessInfo);

    mFinalizedStep = false;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    FluidElement::FinalizeSolutionStep(rCurrentProcessInfo);

    mFinalizedStep = true;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::SetElementVariables(ElementVariables& rVariables,
                                                                  ConstitutiveLaw::Parameters& rValues,
                                                                  const int & rPointNumber)
{

    //to take in account previous step for output print purposes
    unsigned int step = 0;
    if( mFinalizedStep ){
      step = 1;
      this->GetHistoricalVariables(rVariables,rPointNumber);
    }
  
    if(rVariables.detF<0){
        
	std::cout<<" Element: "<<this->Id()<<std::endl;

	unsigned int number_of_nodes = GetGeometry().PointsNumber();

	for ( unsigned int i = 0; i < number_of_nodes; i++ )
	  {
	    array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
	    array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	    array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	    array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
	    std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<" (Cur: "<<CurrentPosition<<") "<<std::endl;
	    std::cout<<" ---Disp: "<<CurrentDisplacement<<" (Pre: "<<PreviousDisplacement<<")"<<std::endl;
	  }

        KRATOS_ERROR << " LARGE DISPLACEMENT SEGREGATED FLUID ELEMENT INVERTED: |F|<0  detF = " << rVariables.detF << std::endl;
    }



    //Compute strain rate measures if they are required by the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    mConstitutiveLawVector[rPointNumber]->GetLawFeatures(LawFeatures);
    
    bool strain_rate_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
      if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Velocity_Gradient)
	strain_rate_measure = true;
    }

    if( strain_rate_measure ){     
      //Compute symmetric spatial velocity gradient [DN_DX = dN/dx_n*1] stored in a vector
      this->CalculateVelocityGradientVector( rVariables.StrainVector, rVariables.DN_DX, step );
      Flags &ConstitutiveLawOptions=rValues.GetOptions();
      ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    }

    //Compute F and detF (from 0 to n+1) : store it in H variable and detH
    rVariables.detH = rVariables.detF * rVariables.detF0;
    noalias(rVariables.H) = prod( rVariables.F, rVariables.F0 );
    
    rValues.SetDeterminantF(rVariables.detH);
    rValues.SetDeformationGradientF(rVariables.H);
    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY
       
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix(); 

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          // operation performed: add Km to the rLefsHandSideMatrix
          this->CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
          
          // operation performed: add Kg to the rLefsHandSideMatrix
          // this->CalculateAndAddKvvg( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

          // operation performed: add Bulk term matrix to rLeftHandSideMatrix
          this->CalculateAndAddKbulk( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

          break;
        }
      case PRESSURE_STEP:
        {
          // operation performed: add Kpp to the rLefsHandSideMatrix
          this->CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
    
    //KRATOS_WATCH( rLeftHandSideMatrix )
  
    KRATOS_CATCH( "" )  
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY
       
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 

    switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
    {
      case VELOCITY_STEP:
        {
          // operation performed: add InternalForces to the rRightHandSideVector
          this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
          
          // operation performed: add ExternalForces to the rRightHandSideVector
          this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );
          

          break;
        }
      case PRESSURE_STEP:
        {
          // operation performed: add PressureForces to the rRightHandSideVector
          this->CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight );
          
          break;
        }
      default:
        KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
    }
    
    //KRATOS_WATCH( rRightHandSideVector )
  
    KRATOS_CATCH( "" ) 
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
                                                                            ElementVariables & rVariables,
                                                                            double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& rGeometry = GetGeometry();    
    unsigned int number_of_nodes = rGeometry.PointsNumber();      
        
    // Add Boundary Vector
    this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);


    // Add Stabilization Vector
    for (SizeType i = 0; i < NumNodes; ++i)
    {         
      // RHS contribution
      // Velocity divergence   
      rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
      this->AddStabilizationNodalTermsRHS(rRightHandSideVector,Tau,Density,GaussWeight,rDN_DX,i);
    }

    // Add stabilization and boundary term RHS
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,PressureValuesForRHS);

    
    // Add Bulk Vector
    double MassFactor = GetGeometry().DomainSize() * GetProperties()[DENSITY] / double(number_of_nodes);
    const double& BulkModulus = GetProperties()[BULK_MODULUS];
    const double& Density     = GetProperties()[DENSITY];
   
    double BulkFactor = MassFactor*Tau*Density/(BulkModulus*rVariables.TimeStep*rVariables.TimeStep);
    
    // for pressure velocity only or x2 if adding pressure acceleration
    // BulkFactor *= 2;

    double coefficient = rGeometry.IntegrationPointsNumber() + dimension + 1;
    BulkFactor /= coefficient;
    
    for( unsigned int i=0; i<number_of_nodes; ++i)
    {
      rRightHandSideVector[i] -= BulkFactor * rGeometry[i].FastGetSolutionStepValue(PRESSURE);
      rRightHandSideVector[i] -= BulkFactor * rVariables.TimeStep * rGeometry[i].FastGetSolutionStepValue(PRESSURE_VELOCITY);
    }
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

unsigned int UpdatedLagrangianSegregatedFluidElement::GetDofsSize()
{
  KRATOS_TRY
     
  const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
  const unsigned int number_of_nodes = GetGeometry().PointsNumber();
  
  unsigned int size = 0;
  switch(StepType(rCurrentProcessInfo[SEGREGATED_STEP]))
  {
    case VELOCITY_STEP:
      size = number_of_nodes * dimension; //size for velocity
      break;
    case PRESSURE_STEP:
      size = number_of_nodes; //size for pressure
      break;
    default:
      KRATOS_ERROR << "Unexpected value for SEGREGATED_STEP index: " << rCurrentProcessInfo[SEGREGATED_STEP] << std::endl;
      
  }
  return size;   
  
  KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKvvg(MatrixType& rLeftHandSideMatrix,
                                                                  ElementVariables& rVariables,
                                                                  double& rIntegrationWeight)

{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
    MathUtils<double>::ExpandAndAddReducedMatrix( rLeftHandSideMatrix, ReducedKg, dimension );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKbulk(MatrixType& rLeftHandSideMatrix,
                                                                   ElementVariables& rVariables,
                                                                   double& rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double StiffnessFactor = 0.0;    
    for ( unsigned int i = 0; i < rLeftHandSideVector.size1(); i++ )
    {
      StiffnessFactor += rLeftHandSideVector[i][i];
    }

    StiffnessFactor /= double(number_of_nodes,dimension);

    double MassFactor = GetGeometry().DomainSize() * GetProperties()[DENSITY] / double(number_of_nodes);

    double BulkFactor = 0.0;
    if(StiffnessFactor!=0 && MassFactor!=0)
      BulkFactor = MassFactor*4.0/(rVariables.TimeStep*StiffnessFactor);
    
    BulkFactor += GetProperties()[BULK_MODULUS]*rVariables.TimeStep; //add fluid bulk modulus for quasi-incompressibility (must be implemented in a ConstituiveModel for Quasi-Incompressible Newtonian Fluid.    
    for ( unsigned int i = 0; i < dimension; i++ )
    {
      for ( unsigned int j = i+1; j < dimension; j++ )
      {
        rVariables.ConstitutiveMatrix[i,j] += BulkFactor;
      }
    }

    for ( unsigned int i = 0; i < rVariables.ConstitutiveMatrix.size1(); i++ )
    {
      rVariables.ConstitutiveMatrix[i,i] += BulkFactor;
    }

    
    this->CalculateAndAddKvvm( rLeftHandSideMatrix, rVariables, rIntegrationWeight );
    
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianSegregatedFluidElement::CalculateAndAddKpp(MatrixType& rLeftHandSideMatrix,
                                                                 ElementVariables& rVariables,
                                                                 double& rIntegrationWeight)

{
    KRATOS_TRY

    GeometryType& rGeometry = GetGeometry();    
    unsigned int number_of_nodes = rGeometry.PointsNumber();

    // Get mean velocity
    double mean_velocity = 0;
    for( unsigned int i=0; i<number_of_nodes; ++i)
    {
      const array_1d<double,3>& Velocity = rGeometry[i].FastGetSolutionStepValue(VELOCITY);
      mean_velocity += norm_2(Velocity);
    }
    mean_velocity /= double(number_of_nodes);

    // Calculate FIC stabilization coefficient
    double Tau = 0;
    if( mean_velocity != 0 ){

      // Get element properties
      const double& Density   = GetProperties()[DENSITY];
      const double& Viscosity = GetProperties()[DYNAMIC_VISCOSITY];
      // Get element size
      double element_size = rGeometry.AverageEdgeLength();

      Tau = (element_size * element_size * rVariables.TimeStep) / ( Density * mean_velocity * rVariables.TimeStep * element_size + Density * element_size * element_size +  8.0 * Viscosity * rVariables.TimeStep );
      
    }

    // Add Boundary Matrix
    double BoundFactor = Tau*4.0*rIntegrationWeight/(element_size*element_size);
     
    DenseMatrix<unsigned int> NodesInFaces;
    rGeometry.NodesInFaces(NodesInFaces);
   
    for( unsigned int i=0; i<NodesInFaces.size2(); ++i ){
      bool free_surface = true;
      for( unsigned int j=1; j<NodesInFaces.size1(); ++j ){
        if( rGeometry[NodesInFaces(j,i)].IsNot(FREE_SURFACE) ){
          free_surface = false;
          break;
        }
      }
      if( free_surface ){
        for( unsigned int j=1; j<NodesInFaces.size1(); ++j ){
          if( rGeometry[NodesInFaces(j,i)].IsNot(INLET) ){
            rLeftHandSideMatrix(NodesInFaces(j,i),NodesInFaces(j,i)) += BoundFactor / double(NodesInFaces.size1());
          }
        }
      }       
    }


    // Add Stabilized Laplacian Matrix
    double StabilizationFactor = Tau * rIntegrationWeight;
    for( unsigned int i=0; i<number_of_nodes; ++i)
    {
      for( unsigned int j=0; j<number_of_nodes; ++j)
      {
        for ( unsigned int k = 0; k < dimension; ++k )
        {
          rLeftHandSideMatrix(i,j) += StabilizationFactor * rVariables.DN_DX(i,k) * rVariables.DN_DX(j,k);
        }
      }
    }

    // Add Bulk Matrix
    double MassFactor = GetGeometry().DomainSize() * GetProperties()[DENSITY] / double(number_of_nodes);
    const double& BulkModulus = GetProperties()[BULK_MODULUS];
    const double& Density     = GetProperties()[DENSITY];
   
    double BulkFactor = MassFactor*Tau*Density/(BulkModulus*rVariables.TimeStep*rVariables.TimeStep);
    
    // for pressure velocity only or x2 if adding pressure acceleration
    // BulkFactor *= 2;

    double coefficient = rGeometry.IntegrationPointsNumber() + dimension + 1;
    for( unsigned int i=0; i<number_of_nodes; ++i)
    {
      rLeftHandSideMatrix(i,i) += BulkFactor / coefficient;
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
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();
    if( dimension == 2 )
    {
      if( LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRAIN_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::PLANE_STRESS_LAW) && LawFeatures.mOptions.IsNot(ConstitutiveLaw::AXISYMMETRIC_LAW) )
	KRATOS_ERROR <<  "wrong constitutive law used. This is a 2D element. Expected plane state or axisymmetric :: element id = " << this->Id() << std::endl;
    }

    // Check that the element nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
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
}

void UpdatedLagrangianSegregatedFluidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, FluidElement )
}


} // Namespace Kratos


