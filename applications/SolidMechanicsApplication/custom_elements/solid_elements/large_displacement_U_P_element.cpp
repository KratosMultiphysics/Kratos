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
#include "custom_elements/solid_elements/large_displacement_U_P_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementUPElement::LargeDisplacementUPElement()
    : LargeDisplacementElement()
{
  //DO NOT CALL IT: only needed for Register and Serialization!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementUPElement::LargeDisplacementUPElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : LargeDisplacementElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementUPElement::LargeDisplacementUPElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : LargeDisplacementElement( NewId, pGeometry, pProperties )
{
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

LargeDisplacementUPElement::LargeDisplacementUPElement( LargeDisplacementUPElement const& rOther)
    :LargeDisplacementElement(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

LargeDisplacementUPElement&  LargeDisplacementUPElement::operator=(LargeDisplacementUPElement const& rOther)
{
    LargeDisplacementElement::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer LargeDisplacementUPElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new LargeDisplacementUPElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer LargeDisplacementUPElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_THROW_ERROR( std::logic_error, "calling the default constructor for a large displacement 3D element ... illegal operation!!", "" )

    LargeDisplacementUPElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_THROW_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() )
      }
    
    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());
       
    return Element::Pointer( new LargeDisplacementUPElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

LargeDisplacementUPElement::~LargeDisplacementUPElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************



void LargeDisplacementUPElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
    }
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }

    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void LargeDisplacementUPElement::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
        }
        else
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
        }

    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void LargeDisplacementUPElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            rValues[index + 3] = 0;
        }
        else
        {
            rValues[index + 2] = 0;
        }
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void LargeDisplacementUPElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       dofs_size       = number_of_nodes * dimension + number_of_nodes;

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0;
        }
        else
        {
            rValues[index + 2] = 0;
        }
    }

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::InitializeElementVariables (ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    LargeDisplacementElement::InitializeElementVariables(rVariables,rCurrentProcessInfo);

    //stabilization factor
    double StabilizationFactor = 1.0;
    if( GetProperties().Has(STABILIZATION_FACTOR) ){
      StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    }
    else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR) ){
      StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR];
    }
    GetProperties().SetValue(STABILIZATION_FACTOR, StabilizationFactor);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(LargeDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS
    }
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, double& rIntegrationWeight)
{

    //contributions of the stiffness matrix calculated on the reference configuration
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    // operation performed: add Km to the rLefsHandSideMatrix

    //respect to the current configuration n+1
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kg to the rLefsHandSideMatrix
    CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kpp Stab to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{

    //contribution of the internal and external forces
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 
      
    // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, rIntegrationWeight);

    //KRATOS_WATCH( rRightHandSideVector )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        ElementVariables& rVariables,
        Vector& rVolumeForce,
        double& rIntegrationWeight)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // VectorType Fh=rRightHandSideVector;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[indexup + j] += rIntegrationWeight * rVariables.N[i] * rVolumeForce[j];
        }

    }

    // std::cout<<std::endl;
    // std::cout<<" Fext "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        ElementVariables & rVariables,
        double& rIntegrationWeight
                                                              )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // VectorType Fh=rRightHandSideVector;

    Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexup = dimension * i + i;
        unsigned int indexu  = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
        }
    }

    // std::cout<<std::endl;
    // std::cout<<"["<<this->Id()<<"] StressVector "<<rVariables.StressVector<<std::endl;
    // std::cout<<" Fint "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

double& LargeDisplacementUPElement::CalculatePUCoefficient(double& rCoefficient, ElementVariables & rVariables)
{
  KRATOS_TRY
   
    //Mechanical volumetric:
    
    //Constitutive A:
    //rCoefficient = 0.5*(rVariables.detF0*rVariables.detF0-1)/rVariables.detF0); //(J²-1)/2
    
    //Constitutive B:
    rCoefficient = (std::log(rVariables.detF0)/rVariables.detF0);  //(ln(J))
  
    //Thermal volumetric:

    double ThermalExpansionCoefficient = 0;
    if( GetProperties().Has(THERMAL_EXPANSION_COEFFICIENT) ){
      ThermalExpansionCoefficient = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    }
    
    double DeltaTemperature     = 0;
    double ReferenceTemperature = 0;
    double CurrentTemperature   = 0;

    if( GetProperties().Has(REFERENCE_TEMPERATURE) )
      ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];

    
    int count = 0;
    for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
      {
	if(this->GetGeometry()[j].SolutionStepsDataHas( TEMPERATURE ) == true)
	  {
	    CurrentTemperature += rVariables.N[j] * GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE);
	    count++;
	  }	
      }

    if( count == 0 ){
      CurrentTemperature = ReferenceTemperature;
      DeltaTemperature  = 0;
    }
    else{
      DeltaTemperature = CurrentTemperature - ReferenceTemperature;
    }
    
    rCoefficient += 3.0 * ThermalExpansionCoefficient * ( (1.0 - std::log(rVariables.detF0)) / (rVariables.detF0 * rVariables.detF0) ) * DeltaTemperature;
    
    return rCoefficient;

  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

double& LargeDisplacementUPElement::CalculatePUDeltaCoefficient(double &rDeltaCoefficient, ElementVariables & rVariables)
{

  KRATOS_TRY

    //Mechanical volumetric:

    //Constitutive A:
    //rDeltaCoefficient = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)/2

    //Constitutive B:
    rDeltaCoefficient = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J))

    //Thermal volumetric:

    double ThermalExpansionCoefficient = 0;
    if( GetProperties().Has(THERMAL_EXPANSION_COEFFICIENT) ){
      ThermalExpansionCoefficient = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    }
    
    double DeltaTemperature     = 0;
    double ReferenceTemperature = 0;
    double CurrentTemperature   = 0;

    if( GetProperties().Has(REFERENCE_TEMPERATURE) )
      ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
    
    int count = 0;
    for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
      {
	if(this->GetGeometry()[j].SolutionStepsDataHas( TEMPERATURE ) == true)
	  {
	    CurrentTemperature += rVariables.N[j] * GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE);
	    count++;
	  }	
      }

    if( count == 0 ){
      CurrentTemperature = ReferenceTemperature;
      DeltaTemperature  = 0;
    }
    else{
      DeltaTemperature = CurrentTemperature - ReferenceTemperature;
    }
    
    rDeltaCoefficient += 3 * ThermalExpansionCoefficient * ( (2 * std::log(rVariables.detF0) - 3.0) / (rVariables.detF0 * rVariables.detF0 * rVariables.detF0) ) * DeltaTemperature;
    
    return rDeltaCoefficient;


  KRATOS_CATCH( "" )

}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
							       ElementVariables & rVariables,
							       double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    // VectorType Fh=rRightHandSideVector;

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    //double consistent=1;

    double Coefficient = 0;
    Coefficient = this->CalculatePUCoefficient( Coefficient, rVariables ); 

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);

	    // consistent=1;
	    // if(i==j)
	    //   consistent=2;

	    // if( dimension == 2 ){ //consistent 2D

	    //   rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/12.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D

	    // }
	    // else{
	      
	    //   rRightHandSideVector[indexp] += consistent * (1.0/BulkModulus) * (1.0/20.0) * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //3D
	    // }

	    rRightHandSideVector[indexp] += (1.0/BulkModulus) * rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) ; //2D-3D

 
        }

        rRightHandSideVector[indexp] -=  Coefficient * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);

        indexp += (dimension + 1);
    }


    // std::cout<<std::endl;
    // std::cout<<" Coefficient " <<Coefficient<<" F0 "<<rVariables.detF0<<std::endl;
    // std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
        ElementVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    unsigned int indexp = dimension;

    // VectorType Fh=rRightHandSideVector;
    // std::cout<<" Element "<<this->Id()<<" "<<std::endl;

    //use of this variable for the complete parameter:
    double AlphaStabilization  = 1.0; 
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    double consistent = 1;

    double FactorValue = 8.0; //JMR deffault value	 
    if( dimension == 3 )
      FactorValue = 10.0; //JMC deffault value

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);

	    if( dimension == 2 ){ //consistent 2D

	      consistent=(-1)*AlphaStabilization*FactorValue/(36.0*LameMu);
	      if(i==j)
                consistent=2*AlphaStabilization*FactorValue/(36.0*LameMu);

	      rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
	    }
	    else{

	      consistent=(-1)*AlphaStabilization*FactorValue/(80.0*LameMu);
	      if(i==j)
                consistent=3*AlphaStabilization*FactorValue/(80.0*LameMu);

	      rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D

	    }

	    // std::cout<<" Pressure "<<Pressure<<std::endl;
        }


        indexp += (dimension + 1);
    }


    // std::cout<<std::endl;
    // std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
    // std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary

    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

     // MatrixType Kh=rLeftHandSideMatrix;

    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    // std::cout<<std::endl;
    // std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}




//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    int size = number_of_nodes * dimension;

    Matrix StressTensor = MathUtils<double>::StressVectorToTensor( rVariables.StressVector );
    Matrix ReducedKg = prod( rVariables.DN_DX,  rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized

    Matrix Kuu = zero_matrix<double> (size);
    MathUtils<double>::ExpandAndAddReducedMatrix( Kuu, ReducedKg, dimension );

    // MatrixType Kh=rLeftHandSideMatrix;

    //assemble into rLeftHandSideMatrix the geometric uu contribution:
    unsigned int indexi = 0;
    unsigned int indexj = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                    indexj++;
                }
            }
            indexi++;
        }
    }

    // std::cout<<std::endl;
    // std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //MatrixType Kh=rLeftHandSideMatrix;
    //contributions to stiffness matrix calculated on the reference configuration
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexp  = dimension;
        unsigned int indexup = dimension * i + i;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight * rVariables.detF;
            }
            indexp += (dimension + 1);
        }
    }

     // std::cout<<std::endl;
     // std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables,
        double& rIntegrationWeight)

{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexp = dimension;

    double DeltaCoefficient = 0;
    DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables ); 


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
	      rLeftHandSideMatrix(indexp,indexup+k) +=  DeltaCoefficient  * rVariables.N[i] * rVariables.DN_DX ( j , k ) * rIntegrationWeight * rVariables.detF;

                //std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
            }
        }
        indexp += (dimension + 1);
    }


    // std::cout<<std::endl;
    // std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
        ElementVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

    // MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    //double consistent = 1.0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {
	     
	    // consistent=1;
	    // if(indexpi==indexpj)
	    //   consistent=2;

	    // if( dimension == 2 ){ //consistent 2D
	      
	    //   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

	    // }
	    // else{

	    //   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/20.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D
	      
	    // }

	    rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D-3D

            indexpj += (dimension + 1);
	  }

        indexpi += (dimension + 1);
    }

    // std::cout<<std::endl;
    // std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        ElementVariables & rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    //repasar

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // MatrixType Kh=rLeftHandSideMatrix;

    //contributions to stiffness matrix calculated on the reference configuration
    unsigned int indexpi = dimension;
    
    double AlphaStabilization  = 1.0; 
    double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
    AlphaStabilization *= StabilizationFactor;

    const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
    const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

    double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));

    //Experimental
    // if(LameMu < rVariables.ConstitutiveMatrix(2,2))
    //   LameMu = rVariables.ConstitutiveMatrix(2,2);

    double consistent = 1.0;

    double FactorValue = 8.0; //JMR deffault value	
    if( dimension == 3 )
      FactorValue = 10.0; //JMC deffault value

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    if( dimension == 2 ){ //consistent 2D
	      
	      consistent=(-1)*AlphaStabilization*FactorValue/(36.0*LameMu);
	      if(indexpi==indexpj)
                consistent=2*AlphaStabilization*FactorValue/(36.0*LameMu);

	      rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
	    
	    }
	    else{

	      consistent=(-1)*AlphaStabilization*FactorValue/(80.0*LameMu);
	      if(indexpi==indexpj)
                consistent=3*AlphaStabilization*FactorValue/(80.0*LameMu);
	      
	      rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D

	    }


            indexpj += (dimension + 1);
	  }

        indexpi += (dimension + 1);
      }

    // std::cout<<std::endl;
    // std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

    KRATOS_CATCH( "" )
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

    // Not Lumped Mass Matrix (numerical integration):

    //reading integration points
    IntegrationMethod CurrentIntegrationMethod = mThisIntegrationMethod; //GeometryData::GI_GAUSS_2; //GeometryData::GI_GAUSS_1;

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( CurrentIntegrationMethod  );

    ElementVariables Variables;
    this->InitializeElementVariables(Variables,rCurrentProcessInfo);

    
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
      //compute element kinematics
      this->CalculateKinematics( Variables, PointNumber );

      //getting informations for integration
      double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
      
      IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

      //compute point volume change
      double PointVolumeChange = 0;
      PointVolumeChange = this->CalculateVolumeChange( PointVolumeChange, Variables );
	
      double CurrentDensity = PointVolumeChange * GetProperties()[DENSITY];

      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      	{
      	  unsigned int indexupi = dimension * i + i;

      	  for ( unsigned int j = 0; j < number_of_nodes; j++ )
      	    {
      	      unsigned int indexupj = dimension * j + j;

      	      for ( unsigned int k = 0; k < dimension; k++ )
      		{
      		  rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * CurrentDensity * IntegrationWeight;
      		}
      	    }
      	}

    }

    // Lumped Mass Matrix:

    // double TotalMass = 0;

    // this->CalculateTotalMass( TotalMass, rCurrentProcessInfo );

    // if ( dimension == 2 ){
    //   if ( this->GetProperties().Has( THICKNESS ) )
    // 	TotalMass *= GetProperties()[THICKNESS];
    // }

    // Vector LumpFact(number_of_nodes);
    // noalias(LumpFact) = ZeroVector(number_of_nodes);

    // LumpFact = GetGeometry().LumpingFactors( LumpFact );

    // for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // {
    //     double temp = LumpFact[i] * TotalMass;

    //     unsigned int indexup = i * dimension + i;

    //     for ( unsigned int j = 0; j < dimension; j++ )
    //     {
    //         rMassMatrix( indexup+j , indexup+j ) = temp;
    //     }
    // }

    //std::cout<<std::endl;
    //std::cout<<" Mass Matrix "<<rMassMatrix<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rDampingMatrix.size1() != MatSize )
        rDampingMatrix.resize( MatSize, MatSize, false );

    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );

    //1.-Calculate StiffnessMatrix:

    MatrixType LHSMatrix  = Matrix();

    this->CalculateLeftHandSide( LHSMatrix, rCurrentProcessInfo );

    MatrixType StiffnessMatrix  = Matrix();

    if ( StiffnessMatrix.size1() != MatSize )
        StiffnessMatrix.resize( MatSize, MatSize, false );

    noalias(StiffnessMatrix) = ZeroMatrix( MatSize, MatSize );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexup = i * dimension + i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
	  StiffnessMatrix( indexup+j , indexup+j ) = LHSMatrix( indexup+j , indexup+j );
        }
    }

    //2.-Calculate MassMatrix:

    MatrixType MassMatrix  = Matrix();

    this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );
    
    
    //3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
    double alpha = 0;
    if( GetProperties().Has(RAYLEIGH_ALPHA) ){
      alpha = GetProperties()[RAYLEIGH_ALPHA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
      alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
    }

    double beta  = 0;
    if( GetProperties().Has(RAYLEIGH_BETA) ){
      beta = GetProperties()[RAYLEIGH_BETA];
    }
    else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
      beta = rCurrentProcessInfo[RAYLEIGH_BETA];
    }

    //4.-Compose the Damping Matrix:
   
    //Rayleigh Damping Matrix: alpha*M + beta*K
    rDampingMatrix  = alpha * MassMatrix;
    rDampingMatrix += beta  * StiffnessMatrix;


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
  KRATOS_TRY
    
  this->CalculateMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

  //KRATOS_WATCH( rLeftHandSideMatrix )
  
  KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void LargeDisplacementUPElement::CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector, ElementVariables& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
  KRATOS_TRY
      
  //mass matrix
  MatrixType LeftHandSideMatrix = Matrix();
  this->CalculateMassMatrix(LeftHandSideMatrix, rCurrentProcessInfo);

  //acceleration vector
  Vector CurrentAccelerationVector( LeftHandSideMatrix.size1() );
  noalias(CurrentAccelerationVector) = ZeroVector( LeftHandSideMatrix.size1() );
  this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);
  
  double AlphaM = 0.0;
  if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
    AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    Vector PreviousAccelerationVector(LeftHandSideMatrix.size1());
    noalias(PreviousAccelerationVector) = ZeroVector( LeftHandSideMatrix.size1() );
    this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
    CurrentAccelerationVector *= (1.0-AlphaM);
    CurrentAccelerationVector +=  AlphaM * (PreviousAccelerationVector);
  }
   
  noalias(rRightHandSideVector) = prod( LeftHandSideMatrix, CurrentAccelerationVector );
  
  //KRATOS_WATCH( rRightHandSideVector )
  
  KRATOS_CATCH( "" )    
}

//************************************************************************************
//************************************************************************************

int  LargeDisplacementUPElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = LargeDisplacementElement::Check(rCurrentProcessInfo);
    
    // Check compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    if(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW))
      KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " Large Displacements U_P" )

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(PRESSURE);

    return ErrorCode;

    KRATOS_CATCH( "" );
}

//************************************************************************************
//************************************************************************************


void LargeDisplacementUPElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeDisplacementElement )
}

void LargeDisplacementUPElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeDisplacementElement )
}


} // Namespace Kratos


