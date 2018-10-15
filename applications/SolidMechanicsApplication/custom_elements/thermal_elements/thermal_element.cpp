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
#include "custom_elements/thermal_elements/thermal_element.hpp"

#include "solid_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( ThermalElement, COMPUTE_RHS_VECTOR,       0 );
KRATOS_CREATE_LOCAL_FLAG( ThermalElement, COMPUTE_LHS_MATRIX,       1 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalElement::ThermalElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalElement::ThermalElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ThermalElement::ThermalElement( ThermalElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ThermalElement&  ThermalElement::operator=(ThermalElement const& rOther)
{
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer ThermalElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<ThermalElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ThermalElement::~ThermalElement()
{
}

//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

ThermalElement::IntegrationMethod ThermalElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}


//************************************************************************************
//************************************************************************************

void ThermalElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
    }

}

//************************************************************************************
//************************************************************************************

void ThermalElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    unsigned int number_of_nodes  = GetGeometry().size();

    if ( rResult.size() != number_of_nodes )
        rResult.resize( number_of_nodes, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	rResult[i] = GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void ThermalElement::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = GetGeometry()[i].GetSolutionStepValue( TEMPERATURE, Step );
    }
}

//************************************VELOCITY****************************************
//************************************************************************************

void ThermalElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = 0;
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void ThermalElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    if ( rValues.size() != MatSize ) rValues.resize( MatSize, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        rValues[i] = 0;
    }
}

//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void ThermalElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                  std::vector<double>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void ThermalElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
                                                  std::vector<Vector>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//*********************************SET MATRIX VALUE***********************************
//************************************************************************************

void ThermalElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                                  std::vector<Matrix>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************

void ThermalElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
                                                  std::vector<double>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//**********************************GET VECTOR VALUE**********************************
//************************************************************************************

void ThermalElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
                                                  std::vector<Vector>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void ThermalElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
                                                  std::vector<Matrix>& rValues,
                                                  const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

void ThermalElement::Initialize()
{
    KRATOS_TRY

    // std::cout<<" Thermal Initialization "<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ThermalElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ThermalElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY


    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{

  const unsigned int number_of_nodes = GetGeometry().size();
  const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();


  rVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

  rVariables.DN_DX.resize( number_of_nodes, dimension );

  //reading shape functions
  rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

  //reading shape functions local gradients
  rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

  //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
  rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

  //Calculate Delta Position
  rVariables.DeltaPosition = CalculateDeltaPosition(rVariables.DeltaPosition);

  //set variables including all integration points values

  //calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d£]
  rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod, rVariables.DeltaPosition );

  //Thermal properties
  rVariables.HeatCapacity     = GetProperties()[HEAT_CAPACITY] * GetProperties()[DENSITY]; //unity factor: [kg/mm3 = 1e-3 Ns2/mm4]
  rVariables.HeatConductivity = GetProperties()[HEAT_CONDUCTIVITY];

  this->CalculateThermalProperties(rVariables);

}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void ThermalElement::CalculateKinematics(GeneralVariables& rVariables,
                                         const double& rPointNumber)

{
    KRATOS_TRY

    //Get the parent coodinates derivative [dN/d£]
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

    //Get the shape functions for the order of the integration method [N]
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    Matrix InvJ;
    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

    //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    Matrix Invj;
    //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n+1]
    MathUtils<double>::InvertMatrix( rVariables.j[rPointNumber], Invj, rVariables.detJ); //overwrites detJ

    //Compute cartesian derivatives
    rVariables.DN_DX = prod( DN_De[rPointNumber] , Invj ); //overwrites DX now is the current position

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber );

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************

Matrix& ThermalElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const GeometryType& rGeometry = GetGeometry();

    ElementUtilities::CalculateDeltaPosition(rDeltaPosition,rGeometry);

    return rDeltaPosition;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateThermalProperties(GeneralVariables& rVariables)
{
  bool TemperatureDependentLaw = false;
  if( GetProperties().Has(TEMPERATURE_DEPENDENT) ){
    TemperatureDependentLaw = GetProperties()[TEMPERATURE_DEPENDENT];
  }

  if( TemperatureDependentLaw ){
    const unsigned int number_of_nodes = GetGeometry().size();

    //double ReferenceTemperature = GetProperties()[REFERENCE_TEMPERATURE];
    double ElementTemperature   = 0;

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
	double NodalTemperature = GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE);
	ElementTemperature += rVariables.N[j] * NodalTemperature;
      }

    if( (ElementTemperature - 273.15) < 0 ){
      ElementTemperature = 273.15;
      std::cout<<" Temperature error on linear thermal properties "<<std::endl;
    }

    double HeatCapacity_a = GetProperties()[HEAT_CAPACITY_A] * GetProperties()[DENSITY];
    double HeatCapacity_b = GetProperties()[HEAT_CAPACITY_B] * GetProperties()[DENSITY];

    rVariables.HeatCapacity     = HeatCapacity_a * (ElementTemperature - 273.15) + HeatCapacity_b;

    double HeatConductivity_a = GetProperties()[HEAT_CONDUCTIVITY_A];
    double HeatConductivity_b = GetProperties()[HEAT_CONDUCTIVITY_B];

    rVariables.HeatConductivity = HeatConductivity_a * (ElementTemperature - 273.15) + HeatConductivity_b;

    //std::cout<<" HeatCapacity "<<rVariables.HeatCapacity<<" HeatConductivity "<<rVariables.HeatConductivity<<std::endl;
  }
}


//************************************************************************************
//************************************************************************************

void ThermalElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					      VectorType& rRightHandSideVector,
					      Flags& rCalculationFlags)
{

  const unsigned int number_of_nodes = GetGeometry().size();

  //resizing as needed the LHS
  unsigned int MatSize = number_of_nodes ;

  if ( rCalculationFlags.Is(ThermalElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
      if ( rLeftHandSideMatrix.size1() != MatSize )
	rLeftHandSideMatrix.resize( MatSize, MatSize, false );

      noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(ThermalElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required

    {
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }
}


//************************************************************************************
//************************************************************************************


void ThermalElement::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix,
					       VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo,
					       Flags& rCalculationFlags )
{
    KRATOS_TRY

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //mechanical evaluation of the material response (add a Element::Pointer as a member to increase speed)
    bool thermo_mechanical = false;

    std::vector<Vector> StressVector;

    std::vector<ConstitutiveLaw::Pointer> ConstitutiveLawVector;

    if( this->Has(MASTER_ELEMENTS) ){

      thermo_mechanical = true;

      if( GetValue(MASTER_ELEMENTS).size() != 1 )
        KRATOS_ERROR<< " Multiple Thermal Element MASTER ELEMENTS "<<GetValue(MASTER_ELEMENTS).size()<<std::endl;

      Element::ElementType& MechanicalElement = GetValue(MASTER_ELEMENTS).back();

      MechanicalElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo);

      MechanicalElement.GetValueOnIntegrationPoints(CONSTITUTIVE_LAW, ConstitutiveLawVector, rCurrentProcessInfo);

    }
    else {

      std::cout<<" NO master element for this thermal element "<<this->Id()<<std::endl;

      thermo_mechanical = false;

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      unsigned int voigtsize = 3;
      if( dimension == 3 )
	{
	  voigtsize  = 6;
	}

      ConstitutiveLawVector.resize( integration_points.size() );
      StressVector.resize( integration_points.size() );

      for ( unsigned int i = 0; i < integration_points.size(); i++ )
	{
	  StressVector[i] = ZeroVector(voigtsize);
	}
    }

    //auxiliary terms
    double HeatSource;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //COMPUTE kinematics B,F,DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

	//calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );

	if ( rCalculationFlags.Is(ThermalElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
          Variables.DeltaPlasticDissipation=0;
          if( thermo_mechanical ){
            //std::cout<<" thermo_mechanical element LHS ["<<this->Id()<<"]"<<std::endl;
            if( ConstitutiveLawVector[PointNumber]->Has(DELTA_PLASTIC_DISSIPATION) ){
              ConstitutiveLawVector[PointNumber]->GetValue(DELTA_PLASTIC_DISSIPATION,Variables.DeltaPlasticDissipation);
              //std::cout<<" with DELTA_PLASTIC_DISSIPATION= "<<Variables.DeltaPlasticDissipation<<std::endl;
            }
          }

          this->CalculateAndAddLHS ( rLeftHandSideMatrix, Variables, IntegrationWeight );
        }

        if ( rCalculationFlags.Is(ThermalElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
          //contribution to external forces
          if( GetProperties().Has(HEAT_SOURCE) ){
            HeatSource = GetProperties()[HEAT_SOURCE];
            // std::cout<<" HeatSource "<<HeatSource<<std::endl;
          }
          
          Variables.PlasticDissipation=0;
          if( thermo_mechanical ){
            //std::cout<<" thermo_mechanical element RHS ["<<this->Id()<<"]"<<std::endl;
            if( ConstitutiveLawVector[PointNumber]->Has(PLASTIC_DISSIPATION) ){
              ConstitutiveLawVector[PointNumber]->GetValue(PLASTIC_DISSIPATION,Variables.PlasticDissipation);
              //std::cout<<" with PLASTIC_DISSIPATION= "<<Variables.PlasticDissipation<<std::endl;
            }
          }

          this->CalculateAndAddRHS ( rRightHandSideVector, Variables, HeatSource, IntegrationWeight );
        }



    }

    // if( Variables.PlasticDissipation!=0 || Variables.DeltaPlasticDissipation!=0 )
    //   {
    // 	std::cout<<" Thermal Element "<<this->Id()<<" [ dissipation:"<<Variables.PlasticDissipation<<", delta_dissipation:"<<Variables.DeltaPlasticDissipation<<"] "<<std::endl;
    // 	std::cout<<" K "<<rLeftHandSideMatrix<<std::endl;
    // 	std::cout<<" F "<<rRightHandSideVector<<std::endl;
    //   }


    KRATOS_CATCH( "" )
}

//***********************COMPUTE LOCAL SYSTEM CONTRIBUTIONS***************************
//************************************************************************************

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{
  //constant during the step
  this->CalculateAndAddKthermal( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  //constant during the step
  this->CalculateAndAddMthermal( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

  //changes during the step
  this->CalculateAndAddHthermal( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

}


//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, double& rHeatSource, double& rIntegrationWeight)
{

  // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
  this->CalculateAndAddExternalForces( rVariables, rRightHandSideVector, rHeatSource, rIntegrationWeight );

  // operation performed: rRightHandSideVector -= ThermalForce*IntegrationWeight
  this->CalculateAndAddThermalForces( rVariables, rRightHandSideVector, rIntegrationWeight );

}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ThermalElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( LeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}


//************************************************************************************
//************************************************************************************


void ThermalElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ThermalElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( rLeftHandSideMatrix, RightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}

//************************************************************************************
//************************************************************************************


void ThermalElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{

    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(ThermalElement::COMPUTE_LHS_MATRIX);
    CalculationFlags.Set(ThermalElement::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}


//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************


double& ThermalElement::CalculateIntegrationWeight( double& rIntegrationWeight )
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 ){
      if ( this->GetProperties().Has( THICKNESS ) )
        rIntegrationWeight *= GetProperties()[THICKNESS];
    }

    return rIntegrationWeight;
}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddExternalForces(GeneralVariables& rVariables,
                                                   VectorType& rRightHandSideVector,
                                                   double& rHeatSource,
                                                   double& rIntegrationWeight)
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    VectorType Fh=rRightHandSideVector;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      rRightHandSideVector[i] += rIntegrationWeight * rVariables.N[i] * rHeatSource;
    }

    // std::cout<<std::endl;
    // std::cout<<" Fext "<<rRightHandSideVector-Fh<<std::endl;


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

inline void ThermalElement::CalculateAndAddThermalForces(GeneralVariables & rVariables,
							 VectorType& rRightHandSideVector,
							 double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    VectorType Fh=rRightHandSideVector;


    double DeltaTemperature = 0;
    double NodalTemperature = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	DeltaTemperature = GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE) - GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1);

	rRightHandSideVector[i] += rVariables.N[i] * rVariables.PlasticDissipation * rIntegrationWeight;
	rRightHandSideVector[i] -= rVariables.N[i] * rVariables.HeatCapacity * ( DeltaTemperature ) * (1.0/rVariables.DeltaTime) * rIntegrationWeight;

	for ( unsigned int j = 0; j < number_of_nodes; j++ )
	{
	  NodalTemperature = GetGeometry()[j].FastGetSolutionStepValue(TEMPERATURE);

	  for ( unsigned int k= 0; k < dimension; k++ )
	    {

	      rRightHandSideVector[i] -= rVariables.HeatConductivity * NodalTemperature * ( rVariables.DN_DX(i,k)*rVariables.DN_DX(j,k) ) * rIntegrationWeight;

	    }


	}

    }

    // std::cout<<std::endl;
    // std::cout<<" Fint ["<<GetValue(MASTER_ELEMENTS).back().Id()<<"] "<<rRightHandSideVector-Fh<<" [ Dissipation: "<<rVariables.PlasticDissipation<<" , HeatConductivity: "<<rVariables.HeatConductivity<<" , HeatCapacity: "<<rVariables.HeatCapacity<<" ]"<<std::endl;

    // if( rVariables.PlasticDissipation > 0 ){
    //   std::cout<<" ThermalElement Id["<<this->Id()<<"] Dissipation "<<rVariables.PlasticDissipation<<" Fint "<<rRightHandSideVector-Fh<<std::endl;
    //   for ( unsigned int i = 0; i < number_of_nodes; i++ )
    // 	{
    // 	  std::cout<<" Node["<<GetGeometry()[i].Id()<<"] DT: "<<GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE) - GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1)<<" T "<<GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE)<<" Capacity "<<rVariables.HeatCapacity * ( GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE) - GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE,1) ) * (1.0/rVariables.DeltaTime)<<" Conductivity "<<rVariables.HeatConductivity<<std::endl;
    // 	}
    // }

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddHthermal(MatrixType& rH,
					     GeneralVariables& rVariables,
					     double& rIntegrationWeight)

{
    KRATOS_TRY

    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    MatrixType Hh=rH;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
        rH(i,j) += rVariables.DeltaPlasticDissipation * rIntegrationWeight * rVariables.N[i] * rVariables.N[j];
      }

    }

   // std::cout<<std::endl;
   // std::cout<<" Htherm "<<rH-Hh<<" [ DeltaDissipation : "<<rVariables.DeltaPlasticDissipation<<" ]"<<std::endl;

    KRATOS_CATCH( "" )
 }


//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddKthermal(MatrixType& rK,
					     GeneralVariables& rVariables,
					     double& rIntegrationWeight)

{
    KRATOS_TRY

    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //double HeatConductivity = GetProperties()[HEAT_CONDUCTIVITY];


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
        for ( unsigned int k = 0; k < dimension; k++ )
        {
          rK(i,j)+= rVariables.HeatConductivity * ( rVariables.DN_DX(i,k)*rVariables.DN_DX(j,k) ) * rIntegrationWeight;
        }
      }
    }

    // std::cout<<std::endl;
    // std::cout<<" Ktherm "<<rK<<" [ head conductivity : "<<rVariables.HeatConductivity<<" ] "<<std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateAndAddMthermal(MatrixType& rM,
					     GeneralVariables& rVariables,
					     double& rIntegrationWeight)

{
    KRATOS_TRY

    //assemble into rk the material uu contribution:
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //double HeatCapacity = GetProperties()[HEAT_CAPACITY] * GetProperties()[DENSITY];

    double consistent = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
      for ( unsigned int j = 0; j < number_of_nodes; j++ )
      {
        consistent = 0;
        if(i==j){
          if( dimension == 2)
            consistent = 3;
          else
            consistent = 4;
        }

        rM(i,j)+= consistent * rVariables.HeatCapacity * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight * (1.0/rVariables.DeltaTime);
      }

    }

    // std::cout<<std::endl;
    // std::cout<<" Mtherm "<<rM<<" [ heatcapacity : "<<rVariables.HeatCapacity<<" ]"<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void ThermalElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  ThermalElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY

     // Check that all required variables have been registered
     KRATOS_CHECK_VARIABLE_KEY(TEMPERATURE);
     
     KRATOS_CHECK_VARIABLE_KEY(HEAT_CAPACITY);
     KRATOS_CHECK_VARIABLE_KEY(HEAT_CONDUCTIVITY);
     
     KRATOS_CHECK_VARIABLE_KEY(HEAT_SOURCE);
     KRATOS_CHECK_VARIABLE_KEY(PLASTIC_DISSIPATION);
     KRATOS_CHECK_VARIABLE_KEY(DELTA_PLASTIC_DISSIPATION);
     
     return 0;

     KRATOS_CATCH( "" );
}


void ThermalElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);

}

void ThermalElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);

}



} // Namespace Kratos
