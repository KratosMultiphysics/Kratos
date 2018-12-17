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
#include "custom_elements/solid_elements/solid_element.hpp"
#include "solid_mechanics_application_variables.h"


namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SolidElement, COMPUTE_RHS_VECTOR,                 0 );
KRATOS_CREATE_LOCAL_FLAG( SolidElement, COMPUTE_LHS_MATRIX,                 1 );
KRATOS_CREATE_LOCAL_FLAG( SolidElement, FINALIZED_STEP,                     2 );


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SolidElement::SolidElement( )
    :Element( )
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SolidElement::SolidElement( IndexType NewId, GeometryType::Pointer pGeometry )
    :Element( NewId, pGeometry )
{
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SolidElement::SolidElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    :Element( NewId, pGeometry, pProperties )
{
    this->Set(SOLID);
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SolidElement::SolidElement( SolidElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SolidElement&  SolidElement::operator=(SolidElement const& rOther)
{
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size() );

    for(SizeType i=0; i<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SolidElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    KRATOS_ERROR << " calling the default method Create for a solid element " << std::endl;
    return Kratos::make_shared< SolidElement >(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SolidElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    KRATOS_ERROR << " calling the default method Clone for a solid element " << std::endl;

    SolidElement NewElement( NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );


    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;


    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());

	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR << " constitutive law not has the correct size solid element " << std::endl;

      }

    for(SizeType i=0; i<mConstitutiveLawVector.size(); i++)
      {
	NewElement.mConstitutiveLawVector[i] = mConstitutiveLawVector[i]->Clone();
      }

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Kratos::make_shared< SolidElement >(NewElement);
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SolidElement::~SolidElement()
{
}


//************************************************************************************
//************************************************************************************

SolidElement::IntegrationMethod SolidElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}


//************************************************************************************
//************************************************************************************

void SolidElement::IncreaseIntegrationMethod(IntegrationMethod& rThisIntegrationMethod, unsigned int increment) const
{
  int IntMethod = int(rThisIntegrationMethod);
  IntMethod += increment;
  rThisIntegrationMethod = IntegrationMethod(IntMethod);
}

//************************************************************************************
//************************************************************************************

void SolidElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    for ( SizeType i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}


//************************************************************************************
//************************************************************************************

void SolidElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = GetDofsSize();

    if ( rResult.size() != dofs_size )
        rResult.resize( dofs_size, false );

    SizeType index = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void SolidElement::GetValuesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    SizeType index = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );

    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void SolidElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    SizeType index = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void SolidElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType dofs_size = GetDofsSize();

    if ( rValues.size() != dofs_size )
      rValues.resize( dofs_size, false );

    SizeType index = 0;
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

}

//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void SolidElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void SolidElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }

}

//*********************************SET MATRIX VALUE***********************************
//************************************************************************************


void SolidElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( SizeType PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }

}


//********************************SET CONSTITUTIVE VALUE******************************
//************************************************************************************

void SolidElement::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
  if(rVariable == CONSTITUTIVE_LAW)  //returns a vector of pointers, it do not clones the constitutives laws (it is not a copy)
    {
        if ( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize(rValues.size());

            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod ) )
                KRATOS_ERROR << " constitutive law not has the correct size solid element " << std::endl;
        }

        for(SizeType i=0; i<rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i];
        }
    }

}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************


void SolidElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
						std::vector<double>& rValues,
						const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == VON_MISES_STRESS || rVariable == NORM_ISOCHORIC_STRESS || rVariable == PRESSURE || DAMAGE_VARIABLE)
    {
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
    }
    else
    {
      const SizeType& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

      if ( rValues.size() != integration_points_number )
      {
        rValues.resize( integration_points_number );
      }

      for ( SizeType ii = 0; ii < integration_points_number; ii++ )
      {
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
      }
    }
}

//**********************************GET VECTOR VALUE**********************************
//************************************************************************************


void SolidElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
						std::vector<Vector>& rValues,
						const ProcessInfo& rCurrentProcessInfo )
{
    const SizeType& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );


    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else if ( rVariable == PK2_STRESS_VECTOR ||  rVariable == CAUCHY_STRESS_VECTOR )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else
    {

        for ( SizeType PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }

    }

}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void SolidElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
						std::vector<Matrix>& rValues,
						const ProcessInfo& rCurrentProcessInfo )
{

    const SizeType& integration_points_number = mConstitutiveLawVector.size();

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    if ( rVariable == PK2_STRESS_TENSOR ||  rVariable == CAUCHY_STRESS_TENSOR )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR ||  rVariable == ALMANSI_STRAIN_TENSOR )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else if ( rVariable == DEFORMATION_GRADIENT )
    {

        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    }
    else
    {

        for ( SizeType PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }

    }


}

//********************************GET CONSTITUTIVE VALUE******************************
//************************************************************************************

void SolidElement::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
						std::vector<ConstitutiveLaw::Pointer>& rValues,
						const ProcessInfo& rCurrentProcessInfo )
{

    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            rValues.resize(mConstitutiveLawVector.size());
        }

        for(SizeType i=0; i<rValues.size(); i++)
        {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }

}


//************************************************************************************
//************************************************************************************

void SolidElement::Initialize()
{
    KRATOS_TRY

    //member variables initialization
    InitializeConstitutiveLaw();

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::InitializeConstitutiveLaw()
{
    KRATOS_TRY


    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
    }

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( SizeType i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else{

      KRATOS_ERROR << " a constitutive law needs to be specified for the element with ID " << std::endl;
    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( SizeType i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::InitializeElementData (ElementDataType& rVariables, const ProcessInfo& rCurrentProcessInfo)
{

    const SizeType number_of_nodes  = GetGeometry().size();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
    const SizeType voigt_size   = dimension * (dimension +1) * 0.5;

    rVariables.Initialize(voigt_size,dimension,number_of_nodes);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //set process info
    rVariables.SetProcessInfo(rCurrentProcessInfo);

    //calculating the current jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n+1/d£]
    rVariables.j = GetGeometry().Jacobian( rVariables.j, mThisIntegrationMethod );

}

//*********************************COMPUTE KINETICS***********************************
//************************************************************************************

void SolidElement::CalculateKinetics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    this->CalculateKinematics(rVariables,rPointNumber);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::TransformElementData(ElementDataType& rVariables, const double& rPointNumber)
{
  // to be used when different reference configuration is used
}

//************************************************************************************
//************************************************************************************

SolidElement::SizeType SolidElement::GetDofsSize()
{
  KRATOS_TRY

  const SizeType dimension        = GetGeometry().WorkingSpaceDimension();
  const SizeType number_of_nodes = GetGeometry().PointsNumber();

  SizeType size = number_of_nodes * dimension; //usual size for displacement based elements

  return size;

  KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
					    VectorType& rRightHandSideVector,
					    Flags& rCalculationFlags)

{

    //resizing as needed the LHS
    const SizeType MatSize = this->GetDofsSize();

    if ( rCalculationFlags.Is(SolidElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(SolidElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

    }
}

//************************************************************************************
//************************************************************************************

void SolidElement::SetElementData(ElementDataType& rVariables,
				       ConstitutiveLaw::Parameters& rValues,
				       const int & rPointNumber)
{

}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateMaterialResponse(ElementDataType& rVariables,
                                             ConstitutiveLaw::Parameters& rValues,
                                             const int & rPointNumber)
{
    KRATOS_TRY

    //set general variables to constitutivelaw parameters
    this->SetElementData(rVariables,rValues,rPointNumber);

    //compute stresses and constitutive parameters
    mConstitutiveLawVector[rPointNumber]->CalculateMaterialResponse(rValues, rVariables.StressMeasure);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
                                             ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
    Vector VolumeForce(dimension);
    noalias(VolumeForce) = ZeroVector(dimension);

    for ( SizeType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematic variables B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //calculate material response
        this->CalculateMaterialResponse(Variables,Values,PointNumber);

	//some transformation of the configuration can be needed (UL element specially)
        this->TransformElementData(Variables,PointNumber);

        //calculating weights for integration on the "reference configuration"
        Variables.IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
        Variables.IntegrationWeight = this->CalculateIntegrationWeight( Variables.IntegrationWeight );


        if( !IsSliver() ){

          if ( rLocalSystem.CalculationFlags.Is(SolidElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
          {
            //contributions to stiffness matrix calculated on the reference config
	    this->CalculateAndAddLHS ( rLocalSystem, Variables, Variables.IntegrationWeight );
          }

          if ( rLocalSystem.CalculationFlags.Is(SolidElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
          {
            //contribution to external forces
            VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );

	    this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, Variables.IntegrationWeight );
          }

        }

	//for debugging purposes
	//this->PrintElementCalculation(rLocalSystem,Variables);

    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateDynamicSystem( LocalSystemComponents& rLocalSystem,
					   ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    IntegrationMethod ThisIntegrationMethod = mThisIntegrationMethod;

    if( rCurrentProcessInfo.Has(COMPUTE_CONSISTENT_MASS_MATRIX) ){

      if(rCurrentProcessInfo[COMPUTE_CONSISTENT_MASS_MATRIX] == true){
	//full quadrature integration:
	this->IncreaseIntegrationMethod(mThisIntegrationMethod,1);
      }
    }

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);


    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    for ( SizeType PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinetic variables B, F, DN_DX ...
        this->CalculateKinetics(Variables, PointNumber);

        //calculating weights for integration on the "reference configuration"
        Variables.IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
        Variables.IntegrationWeight = this->CalculateIntegrationWeight( Variables.IntegrationWeight );


        if ( rLocalSystem.CalculationFlags.Is(SolidElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
	  MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
	  this->CalculateAndAddDynamicLHS ( rLeftHandSideMatrix, Variables, rCurrentProcessInfo, Variables.IntegrationWeight );

        }

        if ( rLocalSystem.CalculationFlags.Is(SolidElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
	  VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();
	  this->CalculateAndAddDynamicRHS ( rRightHandSideVector, Variables, rCurrentProcessInfo, Variables.IntegrationWeight );
        }

	//for debugging purposes
	//this->PrintElementCalculation(rLocalSystem,Variables);

    }

    mThisIntegrationMethod = ThisIntegrationMethod;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void SolidElement::PrintElementCalculation(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables)
{
    KRATOS_TRY

    std::cout<<" Element: "<<this->Id()<<std::endl;
    SizeType number_of_nodes = GetGeometry().PointsNumber();
    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	array_1d<double, 3 > & CurrentPosition      = GetGeometry()[i].Coordinates();
	array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);
	std::cout<<" Previous  Position  node["<<GetGeometry()[i].Id()<<"]: "<<PreviousPosition<<std::endl;
      }

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	array_1d<double, 3> & CurrentPosition  = GetGeometry()[i].Coordinates();
	std::cout<<" Current  Position  node["<<GetGeometry()[i].Id()<<"]: "<<CurrentPosition<<std::endl;
      }

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
	std::cout<<" Previous Displacement  node["<<GetGeometry()[i].Id()<<"]: "<<PreviousDisplacement<<std::endl;
      }

    for ( SizeType i = 0; i < number_of_nodes; i++ )
      {
	array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
	std::cout<<" Current  Displacement  node["<<GetGeometry()[i].Id()<<"]: "<<CurrentDisplacement<<std::endl;
      }

    std::cout<<" Stress "<<rVariables.StressVector<<std::endl;
    std::cout<<" Strain "<<rVariables.StrainVector<<std::endl;
    std::cout<<" F  "<<rVariables.F<<" detF "<<rVariables.detF<<std::endl;
    std::cout<<" F0 "<<rVariables.F0<<" detF0 "<<rVariables.detF0<<std::endl;
    std::cout<<" ConstitutiveMatrix "<<rVariables.ConstitutiveMatrix<<std::endl;
    std::cout<<" K "<<rLocalSystem.GetLeftHandSideMatrix()<<std::endl;
    std::cout<<" f "<<rLocalSystem.GetRightHandSideVector()<<std::endl;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, double& rIntegrationWeight)
{
    KRATOS_TRY

    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();

    // operation performed: add Km to the rLefsHandSideMatrix
    this->CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // operation performed: add Kg to the rLefsHandSideMatrix
    this->CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, ElementDataType& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{
    KRATOS_TRY

    //contribution of the internal and external forces

    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector();

    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
    this->CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddDynamicLHS(MatrixType& rLeftHandSideMatrix, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
  KRATOS_TRY

  //mass matrix
  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  const SizeType dimension = GetGeometry().WorkingSpaceDimension();
  const SizeType MatSize = this->GetDofsSize();

  if(rLeftHandSideMatrix.size1() != MatSize)
    rLeftHandSideMatrix.resize (MatSize, MatSize, false);

  //compute volume change
  double VolumeChange = 1;
  VolumeChange = this->CalculateVolumeChange( VolumeChange, rVariables );

  double CurrentDensity = VolumeChange * GetProperties()[DENSITY];


  unsigned int indexi = 0;
  unsigned int indexj = 0;

  for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      for ( SizeType k = 0; k < dimension; k++ )
	{
	  indexj = 0;
	  for ( SizeType j = 0; j < number_of_nodes; j++ )
	    {

	      rLeftHandSideMatrix(indexi+k,indexj+k) += rVariables.N[i] * rVariables.N[j] * CurrentDensity * rIntegrationWeight;
	      indexj += dimension;
	    }

	}

      indexi += dimension;
    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddDynamicRHS(VectorType& rRightHandSideVector, ElementDataType& rVariables, ProcessInfo& rCurrentProcessInfo, double& rIntegrationWeight)
{
  KRATOS_TRY

  //mass matrix
  const SizeType number_of_nodes = GetGeometry().PointsNumber();
  const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
  const unsigned int MatSize = this->GetDofsSize();

  MatrixType MassMatrix( MatSize, MatSize );
  noalias(MassMatrix) = ZeroMatrix( MatSize, MatSize );

  //compute volume change
  double VolumeChange = 1;
  VolumeChange = this->CalculateVolumeChange( VolumeChange, rVariables );

  double CurrentDensity = VolumeChange * GetProperties()[DENSITY];

  //acceleration vector
  Vector CurrentAccelerationVector(MatSize);
  noalias(CurrentAccelerationVector) = ZeroVector( MatSize );
  this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);

  double AlphaM = 0.0;
  if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
    AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
    Vector PreviousAccelerationVector(MatSize);
    noalias(PreviousAccelerationVector) = ZeroVector( MatSize );
    this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
    CurrentAccelerationVector *= (1.0-AlphaM);
    CurrentAccelerationVector +=  AlphaM * (PreviousAccelerationVector);
  }

  unsigned int indexi = 0;
  unsigned int indexj = 0;

  for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      for ( SizeType k = 0; k < dimension; k++ )
	{
	  indexj = 0;
	  for ( SizeType j = 0; j < number_of_nodes; j++ )
	    {

	      MassMatrix(indexi+k,indexj+k) += rVariables.N[i] * rVariables.N[j] * CurrentDensity * rIntegrationWeight;
	      indexj += dimension;
	    }

	}

      indexi += dimension;
    }


  noalias(rRightHandSideVector) = prod( MassMatrix, CurrentAccelerationVector );

  KRATOS_CATCH("")
}


//***********************************************************************************
//************************************************************************************

double& SolidElement::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 ){
      if ( this->GetProperties().Has( THICKNESS ) )
        rIntegrationWeight *= GetProperties()[THICKNESS];
    }

    return rIntegrationWeight;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************


void SolidElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(rRightHandSideVector);

    //Calculate elemental system
    this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );

    bool test_tangent = false;
    if( test_tangent ){

      //std::cout<<" ["<<this->Id()<<"] MATRIX "<<rLeftHandSideMatrix<<std::endl;

      MatrixType PerturbedLeftHandSideMatrix( rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2() );
      noalias(PerturbedLeftHandSideMatrix) = ZeroMatrix( rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2() );

      this->CalculatePerturbedLeftHandSide( PerturbedLeftHandSideMatrix, rCurrentProcessInfo );

      //std::cout<<" ["<<this->Id()<<"] PERTURBED MATRIX "<<PerturbedLeftHandSideMatrix<<std::endl;

      //std::cout<<" ["<<this->Id()<<"] DIFFERENCES "<<PerturbedLeftHandSideMatrix-rLeftHandSideMatrix<<std::endl;

      //rLeftHandSideMatrix = PerturbedLeftHandSideMatrix;

    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************


void SolidElement::CalculatePerturbedLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create local system components
    LocalSystemComponents LocalSystem;

    //calculation flags
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);
    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_RHS_VECTOR);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, LocalSystem.CalculationFlags );

    //Set Variables to Local system components
    LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
    LocalSystem.SetRightHandSideVector(RightHandSideVector);

    LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX,false);

    DofsVectorType ElementalDofList;
    this->GetDofList( ElementalDofList, rCurrentProcessInfo );

    unsigned int size = ElementalDofList.size();
    for( unsigned int i=0; i<size; i++)
      {
	//Set perturbation in "i" dof component
	double& value    = ElementalDofList[i]->GetSolutionStepValue();
	double  original = value;

	double  deltavalue = 1e-10;
	if( value !=0 )
	  deltavalue = value * 1e-8;

	//Calculate elemental system
	RightHandSideVector.resize(RightHandSideVector.size(),false);
	noalias(RightHandSideVector) = ZeroVector(RightHandSideVector.size());
	value = original + deltavalue;
	this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
	VectorType RightHandSideVectorI = RightHandSideVector;

	//Calculate elemental system
	RightHandSideVector.resize(RightHandSideVector.size(),false);
	noalias(RightHandSideVector) = ZeroVector(RightHandSideVector.size());
	value = original - deltavalue;
	this->CalculateElementalSystem( LocalSystem, rCurrentProcessInfo );
	VectorType RightHandSideVectorII = RightHandSideVector;

	// std::cout<<" i: "<<i<<" RHS "<<RightHandSideVectorI<<std::endl;
	// std::cout<<" ii: "<<i<<" RHS "<<RightHandSideVectorII<<std::endl;

	for( unsigned int j=0; j<size; j++)
	  {
	    rLeftHandSideMatrix(j,i) = (-1) * (RightHandSideVectorI[j] - RightHandSideVectorII[j]) / (2.0*deltavalue);
	  }

	value = original;
      }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    InitializeExplicitContributions();

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                rCurrentProcessInfo );

    this->Set(SolidElement::FINALIZED_STEP,false);

    KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************
void SolidElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    InitializeExplicitContributions();
}

//************************************************************************************
//************************************************************************************

void SolidElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

//************************************************************************************
//************************************************************************************

void SolidElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //create and initialize element variables:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {

        //compute element kinematic variables B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetElementData(Variables,Values,PointNumber);

        //call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponse(Values, Variables.StressMeasure);

        //call the constitutive law to finalize the solution step
        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
								   GetGeometry(),
								   Variables.N,
								   rCurrentProcessInfo );

	//call the element internal variables update
	this->FinalizeStepVariables(Variables,PointNumber);
    }

    this->Set(SolidElement::FINALIZED_STEP,true);

    if(this->Is(SELECTED) && this->Is(ACTIVE)){
      this->Set(SELECTED,false);
      std::cout<<" Undo SELECTED "<<this->Id()<<std::endl;
    }

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void SolidElement::FinalizeStepVariables( ElementDataType & rVariables, const double& rPointNumber )
{
  //to update the internal element variables
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
				       ElementDataType& rVariables,
				       double& rIntegrationWeight)
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    noalias( rLeftHandSideMatrix ) += rIntegrationWeight * prod( trans( rVariables.B ), Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary

    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
				       ElementDataType& rVariables,
				       double& rIntegrationWeight)

{

}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
						 ElementDataType& rVariables,
						 Vector& rVolumeForce,
						 double& rIntegrationWeight)

{
    KRATOS_TRY
    SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;
        for ( SizeType j = 0; j < dimension; j++ )
        {
	  rRightHandSideVector[index + j] += rIntegrationWeight * rVariables.N[i] * rVolumeForce[j];
        }
    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        ElementDataType & rVariables,
        double& rIntegrationWeight
                                                            )
{
    KRATOS_TRY

    VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
    noalias( rRightHandSideVector ) -= InternalForces;

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::InitializeExplicitContributions()
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    for ( SizeType i = 0; i < number_of_nodes; i++ )
    {
      if( GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_FORCE) && GetGeometry()[i].SolutionStepsDataHas(INTERNAL_FORCE) ){

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
        array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);

    	GetGeometry()[i].SetLock();
        ExternalForce.clear();
        InternalForce.clear();
    	GetGeometry()[i].UnSetLock();

      }
    }

    KRATOS_CATCH("")
}

//***********************************************************************************
//***********************************************************************************

void SolidElement::AddExplicitContribution(const VectorType& rRHSVector,
                                           const Variable<VectorType>& rRHSVariable,
                                           Variable<array_1d<double,3> >& rDestinationVariable,
                                           const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
      {

	for(SizeType i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(SizeType j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
      {

	for(SizeType i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &InternalForce = GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
	    for(SizeType j=0; j<dimension; j++)
	      {
		InternalForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(SizeType i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

	    for(SizeType j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH("")
}


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void SolidElement::CalculateKinematics(ElementDataType& rVariables, const double& rPointNumber)
{
    KRATOS_TRY

    KRATOS_ERROR << " calling the default method CalculateKinematics for a solid element " << std::endl;

    KRATOS_CATCH("")
}



//*************************COMPUTE DELTA POSITION*************************************
//************************************************************************************


Matrix& SolidElement::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const GeometryType& rGeometry = GetGeometry();

    ElementUtilities::CalculateDeltaPosition(rDeltaPosition,rGeometry);

    return rDeltaPosition;

    KRATOS_CATCH("")
}

//*************************COMPUTE TOTAL DELTA POSITION*******************************
//************************************************************************************

Matrix& SolidElement::CalculateTotalDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const GeometryType& rGeometry = GetGeometry();

    ElementUtilities::CalculateTotalDeltaPosition(rDeltaPosition,rGeometry);
    return rDeltaPosition;

    KRATOS_CATCH("")
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& SolidElement::CalculateTotalMass( double& rTotalMass, const ProcessInfo& rCurrentProcessInfo  )
{
    KRATOS_TRY

    const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

    //Compute the Volume Change acumulated:
    ElementDataType Variables;
    this->InitializeElementData(Variables,rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    double VolumeChange = 1;

    //reading integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
      {
	//compute element kinematic variables
	this->CalculateKinematics(Variables,PointNumber);

	//getting informations for integration
        Variables.IntegrationWeight = Variables.detJ * integration_points[PointNumber].Weight();

	//compute point volume changes
	VolumeChange = 1;
	VolumeChange = this->CalculateVolumeChange( VolumeChange, Variables );

	rTotalMass += VolumeChange * GetProperties()[DENSITY] * Variables.IntegrationWeight;

      }

    if( dimension == 2 ){
      if ( this->GetProperties().Has( THICKNESS ) )
	rTotalMass *= GetProperties()[THICKNESS];
    }


    return rTotalMass;

    KRATOS_CATCH("")
}

//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

double& SolidElement::CalculateVolumeChange( double& rVolumeChange, ElementDataType& rVariables )
{
    KRATOS_TRY

    rVolumeChange = 1; //no change (reference configuration considered by default)

    return rVolumeChange;

    KRATOS_CATCH("")
}


//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& SolidElement::CalculateVolumeForce( Vector& rVolumeForce, ElementDataType& rVariables )
{
    KRATOS_TRY

    const SizeType number_of_nodes = GetGeometry().PointsNumber();
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    if(rVolumeForce.size() != dimension)
      rVolumeForce.resize(dimension,false);

    noalias(rVolumeForce) = ZeroVector(dimension);

    for ( SizeType j = 0; j < number_of_nodes; j++ )
    {
      if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ){ // it must be checked once at the begining only
	array_1d<double, 3 >& VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
	for( SizeType i = 0; i < dimension; i++ )
	  rVolumeForce[i] += rVariables.N[j] * VolumeAcceleration[i] ;
      }
    }

    double VolumeChange = 1;
    VolumeChange = this->CalculateVolumeChange( VolumeChange, rVariables );

    rVolumeForce *= VolumeChange * GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //1.-Calculate Tangent Inertia Matrix:
    this->CalculateDampingMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

    const unsigned int MatSize = this->GetDofsSize();

    //2.-Calculate Inertial Forces:
    if ( rRightHandSideVector.size() != MatSize )
      rRightHandSideVector.resize( MatSize, false );

    noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

    //acceleration vector
    Vector CurrentVelocityVector(MatSize);
    noalias(CurrentVelocityVector) = ZeroVector( MatSize );
    this->GetFirstDerivativesVector(CurrentVelocityVector, 0);

    noalias(rRightHandSideVector) = prod( rLeftHandSideMatrix, CurrentVelocityVector );


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) ){

      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true){
	ComputeDynamicTangent = true;
      }
    }

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_RHS_VECTOR);
      LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

      const unsigned int MatSize = this->GetDofsSize();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

      //acceleration vector
      Vector CurrentAccelerationVector(MatSize);
      noalias(CurrentAccelerationVector) = ZeroVector( MatSize );
      this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);

      double AlphaM = 0.0;
      if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
	AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
	Vector PreviousAccelerationVector( MatSize );
	noalias(PreviousAccelerationVector) = ZeroVector( MatSize );
	this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
	CurrentAccelerationVector *= (1.0-AlphaM);
	CurrentAccelerationVector +=  AlphaM * (PreviousAccelerationVector);
      }

      noalias(rRightHandSideVector) = prod( rLeftHandSideMatrix, CurrentAccelerationVector );

    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);

      VectorType RightHandSideVector = Vector();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rLeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(RightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( rLeftHandSideMatrix, rCurrentProcessInfo );

    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool ComputeDynamicTangent = false;
    if( rCurrentProcessInfo.Has(COMPUTE_DYNAMIC_TANGENT) )
      if(rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true)
	ComputeDynamicTangent = true;

    if( ComputeDynamicTangent == true ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_RHS_VECTOR);

      MatrixType LeftHandSideMatrix = Matrix();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(LeftHandSideMatrix);
      LocalSystem.SetRightHandSideVector(rRightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      MatrixType LeftHandSideMatrix = Matrix();

      //1.-Calculate Tangent Inertia Matrix:
      this->CalculateMassMatrix( LeftHandSideMatrix, rCurrentProcessInfo );

      const unsigned int MatSize = this->GetDofsSize();

      //2.-Calculate Inertial Forces:
      if ( rRightHandSideVector.size() != MatSize )
	rRightHandSideVector.resize( MatSize, false );

      noalias(rRightHandSideVector) = ZeroVector( MatSize ); //resetting RHS

      //acceleration vector
      Vector CurrentAccelerationVector   = ZeroVector( MatSize );
      this->GetSecondDerivativesVector(CurrentAccelerationVector, 0);

      double AlphaM = 0.0;
      if( rCurrentProcessInfo.Has(BOSSAK_ALPHA) ){
	AlphaM = rCurrentProcessInfo[BOSSAK_ALPHA];
	Vector PreviousAccelerationVector( MatSize );
	noalias(PreviousAccelerationVector) = ZeroVector( MatSize );
	this->GetSecondDerivativesVector(PreviousAccelerationVector, 1);
	CurrentAccelerationVector *= (1.0-AlphaM);
	CurrentAccelerationVector +=  AlphaM * (PreviousAccelerationVector);
      }

      noalias(rRightHandSideVector) = prod( LeftHandSideMatrix, CurrentAccelerationVector );

    }


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    bool ComputeLumpedMassMatrix = false;
    if( rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX) )
      if(rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == true)
	ComputeLumpedMassMatrix = true;


    if( ComputeLumpedMassMatrix == false ){

      //create local system components
      LocalSystemComponents LocalSystem;

      //calculation flags
      LocalSystem.CalculationFlags.Set(SolidElement::COMPUTE_LHS_MATRIX);

      VectorType RightHandSideVector = Vector();

      //Initialize sizes for the system components:
      this->InitializeSystemMatrices( rMassMatrix, RightHandSideVector,  LocalSystem.CalculationFlags );

      //Set Variables to Local system components
      LocalSystem.SetLeftHandSideMatrix(rMassMatrix);
      LocalSystem.SetRightHandSideVector(RightHandSideVector);

      //Calculate elemental system
      CalculateDynamicSystem( LocalSystem, rCurrentProcessInfo );

    }
    else{

      //lumped
      const SizeType dimension  = GetGeometry().WorkingSpaceDimension();
      const SizeType number_of_nodes = GetGeometry().PointsNumber();
      const unsigned int MatSize = this->GetDofsSize();
      if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

      noalias(rMassMatrix) = ZeroMatrix( MatSize, MatSize );

      double TotalMass = 0;
      TotalMass = this->CalculateTotalMass(TotalMass,rCurrentProcessInfo);

      Vector LumpFact(number_of_nodes);
      noalias(LumpFact) = ZeroVector(number_of_nodes);

      LumpFact  = GetGeometry().LumpingFactors( LumpFact );

      for ( SizeType i = 0; i < number_of_nodes; i++ )
	{
	  double temp = LumpFact[i] * TotalMass;

	  for ( SizeType j = 0; j < dimension; j++ )
	    {
	      unsigned int index = i * dimension + j;
	      rMassMatrix( index, index ) = temp;
	    }
	}

    }


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //0.-Initialize the DampingMatrix:

    //resizing as needed the LHS
    const unsigned int MatSize = this->GetDofsSize();

    if ( rDampingMatrix.size1() != MatSize )
        rDampingMatrix.resize( MatSize, MatSize, false );

    noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


    //1.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
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

    if( alpha != 0 || beta != 0){

      //1.-Calculate StiffnessMatrix:

      MatrixType StiffnessMatrix  = Matrix();

      this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

      //2.-Calculate MassMatrix:

      MatrixType MassMatrix  = Matrix();

      this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );

      //3.-Compose the Damping Matrix:

      //Rayleigh Damping Matrix: alpha*M + beta*K
      rDampingMatrix  = alpha * MassMatrix;
      rDampingMatrix += beta  * StiffnessMatrix;

    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void SolidElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );


    if ( rVariable == DAMAGE_VARIABLE)
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateValue(Values,rVariable,rOutput[PointNumber]);

        }
    }
    if ( rVariable == VON_MISES_STRESS )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

            //compute element kinematic variables B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            rOutput[PointNumber] =  ElementUtilities::CalculateVonMises(Variables.StressVector);
        }
    }
    else if ( rVariable == NORM_ISOCHORIC_STRESS )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::ISOCHORIC_TENSOR_ONLY);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            rOutput[PointNumber] =  ElementUtilities::CalculateStressNorm(Variables.StressVector);
        }

    }
    else if (rVariable == PRESSURE)
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

	const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
	    this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

	    //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	    if( dimension == 2)
            {
                rOutput[PointNumber] = 0.5 * (Variables.StressVector[0] + Variables.StressVector[1]);
            }
            else
            {
                rOutput[PointNumber] = 1.0/3.0 * (Variables.StressVector[0] + Variables.StressVector[1] + Variables.StressVector[2]);
            }
        }
    }
    else if (rVariable == STRAIN_ENERGY)
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
	ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN_ENERGY);

	//reading integration points
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

	double StrainEnergy = 0;

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
	    this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

	    //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

	    //get strain energy
	    StrainEnergy = 0;
	    mConstitutiveLawVector[PointNumber]->GetValue(STRAIN_ENERGY,StrainEnergy);

	    Variables.IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
	    Variables.IntegrationWeight = this->CalculateIntegrationWeight( Variables.IntegrationWeight );

	    rOutput[PointNumber] = StrainEnergy * Variables.IntegrationWeight;
        }
    }
    else
    {

        for ( unsigned int ii = 0; ii < integration_points_number; ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
	  {
            //compute element kinematic variables B, F, DN_DX ...
	    this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

	    //call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR)
              mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);
            else
	      mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);

	    if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
                rOutput[PointNumber].resize( Variables.StressVector.size(), false );

            rOutput[PointNumber] = Variables.StressVector;
        }

    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const SizeType dimension        = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );


    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR )
    {
        std::vector<Vector> StressVector;
        if( rVariable == CAUCHY_STRESS_TENSOR )
            this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);

        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR )
    {

        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            this->CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else
            this->CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else if ( rVariable == CONSTITUTIVE_MATRIX )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetElementData(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
	    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

            if( rOutput[PointNumber].size2() != Variables.ConstitutiveMatrix.size2() )
                rOutput[PointNumber].resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );

            rOutput[PointNumber] = Variables.ConstitutiveMatrix;

        }

    }
    else if ( rVariable == DEFORMATION_GRADIENT )
    {
        //create and initialize element variables:
        ElementDataType Variables;
        this->InitializeElementData(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematic variables B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            if( rOutput[PointNumber].size2() != Variables.F.size2() )
                rOutput[PointNumber].resize( Variables.F.size1() , Variables.F.size2() , false );

            rOutput[PointNumber] = Variables.F;

        }
    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }


    KRATOS_CATCH("")
}



//************************************************************************************
//************************************************************************************

int  SolidElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    // Perform base element checks
    int ErrorCode = 0;
    ErrorCode = Element::Check(rCurrentProcessInfo);

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);

    KRATOS_CHECK_VARIABLE_KEY(DENSITY);
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION);

    KRATOS_CHECK_VARIABLE_KEY(VON_MISES_STRESS);
    KRATOS_CHECK_VARIABLE_KEY(NORM_ISOCHORIC_STRESS);
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_TENSOR);
    KRATOS_CHECK_VARIABLE_KEY(CAUCHY_STRESS_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_TENSOR);
    KRATOS_CHECK_VARIABLE_KEY(PK2_STRESS_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_TENSOR);
    KRATOS_CHECK_VARIABLE_KEY(GREEN_LAGRANGE_STRAIN_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_TENSOR);
    KRATOS_CHECK_VARIABLE_KEY(ALMANSI_STRAIN_VECTOR);
    KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_MATRIX);
    KRATOS_CHECK_VARIABLE_KEY(DEFORMATION_GRADIENT);
    KRATOS_CHECK_VARIABLE_KEY(STRAIN_ENERGY);


    // Check that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
      KRATOS_ERROR << "constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else{
      const SizeType dimension  = GetGeometry().WorkingSpaceDimension();

      if( dimension == 3 &&  this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
	KRATOS_ERROR <<  "wrong constitutive law used. This is a 3D element. Expected strain size is 6 :: element id " << this->Id() << std::endl;

      // Check constitutive law
      this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );

    }

    return ErrorCode;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void SolidElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
}

void SolidElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
}


} // Namespace Kratos
