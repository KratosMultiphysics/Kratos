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
#include "includes/define.h"
#include "custom_elements/small_displacement_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "solid_mechanics_application.h"


namespace Kratos
{

/**
 * Flags related to the element computation
 */
KRATOS_CREATE_LOCAL_FLAG( SmallDisplacementElement, COMPUTE_RHS_VECTOR,       0 );
KRATOS_CREATE_LOCAL_FLAG( SmallDisplacementElement, COMPUTE_LHS_MATRIX,       1 );


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

SmallDisplacementElement::SmallDisplacementElement( SmallDisplacementElement const& rOther)
    :Element(rOther)
    ,mThisIntegrationMethod(rOther.mThisIntegrationMethod)
    ,mConstitutiveLawVector(rOther.mConstitutiveLawVector)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

SmallDisplacementElement&  SmallDisplacementElement::operator=(SmallDisplacementElement const& rOther)
{
    Element::operator=(rOther);

    mThisIntegrationMethod = rOther.mThisIntegrationMethod;

    mConstitutiveLawVector.clear();
    mConstitutiveLawVector.resize( rOther.mConstitutiveLawVector.size());

    for(unsigned int i=0; i<<mConstitutiveLawVector.size(); i++)
    {
        mConstitutiveLawVector[i] = rOther.mConstitutiveLawVector[i];
    }


    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer SmallDisplacementElement::Create( IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallDisplacementElement( NewId, GetGeometry().Create( rThisNodes ), pProperties ) );
}


//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer SmallDisplacementElement::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    SmallDisplacementElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    //-----------//

    NewElement.mThisIntegrationMethod = mThisIntegrationMethod;

    if ( NewElement.mConstitutiveLawVector.size() != mConstitutiveLawVector.size() )
      {
	NewElement.mConstitutiveLawVector.resize(mConstitutiveLawVector.size());
	
	if( NewElement.mConstitutiveLawVector.size() != NewElement.GetGeometry().IntegrationPointsNumber() )
	  KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", NewElement.mConstitutiveLawVector.size() );
      }
    
       
    return Element::Pointer( new SmallDisplacementElement(NewElement) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************

SmallDisplacementElement::~SmallDisplacementElement()
{
}


//************* GETTING METHODS
//************************************************************************************
//************************************************************************************

SmallDisplacementElement::IntegrationMethod SmallDisplacementElement::GetIntegrationMethod() const
{
    return mThisIntegrationMethod;
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rElementalDofList.resize( 0 );
    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//*********************************DISPLACEMENT***************************************
//************************************************************************************

void SmallDisplacementElement::GetValuesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        array_1d<double, 3>& Displacement = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT, Step );
        for ( unsigned int j = 0; j < Displacement.size(); j++ )
        {
            rValues[ index + j ] = Displacement[j];
        }
    }
}


//************************************VELOCITY****************************************
//************************************************************************************

void SmallDisplacementElement::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        array_1d<double, 3>& Velocity = GetGeometry()[i].GetSolutionStepValue( VELOCITY, Step );
        for ( unsigned int j = 0; j < Velocity.size(); j++ )
        {
            rValues[ index + j ] = Velocity[j];
        }
    }
}

//*********************************ACCELERATION***************************************
//************************************************************************************

void SmallDisplacementElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension;

    if ( rValues.size() != element_size ) rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        array_1d<double, 3>& Acceleration = GetGeometry()[i].GetSolutionStepValue( ACCELERATION, Step );
        for ( unsigned int j = 0; j < Acceleration.size(); j++ )
        {
            rValues[ index + j ] = Acceleration[j];
        }
    }

}


//*********************************SET DOUBLE VALUE***********************************
//************************************************************************************

void SmallDisplacementElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }
}

//*********************************SET VECTOR VALUE***********************************
//************************************************************************************

void SmallDisplacementElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }

}

//*********************************SET MATRIX VALUE***********************************
//************************************************************************************


void SmallDisplacementElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
    }

}

//*********************************GET DOUBLE VALUE***********************************
//************************************************************************************


void SmallDisplacementElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
        std::vector<double>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rVariable == VON_MISES_STRESS )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    if ( rValues.size() != integration_points_number )
        rValues.resize( integration_points_number );

    for ( unsigned int ii = 0; ii < integration_points_number; ii++ )
        rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
}

//********************************SET CONSTITUTIVE VALUE******************************
//************************************************************************************

void SmallDisplacementElement::SetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    if(rVariable == CONSTITUTIVE_LAW)
    {
        if ( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize(rValues.size());

            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod  ) )
                KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
        }

        for(unsigned int i=0; i<rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i]->Clone();
        }
    }

    if(rVariable == CONSTITUTIVE_LAW_POINTER)
    {
        if ( mConstitutiveLawVector.size() != rValues.size() )
        {
            mConstitutiveLawVector.resize(rValues.size());

            if( mConstitutiveLawVector.size() != GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod  ) )
                KRATOS_ERROR( std::logic_error, "constitutive law not has the correct size ", mConstitutiveLawVector.size() );
        }

        for(unsigned int i=0; i<rValues.size(); i++)
        {
            mConstitutiveLawVector[i] = rValues[i];
        }
    }

}

//************************************************************************************
//************************************************************************************


//**********************************GET VECTOR VALUE**********************************
//************************************************************************************


void SmallDisplacementElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,
        std::vector<Vector>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

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
    else
    {

        for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }

    }

}

//***********************************GET MATRIX VALUE*********************************
//************************************************************************************

void SmallDisplacementElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    const unsigned int& integration_points_number = mConstitutiveLawVector.size();

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
    else
    {

        for ( unsigned int PointNumber = 0;  PointNumber < integration_points_number; PointNumber++ )
        {
            rValues[PointNumber] = mConstitutiveLawVector[PointNumber]->GetValue( rVariable, rValues[PointNumber] );
        }

    }


}

//********************************GET CONSTITUTIVE VALUE******************************
//************************************************************************************

void SmallDisplacementElement::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo )
{

    if(rVariable == CONSTITUTIVE_LAW || rVariable == CONSTITUTIVE_LAW_POINTER)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            rValues.resize(mConstitutiveLawVector.size());
        }

        for(unsigned int i=0; i<rValues.size(); i++)
        {
            rValues[i] = mConstitutiveLawVector[i];
        }
    }

}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void SmallDisplacementElement::Initialize()
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
    {
        mConstitutiveLawVector.resize( integration_points.size() );
    }

    //Material initialisation
    InitializeMaterial();


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::SetGeneralVariables(GeneralVariables& rVariables,
        ConstitutiveLaw::Parameters& rValues,
        const int & rPointNumber)
{

    rValues.SetStrainVector(rVariables.StrainVector);
    rValues.SetStressVector(rVariables.StressVector);
    rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rValues.SetShapeFunctionsDevivatives(rVariables.DN_DX);
    rValues.SetShapeFunctionsValues(rVariables.N);

    if(rVariables.detJ<0)
        KRATOS_ERROR(std::invalid_argument,"detJ < 0","");

    rValues.SetDeterminantF0(rVariables.detF0);
    rValues.SetDeformationGradientF0(rVariables.F0);
    rValues.SetDeterminantF(rVariables.detF);
    rValues.SetDeformationGradientF(rVariables.F);

}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::InitializeGeneralVariables (GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    unsigned int voigtsize  = 3;
    if( dimension == 3 )
    {
        voigtsize  = 6;
    }

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.B.resize( voigtsize, number_of_nodes * dimension );

    rVariables.H.resize( dimension, dimension );

    rVariables.ConstitutiveMatrix.resize( voigtsize, voigtsize );

    rVariables.StrainVector.resize( voigtsize );

    rVariables.StressVector.resize( voigtsize );

    rVariables.DN_DX.resize( number_of_nodes, dimension );

    //needed parameters for consistency with the general constitutive law: small displacements
    rVariables.detF0 = 1;
    rVariables.detF  = 1;
    rVariables.F     = identity_matrix<double>(dimension);
    rVariables.F0    = identity_matrix<double>(dimension);

    //set variables including all integration points values

    //reading shape functions
    rVariables.SetShapeFunctions(GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ));

    //reading shape functions local gradients
    rVariables.SetShapeFunctionsGradients(GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod ));

    //calculating the jacobian from cartesian coordinates to parent coordinates for all integration points
    rVariables.J = GetGeometry().Jacobian( rVariables.J, mThisIntegrationMethod );


}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        Flags& rCalculationFlags)

{

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rCalculationFlags.Is(SmallDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(SmallDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }
}



//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateElementalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        Flags& rCalculationOptions)
{
    KRATOS_TRY

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    //reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    //auxiliary terms
    Vector VolumeForce;

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values,PointNumber);

        //compute stresses and constitutive parameters
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);

        //calculating weights for integration on the "reference configuration"
        double IntegrationWeight = integration_points[PointNumber].Weight() * Variables.detJ;
        IntegrationWeight = this->CalculateIntegrationWeight( IntegrationWeight );


        //if ( dimension == 2 ) IntegrationWeight *= GetProperties()[THICKNESS];

        if ( rCalculationOptions.Is(SmallDisplacementElement::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        {
            //contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddLHS ( rLeftHandSideMatrix, Variables, IntegrationWeight );

        }

        if ( rCalculationOptions.Is(SmallDisplacementElement::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        {
            //contribution to external forces
            VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables.N );

            this->CalculateAndAddRHS ( rRightHandSideVector, Variables, VolumeForce, IntegrationWeight );

        }



    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
{

    //contributions to stiffness matrix calculated on the reference config

    // operation performed: add Km to the rLefsHandSideMatrix
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );


    //KRATOS_WATCH(rLeftHandSideMatrix)
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
{

    //contribution to external forces

    // operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables,  rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight );

    //KRATOS_WATCH(rRightHandSideVector)
}

//************************************************************************************
//************************************************************************************

double& SmallDisplacementElement::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(SmallDisplacementElement::COMPUTE_RHS_VECTOR);

    MatrixType LeftHandSideMatrix = Matrix();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( LeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( LeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );
}


//************************************************************************************
//************************************************************************************


void SmallDisplacementElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
{

    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(SmallDisplacementElement::COMPUTE_LHS_MATRIX);

    VectorType RightHandSideVector = Vector();

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, RightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( rLeftHandSideMatrix, RightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    Flags CalculationFlags;
    CalculationFlags.Set(SmallDisplacementElement::COMPUTE_LHS_MATRIX);
    CalculationFlags.Set(SmallDisplacementElement::COMPUTE_RHS_VECTOR);

    //Initialize sizes for the system components:
    this->InitializeSystemMatrices( rLeftHandSideMatrix, rRightHandSideVector, CalculationFlags );

    //Calculate elemental system
    CalculateElementalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculationFlags );

}


////************************************************************************************
////************************************************************************************

void SmallDisplacementElement::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    ClearNodalForces();

    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        mConstitutiveLawVector[i]->InitializeSolutionStep( GetProperties(),
                GetGeometry(),
                row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ),
                rCurrentProcessInfo );
}


////************************************************************************************
////************************************************************************************
void SmallDisplacementElement::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    ClearNodalForces();
}

////************************************************************************************
////************************************************************************************

void SmallDisplacementElement::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{

}

////************************************************************************************
////************************************************************************************

void SmallDisplacementElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{

    //create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    //create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();

    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {

        //compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables,PointNumber);

        //set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables,Values,PointNumber);

        //call the constitutive law to update material variables
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy(Values);

        //call the constitutive law to finalize the solution step
        mConstitutiveLawVector[PointNumber]->FinalizeSolutionStep( GetProperties(),
                GetGeometry(),
                Variables.N,
                rCurrentProcessInfo );
    }
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::InitializeMaterial()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_ERROR( std::logic_error, "a constitutive law needs to be specified for the element with ID ", this->Id() )
        KRATOS_CATCH( "" )
    }


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        double& rIntegrationWeight)

{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double Fext=0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);

        Fext = 0;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            Fext = rIntegrationWeight * rVariables.N[i] * rVolumeForce[j];
            rRightHandSideVector[index + j] += Fext;
            ExternalForce[j] +=Fext;
        }
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        double& rIntegrationWeight)
{
    KRATOS_TRY

    VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), rVariables.StressVector );
    noalias( rRightHandSideVector ) -= InternalForces;

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int indexu  = dimension * i;
        array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_INTERNAL);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            InternalForce[j] -= InternalForces [indexu+j];
        }
    }

    // std::cout<<std::endl;
    // std::cout<<" Fint "<<InternalForces<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        double& rIntegrationWeight
                                                  )
{
    KRATOS_TRY

    //contributions to stiffness matrix calculated on the reference config
    noalias( rLeftHandSideMatrix ) += prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); //to be optimized to remove the temporary

    // std::cout<<std::endl;
    // std::cout<<" Kmat "<<rLeftHandSideMatrix<<std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::ClearNodalForces()
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {

        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);
        array_1d<double, 3 > & InternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_INTERNAL);
        array_1d<double, 3 > & DynamicForce  = GetGeometry()[i].FastGetSolutionStepValue(FORCE_DYNAMIC);

        ExternalForce.clear();
        InternalForce.clear();
        DynamicForce.clear();

    }

    KRATOS_CATCH( "" )
}

//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void SmallDisplacementElement::CalculateKinematics(GeneralVariables& rVariables,
        const double& rPointNumber)

{
    KRATOS_TRY

    const GeometryType::ShapeFunctionsGradientsType& DN_De = rVariables.GetShapeFunctionsGradients();
    const Matrix& Ncontainer = rVariables.GetShapeFunctions();

    //Calculating the inverse of the jacobian and the parameters needed
    Matrix InvJ;
    MathUtils<double>::InvertMatrix( rVariables.J[rPointNumber], InvJ, rVariables.detJ);

   //Compute cartesian derivatives
    noalias( rVariables.DN_DX ) = prod( DN_De[rPointNumber] , InvJ );

    //Displacement Gradient H
    this->CalculateDisplacementGradient( rVariables.H, rVariables.DN_DX );

    //Set Shape Functions Values for this integration point
    rVariables.N=row( Ncontainer, rPointNumber);

    //Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B,rVariables.DN_DX);

 
    //Compute infinitessimal strain
    this->CalculateInfinitesimalStrain(rVariables.H,rVariables.StrainVector);


    KRATOS_CATCH( "" )
}


//*************************COMPUTE DISPLACEMENT GRADIENT******************************
//************************************************************************************

void SmallDisplacementElement::CalculateDisplacementGradient(Matrix& rH,
        const Matrix& rDN_DX)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rH = zero_matrix<double> ( dimension );

    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
        }
    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {

            array_1d<double, 3 > & Displacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            rH ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            rH ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            rH ( 0 , 2 ) += Displacement[0]*rDN_DX ( i , 2 );
            rH ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            rH ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            rH ( 1 , 2 ) += Displacement[1]*rDN_DX ( i , 2 );
            rH ( 2 , 0 ) += Displacement[2]*rDN_DX ( i , 0 );
            rH ( 2 , 1 ) += Displacement[2]*rDN_DX ( i , 1 );
            rH ( 2 , 2 ) += Displacement[2]*rDN_DX ( i , 2 );
        }
    }
    else
    {

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateInfinitesimalStrain(const Matrix& rH,
        Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
    {

        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );

        rStrainVector[0] = rH( 0, 0 );

        rStrainVector[1] = rH( 1, 1 );

        rStrainVector[2] = (rH( 0, 1 ) + rH( 1, 0 )); // xy

    }
    else if( dimension == 3 )
    {

        //Infinitesimal Strain Calculation
        if ( rStrainVector.size() != 6 ) rStrainVector.resize( 6, false );

        rStrainVector[0] = rH( 0, 0 );

        rStrainVector[1] = rH( 1, 1 );

        rStrainVector[2] = rH( 2, 2 );

        rStrainVector[3] = ( rH( 0, 1 ) + rH( 1, 0 ) ); // xy

        rStrainVector[4] = ( rH( 1, 2 ) + rH( 2, 1 ) ); // yz

        rStrainVector[5] = ( rH( 0, 2 ) + rH( 2, 0 ) ); // xz

    }
    else
    {

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )

}



//****************************COMPUTE VELOCITY GRADIENT*******************************
//************************************************************************************

void SmallDisplacementElement::CalculateVelocityGradient(const Matrix& rDN_DX,
        Matrix& rDF )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    rDF=zero_matrix<double> ( dimension );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        //Displacement from the reference to the current configuration
        array_1d<double, 3 > & CurrentVelocity  = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rDF ( j , k ) += CurrentVelocity[j]*rDN_DX ( i , k );
            }

        }

    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************
void SmallDisplacementElement::CalculateDeformationMatrix(Matrix& rB,
        const Matrix& rDN_DX)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rB.clear();
    
    if( dimension == 2 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            unsigned int index = 2 * i;

            rB( 0, index + 0 ) = rDN_DX( i, 0 );
            rB( 0, index + 1 ) = 0.0;
            rB( 1, index + 0 ) = 0.0;
            rB( 1, index + 1 ) = rDN_DX( i, 1 );
            rB( 2, index + 0 ) = rDN_DX( i, 1 );
            rB( 2, index + 1 ) = rDN_DX( i, 0 );

        }

    }
    else if( dimension == 3 )
    {

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
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

        KRATOS_ERROR( std::invalid_argument, "something is wrong with the dimension", "" );

    }

    KRATOS_CATCH( "" )
}


//************************************CALCULATE TOTAL MASS****************************
//************************************************************************************

double& SmallDisplacementElement::CalculateTotalMass( double& rTotalMass )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rTotalMass = GetGeometry().DomainSize() * GetProperties()[DENSITY];

    if( dimension == 2 )
        rTotalMass *= GetProperties()[THICKNESS];

    return rTotalMass;

    KRATOS_CATCH( "" )
}


//************************************CALCULATE VOLUME ACCELERATION*******************
//************************************************************************************

Vector& SmallDisplacementElement::CalculateVolumeForce( Vector& rVolumeForce, const Vector &rN)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    rVolumeForce = ZeroVector(dimension);
    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        if( GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION) ) //temporary, will be checked once at the beginning only
            rVolumeForce += rN[j] * GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    rVolumeForce *= GetProperties()[DENSITY];

    return rVolumeForce;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int MatSize = dimension * number_of_nodes;

    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );

    rMassMatrix = ZeroMatrix( MatSize, MatSize );

    double TotalMass = 0;
    TotalMass = this->CalculateTotalMass(TotalMass);

    Vector LumpFact = ZeroVector(number_of_nodes);

    LumpFact = GetGeometry().LumpingFactors( LumpFact );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        double temp = LumpFact[i] * TotalMass;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            unsigned int index = i * dimension + j;
            rMassMatrix( index, index ) = temp;
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( rDampMatrix.size1() != MatSize )
        rDampMatrix.resize( MatSize, MatSize, false );

    noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );


    if ( rVariable == VON_MISES_STRESS )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy (Values);

            ComparisonUtils EquivalentStress;
            rOutput[PointNumber] =  EquivalentStress.CalculateVonMises(Variables.StressVector);
        }
    }
    else
    {

        for ( unsigned int ii = 0; ii < integration_points_number; ii++ )
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rOutput[ii] );
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );



    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);


        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

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
    else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;

        }

    }
    else
    {
        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
        {
            rOutput[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable , rOutput[ii] );
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void SmallDisplacementElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

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
            if ( rOutput[PointNumber].size2() != StressVector[PointNumber].size() )
                rOutput[PointNumber].resize( 1, StressVector[PointNumber].size(), false );

            for ( unsigned int ii = 0; ii < StressVector[PointNumber].size(); ii++ )
            {
                rOutput[PointNumber]( 0, ii ) = StressVector[PointNumber][ii];
            }
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR)
    {

        std::vector<Vector> StrainVector;
        if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );
        else
            CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {

            if ( rOutput[PointNumber].size2() != StrainVector[PointNumber].size() )
                rOutput[PointNumber].resize( 1, StrainVector[PointNumber].size(), false );

            for ( unsigned int ii = 0; ii < StrainVector[PointNumber].size(); ii++ )
                rOutput[PointNumber]( 0, ii ) = StrainVector[PointNumber][ii];
        }
    }
    else if ( rVariable == CONSTITUTIVE_MATRIX )
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

        //set constitutive law flags:
        Flags &ConstitutiveLawOptions=Values.GetOptions();

        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::LAST_KNOWN_CONFIGURATION);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
            this->CalculateKinematics(Variables,PointNumber);

            //set general variables to constitutivelaw parameters
            this->SetGeneralVariables(Variables,Values,PointNumber);

            //call the constitutive law to update material variables
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponsePK2(Values);


            if( rOutput[PointNumber].size2() != Variables.ConstitutiveMatrix.size2() )
                rOutput[PointNumber].resize( Variables.ConstitutiveMatrix.size1() , Variables.ConstitutiveMatrix.size2() , false );

            rOutput[PointNumber] = Variables.ConstitutiveMatrix;

        }


    }
    else if ( rVariable == DEFORMATION_GRADIENT )  // VARIABLE SET FOR TRANSFER PURPOUSES
    {
        //create and initialize element variables:
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //reading integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics B, F, DN_DX ...
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





    KRATOS_CATCH( "" )
}



//*************************DECIMAL CORRECTION OF STRAINS******************************
//************************************************************************************

void SmallDisplacementElement::DecimalCorrection(Vector& rVector)
{
    KRATOS_TRY

    for ( unsigned int i = 0; i < rVector.size(); i++ )
    {
        if( rVector[i]*rVector[i]<1e-24 )
        {
            rVector[i]=0;
        }

    }

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
int  SmallDisplacementElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();



    //verify that the variables are correctly initialized

    if ( VELOCITY.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" );

    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" );

    if ( ACCELERATION.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" );

    if ( DENSITY.Key() == 0 )
        KRATOS_ERROR( std::invalid_argument, "DENSITY has Key zero! (check if the application is correctly registered", "" );

    // if ( BODY_FORCE.Key() == 0 )
    //     KRATOS_ERROR( std::invalid_argument, "BODY_FORCE has Key zero! (check if the application is correctly registered", "" );

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() );

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() );
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
    {
        KRATOS_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() );
    }

    //Verify that the body force is defined
    // if ( this->GetProperties().Has( BODY_FORCE ) == false )
    // {
    //     KRATOS_ERROR( std::logic_error, "BODY_FORCE not provided for property ", this->GetProperties().Id() )
    // }

    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {
        // if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 3 )
        //   KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );


        if ( THICKNESS.Key() == 0 )
            KRATOS_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered", "" );

        if ( this->GetProperties().Has( THICKNESS ) == false )
            KRATOS_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );


    }
    else
    {
        if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() != 6 )
            KRATOS_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
    }

    //check constitutive law
    for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
    {
        return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
    }

    //check if it is in the XY plane for 2D case


    return 0;

    KRATOS_CATCH( "" );
}


void SmallDisplacementElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    int IntMethod = int(mThisIntegrationMethod);
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
}

void SmallDisplacementElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
}


} // Namespace Kratos


