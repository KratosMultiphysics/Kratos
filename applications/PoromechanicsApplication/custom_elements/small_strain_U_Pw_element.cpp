//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

/* Project includes */
#include "includes/define.h"
#include "custom_elements/small_strain_U_Pw_element.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
SmallStrainUPwElement::SmallStrainUPwElement() : Element() {}

//----------------------------------------------------------------------------------------

//Constructor 1
SmallStrainUPwElement::SmallStrainUPwElement( IndexType NewId, GeometryType::Pointer pGeometry ) : Element( NewId, pGeometry ) {}

//----------------------------------------------------------------------------------------

//Constructor 2
SmallStrainUPwElement::SmallStrainUPwElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties ) : Element( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}

//----------------------------------------------------------------------------------------

//Destructor
SmallStrainUPwElement::~SmallStrainUPwElement() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Element::Pointer SmallStrainUPwElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new SmallStrainUPwElement( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//----------------------------------------------------------------------------------------

int  SmallStrainUPwElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    //verify that the variables are correctly initialized

    //Solid variables
    if ( DISPLACEMENT.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DISPLACEMENT has Key zero! (check if the application is correctly registered", "" )

    if ( VELOCITY.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "VELOCITY has Key zero! (check if the application is correctly registered", "" )

    if ( ACCELERATION.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "ACCELERATION has Key zero! (check if the application is correctly registered", "" )

    if ( DENSITY_SOLID.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_SOLID has Key zero! (check if the application is correctly registered", "" )

    //Fluid variables
    if ( WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )

    if ( DERIVATIVE_WATER_PRESSURE.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DERIVATIVE_WATER_PRESSURE has Key zero! (check if the application is correctly registered", "" )

    if ( DENSITY_WATER.Key() == 0 )
        KRATOS_THROW_ERROR( std::invalid_argument, "DENSITY_WATER has Key zero! (check if the application is correctly registered", "" )

    //verify that the dofs exist
    for ( unsigned int i = 0; i < this->GetGeometry().size(); i++ )
    {
        if ( this->GetGeometry()[i].SolutionStepsDataHas( DISPLACEMENT ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable DISPLACEMENT on node ", this->GetGeometry()[i].Id() )

        if ( this->GetGeometry()[i].HasDofFor( DISPLACEMENT_X ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Y ) == false || this->GetGeometry()[i].HasDofFor( DISPLACEMENT_Z ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing one of the dofs for the variable DISPLACEMENT on node ", GetGeometry()[i].Id() )


        if ( this->GetGeometry()[i].SolutionStepsDataHas( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing variable WATER_PRESSURE on node ", this->GetGeometry()[i].Id() )

        if ( this->GetGeometry()[i].HasDofFor( WATER_PRESSURE ) == false )
            KRATOS_THROW_ERROR( std::invalid_argument, "missing the dof for the variable WATER_PRESSURE on node ", GetGeometry()[i].Id() )
    }

    //verify that the constitutive law exists
    if ( this->GetProperties().Has( CONSTITUTIVE_LAW_POINTER ) == false )
        KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property ", this->GetProperties().Id() )

	//verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;
    this->GetProperties().GetValue( CONSTITUTIVE_LAW_POINTER )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
	    if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Infinitesimal)
		    correct_strain_measure = true;
    }

    if( correct_strain_measure == false )
	    KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the element type ", " StrainMeasure_Infinitesimal " );

    //verify that the constitutive law has the correct dimension
    if ( dimension == 2 )
    {
        if ( this->GetProperties().Has( THICKNESS ) == false )
            KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() )

        if ( THICKNESS.Key() == 0 )
            KRATOS_THROW_ERROR( std::invalid_argument, "THICKNESS has Key zero! (check if the application is correctly registered)", "" )
    }

	this->GetProperties().GetValue( CONSTITUTIVE_LAW_POINTER )->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );

    return 0;

    KRATOS_CATCH( "" );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::Initialize()
{
    KRATOS_TRY

	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    if ( GetProperties()[CONSTITUTIVE_LAW_POINTER] != NULL )
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW_POINTER]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }
    }
    else
        KRATOS_THROW_ERROR( std::logic_error, "A constitutive law needs to be specified for the element with ID ", this->Id() )

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    rElementalDofList.resize( 0 );
    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ));
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dimension + 1);

    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != MatSize )
        rLeftHandSideMatrix.resize( MatSize, MatSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize );

    //Resetting the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( MatSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int MatSize = number_of_nodes * (dimension + 1);

    //Resetting mass matrix
    if ( rMassMatrix.size1() != MatSize )
        rMassMatrix.resize( MatSize, MatSize, false );
    noalias( rMassMatrix ) = ZeroMatrix( MatSize, MatSize );

    //Defining shape functions and the determinant of the jacobian at all integration points
    Matrix Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    Vector detJcontainer;
    GetGeometry().DeterminantOfJacobian(detJcontainer,mThisIntegrationMethod);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
    double IntegrationCoefficient;
    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_WATER] + (1-Porosity)*GetProperties()[DENSITY_SOLID];
    unsigned int node;
    Matrix Nut = ZeroMatrix( dimension + 1 , number_of_nodes * (dimension + 1) );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //Setting the shape function matrix
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = (dimension+1)*i;
            Nut(0,node)=Ncontainer(PointNumber,i);
            Nut(1,node+1)=Ncontainer(PointNumber,i);
            if(dimension==3)
                Nut(2,node+2)=Ncontainer(PointNumber,i);
        }

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( IntegrationCoefficient, detJcontainer[PointNumber], integration_points[PointNumber].Weight() );

        //Adding contribution to Mass matrix
        noalias(rMassMatrix) += Density*prod(trans(Nut),Nut)*IntegrationCoefficient;
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * (dimension + 1);
    int index;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = i * (dimension + 1);
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
        }
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * (dimension + 1);
    unsigned int index;

    if ( rValues.size() != element_size )
        rValues.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = i * (dimension + 1);
        rValues[index]     = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );
        if ( dimension == 3 )
        {
            rValues[index + 2] = GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION_Z, Step );
            rValues[index + 3] = 0.0;
        }
        else
        {
            rValues[index + 2] = 0.0;
        }

    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dimension + 1);

    //Resetting the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( MatSize );

    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    //Definition of variables
    ElementalVariables Variables;
    this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    //Loop over integration points
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
    {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,PointNumber);

        //set gauss points variables to constitutivelaw parameters
        this->SetElementalVariables(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        mConstitutiveLawVector[PointNumber]->FinalizeMaterialResponseCauchy(ConstitutiveParameters);
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::SetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable,std::vector<Vector>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        mConstitutiveLawVector[PointNumber]->SetValue( rVariable, rValues[PointNumber], rCurrentProcessInfo );
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetValueOnIntegrationPoints( const Variable<double>& rVariable,std::vector<double>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == VON_MISES_STRESS )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0; i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable,std::vector<Vector>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else if ( rVariable == FLUID_FLUX )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = mConstitutiveLawVector.size();

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0;  i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,std::vector<Matrix>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    //Printing stress or strain tensors
    if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );

    else
    {
        const unsigned int& integration_points_number = mConstitutiveLawVector.size();

        if ( rValues.size() != integration_points_number )
            rValues.resize( integration_points_number );

        for ( unsigned int i = 0;  i < integration_points_number; i++ )
            rValues[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rValues[i] );
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::GetValueOnIntegrationPoints( const Variable<ConstitutiveLaw::Pointer>& rVariable,std::vector<ConstitutiveLaw::Pointer>& rValues,const ProcessInfo& rCurrentProcessInfo )
{
    if(rVariable == CONSTITUTIVE_LAW_POINTER)
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize(mConstitutiveLawVector.size());

        for(unsigned int i=0; i<rValues.size(); i++)
            rValues[i] = mConstitutiveLawVector[i];
    }
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number, false );

    if ( rVariable == VON_MISES_STRESS )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //set gauss points variables to constitutivelaw parameters
            this->SetElementalVariables(Variables,ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            PoromechanicsMathUtilities EquivalentStress;
            rOutput[PointNumber] =  EquivalentStress.CalculateVonMises(Variables.StressVector);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < integration_points_number; i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable, rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_VECTOR )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables, rCurrentProcessInfo);

        //Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //set gauss points variables to constitutivelaw parameters
            this->SetElementalVariables(Variables,ConstitutiveParameters);

            //compute constitutive tensor and/or stresses
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

            if ( rOutput[PointNumber].size() != Variables.StressVector.size() )
                rOutput[PointNumber].resize( Variables.StressVector.size(), false );

            rOutput[PointNumber] = Variables.StressVector;
        }
    }
    else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            if ( rOutput[PointNumber].size() != Variables.StrainVector.size() )
                rOutput[PointNumber].resize( Variables.StrainVector.size(), false );

            rOutput[PointNumber] = Variables.StrainVector;
        }
    }
    else if ( rVariable == FLUID_FLUX )
    {
        //Definition of variables
        ElementalVariables Variables;
        this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

        //Loop over integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            //compute element kinematics (Np, gradNpT, |J|, B, strains)
            this->CalculateKinematics(Variables,PointNumber);

            //Compute FluidFlux vector q [m/s]
            const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
            Vector BodyAcceleration = ZeroVector(dimension);
            const unsigned int number_of_nodes = GetGeometry().size();
            unsigned int Index = 0;
            for(unsigned int i = 0; i < number_of_nodes; i++)
            {
                BodyAcceleration[0] += Variables.Np[i]*Variables.BodyAcceleration[Index++];
                BodyAcceleration[1] += Variables.Np[i]*Variables.BodyAcceleration[Index++];
                if (dimension>2)
                    BodyAcceleration[2] += Variables.Np[i]*Variables.BodyAcceleration[Index++];
            }

            Vector GradPressureTerm = prod(trans(Variables.GradNpT),Variables.PressureVector);
            GradPressureTerm -= GetProperties()[DENSITY_WATER]*BodyAcceleration;

            Vector AuxFluidFlux = ZeroVector(dimension);
            AuxFluidFlux = - 1/Variables.DynamicViscosity * prod(Variables.IntrinsicPermeability, GradPressureTerm );

            Vector FluidFlux = ZeroVector(3);
            FluidFlux[0] = AuxFluidFlux[0];
            FluidFlux[1] = AuxFluidFlux[1];
            if(dimension>2)
                FluidFlux[2] = AuxFluidFlux[2];

            if ( rOutput[PointNumber].size() != 3 )
                rOutput[PointNumber].resize( 3, false );

            rOutput[PointNumber] = FluidFlux;
        }
    }
    else
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateOnIntegrationPoints( const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if ( rOutput.size() != integration_points_number )
        rOutput.resize( integration_points_number );

    if ( rVariable == CAUCHY_STRESS_TENSOR )
    {
        std::vector<Vector> StressVector;

        this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, StressVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StressVectorToTensor(StressVector[PointNumber]);
        }
    }
    else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
    {
        std::vector<Vector> StrainVector;

        CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, StrainVector, rCurrentProcessInfo );

        //loop integration points
        for ( unsigned int PointNumber = 0; PointNumber < mConstitutiveLawVector.size(); PointNumber++ )
        {
            if ( rOutput[PointNumber].size2() != dimension )
                rOutput[PointNumber].resize( dimension, dimension, false );

            rOutput[PointNumber] = MathUtils<double>::StrainVectorToTensor(StrainVector[PointNumber]);
        }
    }
    else
    {
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            rOutput[i] = mConstitutiveLawVector[i]->GetValue( rVariable , rOutput[i] );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                                bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    //Definition of variables
    ElementalVariables Variables;
    this->InitializeElementalVariables(Variables,rCurrentProcessInfo);

    //Create constitutive law parameters:
    ConstitutiveLaw::Parameters ConstitutiveParameters(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    if(CalculateLHSMatrixFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    if(CalculateResidualVectorFlag)
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

    //Loop over integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //compute element kinematics (Np, gradNpT, |J|, B, strains)
        this->CalculateKinematics(Variables,PointNumber);

        //set gauss points variables to constitutivelaw parameters
        this->SetElementalVariables(Variables,ConstitutiveParameters);

        //compute constitutive tensor and/or stresses
        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(ConstitutiveParameters);

        //calculating weighting coefficient for integration
        this->CalculateIntegrationCoefficient( Variables.IntegrationCoefficient, Variables.detJ, integration_points[PointNumber].Weight() );

        //Contributions to the left hand side
        if ( CalculateLHSMatrixFlag )
            this->CalculateAndAddLHS(rLeftHandSideMatrix, Variables);

        //Contributions to the right hand side
        if ( CalculateResidualVectorFlag )
            this->CalculateAndAddRHS(rRightHandSideVector, Variables);
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::InitializeElementalVariables (ElementalVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    //Variables at all integration points
    rVariables.NContainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );
    rVariables.DN_DlocalContainer = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
    GetGeometry().Jacobian( rVariables.JContainer, mThisIntegrationMethod );

    //Variables computed at each integration point
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;
    rVariables.B = ZeroMatrix( voigtsize, number_of_nodes * dimension );
    rVariables.StrainVector = ZeroVector( voigtsize );
    rVariables.ConstitutiveMatrix = ZeroMatrix( voigtsize, voigtsize );
    rVariables.StressVector = ZeroVector( voigtsize );

    //Needed parameters for consistency with the general constitutive law
    rVariables.detF  = 1;
    rVariables.F     = identity_matrix<double>(dimension);

    //Nodal variables
    this->InitializeNodalVariables(rVariables);

    //Properties variables
    this->InitializeProperties(rVariables);

    //ProcessInfo variables
    double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    rVariables.NewmarkCoefficient1 = rCurrentProcessInfo[GAMMA_NEWMARK]/(rCurrentProcessInfo[BETA_NEWMARK]*DeltaTime);
    rVariables.NewmarkCoefficient2 = 1/(rCurrentProcessInfo[THETA_NEWMARK]*DeltaTime);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::InitializeNodalVariables (ElementalVariables& rVariables)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    unsigned int Local_i;
    Vector BodyAccelerationAux    = ZeroVector(3);
    rVariables.BodyAcceleration   = ZeroVector(number_of_nodes * dimension);
    rVariables.DisplacementVector = ZeroVector(number_of_nodes * dimension);
    rVariables.VelocityVector     = ZeroVector(number_of_nodes * dimension);
    rVariables.PressureVector     = ZeroVector(number_of_nodes);
    rVariables.PressureDtVector   = ZeroVector(number_of_nodes);
    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        Local_i = i * dimension;
        BodyAccelerationAux = GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);

        rVariables.BodyAcceleration[Local_i]   = BodyAccelerationAux[0];
        rVariables.DisplacementVector[Local_i] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_X);
        rVariables.VelocityVector[Local_i]     = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X);
        rVariables.PressureVector[i]           = GetGeometry()[i].FastGetSolutionStepValue(WATER_PRESSURE);
        rVariables.PressureDtVector[i]         = GetGeometry()[i].FastGetSolutionStepValue(DERIVATIVE_WATER_PRESSURE);

        rVariables.BodyAcceleration[Local_i+1]   = BodyAccelerationAux[1];
        rVariables.DisplacementVector[Local_i+1] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        rVariables.VelocityVector[Local_i+1]     = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y);

        if(dimension == 3)
        {
            rVariables.BodyAcceleration[Local_i+2]   = BodyAccelerationAux[2];
            rVariables.DisplacementVector[Local_i+2] = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
            rVariables.VelocityVector[Local_i+2]     = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::InitializeProperties (ElementalVariables& rVariables)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double BulkModulus = GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));
    double BulkModulusSolid = GetProperties()[BULK_MODULUS_SOLID];
    rVariables.BiotCoefficient = 1-BulkModulus/BulkModulusSolid;
    double Porosity = GetProperties()[POROSITY];
    rVariables.BiotModulusInverse = (rVariables.BiotCoefficient-Porosity)/BulkModulusSolid + Porosity/GetProperties()[BULK_MODULUS_FLUID];
    rVariables.DynamicViscosity = GetProperties()[DYNAMIC_VISCOSITY];
    //Setting the intrinsic permeability matrix
    rVariables.IntrinsicPermeability = ZeroMatrix(dimension, dimension);
    rVariables.IntrinsicPermeability(0,0) = GetProperties()[PERMEABILITY_XX];
    rVariables.IntrinsicPermeability(1,1) = GetProperties()[PERMEABILITY_YY];
    rVariables.IntrinsicPermeability(0,1) = GetProperties()[PERMEABILITY_XY];
    rVariables.IntrinsicPermeability(1,0) = rVariables.IntrinsicPermeability(0,1);
    if(dimension==3)
    {
        rVariables.IntrinsicPermeability(2,2) = GetProperties()[PERMEABILITY_ZZ];
        rVariables.IntrinsicPermeability(2,0) = GetProperties()[PERMEABILITY_ZX];
        rVariables.IntrinsicPermeability(1,2) = GetProperties()[PERMEABILITY_YZ];
        rVariables.IntrinsicPermeability(0,2) = rVariables.IntrinsicPermeability(2,0);
        rVariables.IntrinsicPermeability(2,1) = rVariables.IntrinsicPermeability(1,2);
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateKinematics(ElementalVariables& rVariables, unsigned int PointNumber)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int node;

    //Setting the shape function vector
    rVariables.Np = row( rVariables.NContainer, PointNumber);

    //Calculating the determinant and inverse of the jacobian matrix
    Matrix InvJ;
    MathUtils<double>::InvertMatrix(rVariables.JContainer[PointNumber], InvJ, rVariables.detJ);
    if(rVariables.detJ<0)
        KRATOS_THROW_ERROR( std::invalid_argument," |J|<0 , |J| = ", rVariables.detJ )

    //Calculate shape functions global gradients
    rVariables.GradNpT = prod(rVariables.DN_DlocalContainer[PointNumber] , InvJ);

    //Compute the deformation matrix B
    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 2 * i;

            rVariables.B( 0, node + 0 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 0 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 1 ) = rVariables.GradNpT( i, 0 );
        }
    }
    else
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            node = 3 * i;

            rVariables.B( 0, node + 0 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 1, node + 1 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 2, node + 2 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 3, node + 0 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 3, node + 1 ) = rVariables.GradNpT( i, 0 );
            rVariables.B( 4, node + 1 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 4, node + 2 ) = rVariables.GradNpT( i, 1 );
            rVariables.B( 5, node + 0 ) = rVariables.GradNpT( i, 2 );
            rVariables.B( 5, node + 2 ) = rVariables.GradNpT( i, 0 );
        }
    }

    //Compute infinitessimal strain
    rVariables.StrainVector = prod(rVariables.B,rVariables.DisplacementVector);

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::SetElementalVariables(ElementalVariables& rVariables,ConstitutiveLaw::Parameters& rConstitutiveParameters)
{
    rConstitutiveParameters.SetStrainVector(rVariables.StrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
    rConstitutiveParameters.SetStressVector(rVariables.StressVector);

    //Needed parameters for consistency with the general constitutive law
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rVariables.GradNpT);
    rConstitutiveParameters.SetShapeFunctionsValues(rVariables.Np);

    rConstitutiveParameters.SetDeterminantF(rVariables.detF);
    rConstitutiveParameters.SetDeformationGradientF(rVariables.F);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateIntegrationCoefficient(double& rIntegrationCoefficient, double detJ, double weight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rIntegrationCoefficient = detJ * weight;

    if( dimension == 2 )
        rIntegrationCoefficient *= GetProperties()[THICKNESS];
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    this->CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCouplingMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix,rVariables);

    this->CalculateAndAddPermeabilityMatrix(rLeftHandSideMatrix,rVariables);
}


//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddStiffnessMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix StiffnessMatrix = prod(trans(rVariables.B), Matrix(prod(rVariables.ConstitutiveMatrix, rVariables.B)))*rVariables.IntegrationCoefficient;

    //Distribute stiffness block matrix into the elemental matrix
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i, Global_j, Local_i, Local_j;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1);
        Local_i = i * dimension;

        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1);
            Local_j = j * dimension;

            rLeftHandSideMatrix(Global_i,Global_j)     += StiffnessMatrix(Local_i,Local_j);
            rLeftHandSideMatrix(Global_i,Global_j+1)   += StiffnessMatrix(Local_i,Local_j+1);
            rLeftHandSideMatrix(Global_i+1,Global_j)   += StiffnessMatrix(Local_i+1,Local_j);
            rLeftHandSideMatrix(Global_i+1,Global_j+1) += StiffnessMatrix(Local_i+1,Local_j+1);
            if(dimension == 3)
            {
                rLeftHandSideMatrix(Global_i,Global_j+2)   += StiffnessMatrix(Local_i,Local_j+2);
                rLeftHandSideMatrix(Global_i+1,Global_j+2) += StiffnessMatrix(Local_i+1,Local_j+2);
                rLeftHandSideMatrix(Global_i+2,Global_j+1) += StiffnessMatrix(Local_i+2,Local_j+1);
                rLeftHandSideMatrix(Global_i+2,Global_j)   += StiffnessMatrix(Local_i+2,Local_j);
                rLeftHandSideMatrix(Global_i+2,Global_j+2) += StiffnessMatrix(Local_i+2,Local_j+2);
            }
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddCouplingMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;

    Vector VoigtVector = ZeroVector(voigtsize);
    VoigtVector[0] = 1;
    VoigtVector[1] = 1;
    if(dimension == 3) VoigtVector[2] = 1;

    Matrix CouplingMatrix = rVariables.BiotCoefficient*prod(trans(rVariables.B),Matrix(outer_prod(VoigtVector,rVariables.Np)))*rVariables.IntegrationCoefficient;

    //Distribute coupling block matrix into the elemental matrix
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int Global_i, Global_j, Local_i, Local_j;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1);
        Local_i = i * dimension;

        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1) + dimension;

            rLeftHandSideMatrix(Global_i,Global_j)   -= CouplingMatrix(Local_i,j);
            rLeftHandSideMatrix(Global_i+1,Global_j) -= CouplingMatrix(Local_i+1,j);
            if(dimension == 3)
                rLeftHandSideMatrix(Global_i+2,Global_j) -= CouplingMatrix(Local_i+2,j);
        }
    }

    Matrix CouplingMatrixT = rVariables.NewmarkCoefficient1*trans(CouplingMatrix);

    //Distribute transposed coupling block matrix into the elemental matrix
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1);
            Local_j = j * dimension;

            rLeftHandSideMatrix(Global_i,Global_j)   += CouplingMatrixT(i,Local_j);
            rLeftHandSideMatrix(Global_i,Global_j+1) += CouplingMatrixT(i,Local_j+1);
            if(dimension == 3)
                rLeftHandSideMatrix(Global_i,Global_j+2) += CouplingMatrixT(i,Local_j+2);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddCompressibilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix CompressibilityMatrix = rVariables.NewmarkCoefficient2*rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    //Distribute compressibility block matrix into the elemental matrix
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i, Global_j;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1) + dimension;

            rLeftHandSideMatrix(Global_i,Global_j) += CompressibilityMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddPermeabilityMatrix(MatrixType& rLeftHandSideMatrix, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1/rVariables.DynamicViscosity*
                                prod(rVariables.GradNpT,Matrix(prod(rVariables.IntrinsicPermeability,trans(rVariables.GradNpT))))*
                                rVariables.IntegrationCoefficient;

    //Distribute permeability block matrix into the elemental matrix
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i, Global_j;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
            Global_j = j * (dimension + 1) + dimension;

            rLeftHandSideMatrix(Global_i,Global_j) += PermeabilityMatrix(i,j);
        }
    }
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddRHS(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{

    this->CalculateAndAddStiffnessForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    this->CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    this->CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);
    
    this->CalculateAndAddPermeabilityFlow(rRightHandSideVector, rVariables);

    this->CalculateAndAddFluidBodyFlow(rRightHandSideVector, rVariables);
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddStiffnessForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Vector StiffnessForce = prod(trans(rVariables.B), rVariables.StressVector)*rVariables.IntegrationCoefficient;

    //Distribute stiffness block vector into elemental vector
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i, Local_i;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1);
        Local_i  = i * dimension;

        rRightHandSideVector[Global_i]   -= StiffnessForce[Local_i];
        rRightHandSideVector[Global_i+1] -= StiffnessForce[Local_i+1];
        if(dimension==3)
            rRightHandSideVector[Global_i+2] -= StiffnessForce[Local_i+2];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddMixBodyForce(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i;

    double Porosity = GetProperties()[POROSITY];
    double Density = Porosity*GetProperties()[DENSITY_WATER] + (1-Porosity)*GetProperties()[DENSITY_SOLID];

    Vector BodyAcceleration = ZeroVector(dimension);
    unsigned int Index = 0;
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        BodyAcceleration[0] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
        BodyAcceleration[1] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
        if (dimension>2)
            BodyAcceleration[2] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
    }

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        Global_i = i * (dimension + 1);

        rRightHandSideVector[Global_i]   += rVariables.Np[i] * Density * BodyAcceleration[0] * rVariables.IntegrationCoefficient;
        rRightHandSideVector[Global_i+1] += rVariables.Np[i] * Density * BodyAcceleration[1] * rVariables.IntegrationCoefficient;
        if(dimension == 3)
            rRightHandSideVector[Global_i+2] += rVariables.Np[i] * Density * BodyAcceleration[2] * rVariables.IntegrationCoefficient;
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddCouplingTerms(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int voigtsize  = 3;
    if( dimension == 3 ) voigtsize  = 6;

    Vector VoigtVector = ZeroVector(voigtsize);
    VoigtVector[0] = 1;
    VoigtVector[1] = 1;
    if(dimension == 3) VoigtVector[2] = 1;

    Matrix CouplingMatrix = rVariables.BiotCoefficient*prod(trans(rVariables.B),Matrix(outer_prod(VoigtVector,rVariables.Np)))*rVariables.IntegrationCoefficient;

    Vector CouplingForce = prod(CouplingMatrix,rVariables.PressureVector);

    //Distribute coupling block vector 1 into elemental vector
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int Global_i, Local_i;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1);
        Local_i  = i * dimension;

        rRightHandSideVector[Global_i]   += CouplingForce[Local_i];
        rRightHandSideVector[Global_i+1] += CouplingForce[Local_i+1];
        if(dimension==3)
            rRightHandSideVector[Global_i+2] += CouplingForce[Local_i+2];
    }

    Vector CouplingFlow = prod(trans(CouplingMatrix),rVariables.VelocityVector);

    //Distribute coupling block vector 2 into elemental vector
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        rRightHandSideVector[Global_i] -= CouplingFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddCompressibilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix CompressibilityMatrix = rVariables.BiotModulusInverse*outer_prod(rVariables.Np,rVariables.Np)*rVariables.IntegrationCoefficient;

    Vector CompressibilityFlow = prod(CompressibilityMatrix,rVariables.PressureDtVector);

    //Distribute compressibility block vector into elemental vector
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        rRightHandSideVector[Global_i] -= CompressibilityFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddPermeabilityFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix PermeabilityMatrix = 1/rVariables.DynamicViscosity*prod(rVariables.GradNpT,Matrix(prod(rVariables.IntrinsicPermeability,trans(rVariables.GradNpT))))*rVariables.IntegrationCoefficient;

    Vector PermeabilityFlow = prod(PermeabilityMatrix,rVariables.PressureVector);

    //Distribute permeability block vector into elemental vector
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int Global_i;

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        rRightHandSideVector[Global_i] -= PermeabilityFlow[i];
    }
}

//----------------------------------------------------------------------------------------

void SmallStrainUPwElement::CalculateAndAddFluidBodyFlow(VectorType& rRightHandSideVector, ElementalVariables& rVariables)
{
    Matrix GradNpTPerm = 1/rVariables.DynamicViscosity*GetProperties()[DENSITY_WATER]*
                         prod(rVariables.GradNpT,rVariables.IntrinsicPermeability)*rVariables.IntegrationCoefficient;

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int Global_i;

    Vector BodyAcceleration = ZeroVector(dimension);
    unsigned int Index = 0;
    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        BodyAcceleration[0] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
        BodyAcceleration[1] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
        if (dimension>2)
            BodyAcceleration[2] += rVariables.Np[i]*rVariables.BodyAcceleration[Index++];
    }

    for(unsigned int i = 0; i < number_of_nodes; i++)
    {
        Global_i = i * (dimension + 1) + dimension;

        rRightHandSideVector[Global_i] += inner_prod(row(GradNpTPerm,i),BodyAcceleration);
    }
}

} // Namespace Kratos
