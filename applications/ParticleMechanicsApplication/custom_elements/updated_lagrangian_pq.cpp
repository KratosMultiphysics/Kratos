//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_PQ.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application_variables.h"
#include "includes/checks.h"


#include "custom_utilities/mpm_explicit_utilities.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "custom_utilities/mpm_energy_calculation_utility.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianPQ::UpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry )
    : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
    mFinalizedStep = true;
}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianPQ::UpdatedLagrangianPQ(UpdatedLagrangianPQ const& rOther)
    :UpdatedLagrangian(rOther)
{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

UpdatedLagrangianPQ& UpdatedLagrangianPQ::operator=(UpdatedLagrangianPQ const& rOther)
{
    UpdatedLagrangian::operator=(rOther);
    return *this;
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianPQ::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianPQ( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer UpdatedLagrangianPQ::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< UpdatedLagrangianPQ >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianPQ::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
    UpdatedLagrangianPQ NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new UpdatedLagrangianPQ(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************
UpdatedLagrangianPQ::~UpdatedLagrangianPQ()
{
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::Initialize()
{
    KRATOS_TRY

    // Initialize parameters
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(dimension);

    // Initialize constitutive law and materials
    InitializeMaterial();

    std::cout << "===== element initialize printing geometries and elements =====" << std::endl;
    auto geom_it = pGetGeometry();

    for (size_t i = 0; i < geom_it->PointsNumber(); i++)
    {
        std::cout << "\tpoint " << i + 1 << " ID = " << geom_it->GetPoint(i).Id() << std::endl;

    }

    std::cout << "\n\nshape functions = " << geom_it->ShapeFunctionsValues() << std::endl;
    std::cout << "\n\nshape function grad 0 = " << geom_it->ShapeFunctionLocalGradient(0) << std::endl;
    std::cout << "\n\nshape function grad 1 = " << geom_it->ShapeFunctionLocalGradient(1) << std::endl;
    Matrix jac;
    geom_it->Jacobian(jac, 0);
    std::cout << "\n\nshape function jac 0 = " << jac << std::endl;
    geom_it->Jacobian(jac, 1);
    std::cout << "\n\nshape function jac 1 = " << jac << std::endl;
    std::cout << "\n\nintegration points numbers = " << geom_it->IntegrationPointsNumber() << std::endl;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
        ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
        : false;
    const SizeType def_grad_dim = (is_axisymmetric)
        ? 3
        : dimension;

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.B.resize(strain_size, number_of_nodes * dimension, false );

    rVariables.F.resize(def_grad_dim, def_grad_dim, false );

    rVariables.F0.resize(def_grad_dim, def_grad_dim, false );

    rVariables.FT.resize(def_grad_dim, def_grad_dim, false );

    rVariables.ConstitutiveMatrix.resize(strain_size, strain_size, false );

    rVariables.StrainVector.resize(strain_size, false );

    rVariables.StressVector.resize(strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );

    // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::CalculateElementalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
    const Vector& r_N = row(GetGeometry().ShapeFunctionsValues(), 0);
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    if (!is_explicit)
    {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, rCurrentProcessInfo);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, r_N);

        // Calculate Material Response
        /* NOTE:
        The function below will call CalculateMaterialResponseCauchy() by default and then (may)
        call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
        mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

        /* NOTE:
        The material points will have constant mass as defined at the beginning.
        However, the density and volume (integration weight) are changing every time step.*/
        // Update MP_Density
        mMP.density = (GetProperties()[DENSITY]) / Variables.detFT;
    }

    // The MP_Volume (integration weight) is evaluated
    mMP.volume = mMP.mass / mMP.density;

    if (CalculateStiffnessMatrixFlag && !is_explicit) // if calculation of the matrix is required
    {
        // Contributions to stiffness matrix calculated on the reference configuration
        this->CalculateAndAddLHS(
            rLeftHandSideMatrix,
            Variables,
            mMP.volume,
            rCurrentProcessInfo);
    }

    if (CalculateResidualVectorFlag) // if calculation of the vector is required
    {
        // Contribution to forces (in residual term) are calculated
        Vector volume_force = mMP.volume_acceleration * mMP.mass;
        this->CalculateAndAddRHS(
            rRightHandSideVector,
            Variables,
            volume_force,
            mMP.volume,
            rCurrentProcessInfo);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
    this->CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;
    if (is_explicit)
    {
        Matrix Jacobian;
        GetGeometry().Jacobian(Jacobian, 0);
        Matrix InvJ;
        double detJ;
        MathUtils<double>::InvertMatrix(Jacobian, InvJ, detJ);
        const Matrix& r_DN_De = GetGeometry().ShapeFunctionLocalGradient(0);
        rVariables.DN_DX = prod(r_DN_De, InvJ); // cartesian gradients

        const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
            ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
            : false;

        if (is_axisymmetric) {
            const double current_radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(
                GetGeometry().ShapeFunctionsValues(), GetGeometry());
            MPMExplicitUtilities::CalculateAndAddAxisymmetricExplicitInternalForce(*this,
                rVariables.DN_DX, mMP.cauchy_stress_vector, mMP.volume,
                mConstitutiveLawVector->GetStrainSize(), current_radius, rRightHandSideVector);
        }
        else MPMExplicitUtilities::CalculateAndAddExplicitInternalForce(*this,
            rVariables.DN_DX, mMP.cauchy_stress_vector, mMP.volume,
            mConstitutiveLawVector->GetStrainSize(), rRightHandSideVector);
    }
    else
    {
        // Operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
        this->CalculateAndAddInternalForces(rRightHandSideVector, rVariables, rIntegrationWeight);
    }
}

//************************************************************************************
//*********************Calculate the contribution of external force*******************

void UpdatedLagrangianPQ::CalculateAndAddExternalForces(
    VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dimension * i;

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rRightHandSideVector[index + j] += r_N(0, i) * rVolumeForce[j];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::CalculateExplicitStresses(const ProcessInfo& rCurrentProcessInfo,
    GeneralVariables& rVariables)
{
    KRATOS_TRY

    const bool is_axisymmetric = (rCurrentProcessInfo.Has(IS_AXISYMMETRIC))
        ? rCurrentProcessInfo.GetValue(IS_AXISYMMETRIC)
        : false;
    const bool is_pqmpm = (rCurrentProcessInfo.Has(IS_PQMPM))
        ? rCurrentProcessInfo.GetValue(IS_PQMPM)
        : false;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions = Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    // use element provided strain incremented from velocity gradient
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);

    // Compute explicit element kinematics, strain is incremented here.
    Matrix Jacobian;
    GetGeometry().Jacobian(Jacobian, 0);
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix(Jacobian, InvJ, detJ);
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Matrix r_DN_De = GetGeometry().ShapeFunctionLocalGradient(0);
    rVariables.DN_DX = prod(r_DN_De, InvJ); // cartesian gradients

    if (is_axisymmetric)
    {
        const double current_radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(r_N, GetGeometry());
        MPMExplicitUtilities::CalculateExplicitAsymmetricKinematics(rCurrentProcessInfo, *this, rVariables.DN_DX,
            mMP.almansi_strain_vector, rVariables.F, mConstitutiveLawVector->GetStrainSize(), current_radius);
    }
    else
    {
        MPMExplicitUtilities::CalculateExplicitKinematics(rCurrentProcessInfo, *this, rVariables.DN_DX,
            mMP.almansi_strain_vector, rVariables.F, mConstitutiveLawVector->GetStrainSize());
        // TODO split the above function into 1) calc vel grad, 2) calc strain increment
        // CalcVelGrad()
        // if (is_pqmpm) RecombineSubPointsIntoMasterPoint() This probably needs to be moved outside this utility so 
        // Calc strain increment and def grad of master particle
    }
    rVariables.StressVector = mMP.cauchy_stress_vector;
    rVariables.StrainVector = mMP.almansi_strain_vector;

    // Update gradient deformation
    rVariables.F0 = mDeformationGradientF0; // total member def grad NOT including this increment    
    rVariables.FT = prod(rVariables.F, rVariables.F0); // total def grad including this increment    
    rVariables.detF = MathUtils<double>::Det(rVariables.F); // det of current increment
    rVariables.detF0 = MathUtils<double>::Det(rVariables.F0); // det of def grad NOT including this increment
    rVariables.detFT = MathUtils<double>::Det(rVariables.FT); // det of total def grad including this increment
    mDeformationGradientF0 = rVariables.FT; // update member internal total grad def
    mDeterminantF0 = rVariables.detFT; // update member internal total grad def det

    // Update MP volume
    if (rCurrentProcessInfo.GetValue(IS_COMPRESSIBLE))
    {
        mMP.density = (GetProperties()[DENSITY]) / rVariables.detFT;
        mMP.volume = mMP.mass / mMP.density;
    }

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);

    // Set general variables to constitutivelaw parameters
    const Vector& r_N_vec = row(r_N, 0);
    this->SetGeneralVariables(rVariables, Values, r_N_vec);

    // Calculate Material Response
    /* NOTE:
    The function below will call CalculateMaterialResponseCauchy() by default and then (may)
    call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
    mConstitutiveLawVector->CalculateMaterialResponse(Values, rVariables.StressMeasure);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void UpdatedLagrangianPQ::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType left_hand_side_matrix = Matrix(0, 0);

    const SizeType mat_size = GetNumberOfDofs() * GetGeometry().size();
    if (rRightHandSideVector.size() != mat_size) {
        rRightHandSideVector.resize(mat_size, false);
    }
    rRightHandSideVector = ZeroVector(mat_size);

    CalculateElementalSystem(left_hand_side_matrix, rRightHandSideVector,
        rCurrentProcessInfo, false, true);
}

//*******************************************************************************************
//*******************************************************************************************
void UpdatedLagrangianPQ::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo )
{
    /* NOTE:
    In the InitializeSolutionStep of each time step the nodal initial conditions are evaluated.
    This function is called by the base scheme class.*/
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    mFinalizedStep = false;

    const bool is_explicit_central_difference = (rCurrentProcessInfo.Has(IS_EXPLICIT_CENTRAL_DIFFERENCE))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT_CENTRAL_DIFFERENCE)
        : false;

    // Calculating shape functions
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    array_1d<double,3> nodal_momentum = ZeroVector(3);
    array_1d<double,3> nodal_inertia  = ZeroVector(3);

    // Here MP contribution in terms of momentum, inertia and mass are added
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for (unsigned int j = 0; j < dimension; j++)
        {
            nodal_momentum[j] = r_N(0, i) * mMP.velocity[j] * mMP.mass;
            nodal_inertia[j] = r_N(0, i) * mMP.acceleration[j] * mMP.mass;
        }

        // Add in the predictor velocity increment for central difference explicit
        // This is the 'previous grid acceleration', which is actually
        // be the initial particle acceleration mapped to the grid.
        if (is_explicit_central_difference) {
            const double& delta_time = rCurrentProcessInfo[DELTA_TIME];
            for (unsigned int j = 0; j < dimension; j++) {
                nodal_momentum[j] += 0.5 * delta_time * (r_N(0, i) * mMP.acceleration[j]) * mMP.mass;
            }
        }

        r_geometry[i].SetLock();
        r_geometry[i].FastGetSolutionStepValue(NODAL_MOMENTUM, 0) += nodal_momentum;
        r_geometry[i].FastGetSolutionStepValue(NODAL_INERTIA, 0)  += nodal_inertia;
        r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) += r_N(0, i) * mMP.mass;
        r_geometry[i].UnSetLock();
    }
}


void UpdatedLagrangianPQ::InitializeMaterial()
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        Vector N = row(GetGeometry().ShapeFunctionsValues(), 0);
        mConstitutiveLawVector->InitializeMaterial( 
            GetProperties(), GetGeometry(), N);

        mMP.almansi_strain_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());
        mMP.cauchy_stress_vector = ZeroVector(mConstitutiveLawVector->GetStrainSize());

        // Resize the deformation gradient if we are axisymmetric
        if (mConstitutiveLawVector->GetStrainSize() == 4) mDeformationGradientF0 = IdentityMatrix(3);
    }
    else
        KRATOS_ERROR <<  "A constitutive law needs to be specified for the element with ID: " << this->Id() << std::endl;

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::ResetConstitutiveLaw()
{
    KRATOS_TRY
    GeneralVariables Variables;

    if ( GetProperties()[CONSTITUTIVE_LAW] != NULL )
    {
        mConstitutiveLawVector->ResetMaterial(
            GetProperties(),
            GetGeometry(),
            row(GetGeometry().ShapeFunctionsValues(), 0));
    }

    KRATOS_CATCH( "" )
}


//*************************COMPUTE CURRENT DISPLACEMENT*******************************
//************************************************************************************
/*
This function convert the computed nodal displacement into matrix of (number_of_nodes, dimension)
*/
Matrix& UpdatedLagrangianPQ::CalculateCurrentDisp(Matrix & rCurrentDisp, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    rCurrentDisp = ZeroMatrix(number_of_nodes, dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const array_1d<double, 3 > & current_displacement  = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rCurrentDisp(i,j) = current_displacement[j];
        }
    }

    return rCurrentDisp;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double& UpdatedLagrangianPQ::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( dimension == 2 )
        rIntegrationWeight *= GetProperties()[THICKNESS];

    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::EquationIdVector( EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo ) const
{
    const GeometryType& r_geometry = GetGeometry();
    int number_of_nodes = r_geometry.size();
    int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if ( dimension == 3 )
            rResult[index + 2] = r_geometry[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::GetDofList( DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo ) const
{
    const GeometryType& r_geometry = GetGeometry();
    rElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );

        if ( r_geometry.WorkingSpaceDimension() == 3 )
        {
            rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

}


void UpdatedLagrangianPQ::AddExplicitContribution(const VectorType& rRHSVector, const Variable<VectorType>& rRHSVariable, const Variable<array_1d<double, 3>>& rDestinationVariable, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    if (rRHSVariable == RESIDUAL_VECTOR &&
        rDestinationVariable == FORCE_RESIDUAL) {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int dimension = r_geometry.WorkingSpaceDimension();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();

        for (size_t i = 0; i < number_of_nodes; ++i) {
            size_t index = dimension * i;
            array_1d<double, 3>& r_force_residual = r_geometry[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
            for (size_t j = 0; j < dimension; ++j) {
                r_force_residual[j] += rRHSVector[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianPQ::GetValuesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            values[index + 2] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Z, Step );
    }
}


///@}
///@name Access Get Values
///@{

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == CALCULATE_EXPLICIT_MP_STRESS)
    {
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);
        this->CalculateExplicitStresses(rCurrentProcessInfo, Variables);
        this->FinalizeStepVariables(Variables, rCurrentProcessInfo);
        rValues[0] = true;
    }
    else if (rVariable == EXPLICIT_MAP_GRID_TO_MP)
    {
        MPMExplicitUtilities::UpdateGaussPointExplicit(rCurrentProcessInfo, *this);
        rValues[0] = true;
    }
    else if (rVariable == CALCULATE_MUSL_VELOCITY_FIELD)
    {
        MPMExplicitUtilities::CalculateMUSLGridVelocity(rCurrentProcessInfo, *this);
        rValues[0] = true;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_MATERIAL_ID) {
        rValues[0] = GetProperties().Id();
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_DENSITY) {
        rValues[0] = mMP.density;
    }
    else if (rVariable == MP_MASS) {
        rValues[0] = mMP.mass;
    }
    else if (rVariable == MP_VOLUME) {
        rValues[0] = mMP.volume;
    }
    else if (rVariable == MP_POTENTIAL_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculatePotentialEnergy(*this);
    }
    else if (rVariable == MP_KINETIC_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateKineticEnergy(*this);
    }
    else if (rVariable == MP_STRAIN_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateStrainEnergy(*this);
    }
    else if (rVariable == MP_TOTAL_ENERGY) {
        rValues[0] = MPMEnergyCalculationUtility::CalculateTotalEnergy(*this);
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > >& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        rValues[0] = mMP.xg;
    }
    else if (rVariable == MP_DISPLACEMENT) {
        rValues[0] = mMP.displacement;
    }
    else if (rVariable == MP_VELOCITY) {
        rValues[0] = mMP.velocity;
    }
    else if (rVariable == MP_ACCELERATION) {
        rValues[0] = mMP.acceleration;
    }
    else if (rVariable == MP_VOLUME_ACCELERATION) {
        rValues[0] = mMP.volume_acceleration;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
        rValues[0] = mMP.cauchy_stress_vector;
    }
    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
        rValues[0] = mMP.almansi_strain_vector;
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in CalculateOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

///@}
///@name Access Set Values
///@{

void UpdatedLagrangianPQ::SetValuesOnIntegrationPoints(const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
}

void UpdatedLagrangianPQ::SetValuesOnIntegrationPoints(const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_MASS) {
        mMP.mass = rValues[0];
    }
    else if (rVariable == MP_DENSITY) {
        mMP.density = rValues[0];
    }
    else if (rVariable == MP_VOLUME) {
        mMP.volume = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::SetValuesOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > > rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_COORD || rVariable == MPC_COORD) {
        mMP.xg = rValues[0];
    }
    else if (rVariable == MP_DISPLACEMENT) {
        mMP.displacement = rValues[0];
    }
    else if (rVariable == MP_VELOCITY) {
        mMP.velocity = rValues[0];
    }
    else if (rVariable == MP_ACCELERATION) {
        mMP.acceleration = rValues[0];
    }
    else if (rVariable == MP_VOLUME_ACCELERATION) {
        mMP.volume_acceleration = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

void UpdatedLagrangianPQ::SetValuesOnIntegrationPoints(const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rValues.size() > 1)
        << "Only 1 value per integration point allowed! Passed values vector size: "
        << rValues.size() << std::endl;

    if (rVariable == MP_CAUCHY_STRESS_VECTOR) {
        mMP.cauchy_stress_vector = rValues[0];
    }
    else if (rVariable == MP_ALMANSI_STRAIN_VECTOR) {
        mMP.almansi_strain_vector = rValues[0];
    }
    else
    {
        KRATOS_ERROR << "Variable " << rVariable << " is called in SetValuesOnIntegrationPoints, but is not implemented." << std::endl;
    }
}

///@}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  UpdatedLagrangianPQ::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    Element::Check(rCurrentProcessInfo);

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Verify compatibility with the constitutive law
    ConstitutiveLaw::Features LawFeatures;

    this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

    bool correct_strain_measure = false;
    for(unsigned int i=0; i<LawFeatures.mStrainMeasures.size(); i++)
    {
        if(LawFeatures.mStrainMeasures[i] == ConstitutiveLaw::StrainMeasure_Deformation_Gradient)
            correct_strain_measure = true;
    }

    KRATOS_ERROR_IF(correct_strain_measure == false ) << "Constitutive law is not compatible with the element type: Large Displacements " << std::endl;

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)

    // Verify that the dofs exist
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    if( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false)
    {
        KRATOS_ERROR << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;
    }
    else
    {
        // Verify that the constitutive law has the correct dimension
        if ( dimension == 2 )
        {
            KRATOS_CHECK_VARIABLE_KEY(THICKNESS)
            KRATOS_ERROR_IF_NOT(this->GetProperties().Has( THICKNESS )) << "THICKNESS not provided for element " << this->Id() << std::endl;
        }
        else
        {
            KRATOS_ERROR_IF_NOT(this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize() == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) " << this->Id() << std::endl;
        }

        // Check constitutive law
        this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Check( this->GetProperties(), r_geometry, rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

void UpdatedLagrangianPQ::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )

    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void UpdatedLagrangianPQ::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}


} // Namespace Kratos

