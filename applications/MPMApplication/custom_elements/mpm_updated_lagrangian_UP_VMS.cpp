//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Laura Moreno Martínez
//


// System includes
#include <cmath>
#include <limits>

// Project includes
#include "includes/define.h"
#include "custom_elements/mpm_updated_lagrangian.hpp"
#include "custom_elements/mpm_updated_lagrangian_UP_VMS.hpp"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{


    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt2D3C [3][2] = { {0, 0}, {1, 1}, {0, 1} };


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS()
    : MPMUpdatedLagrangianUP()
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry )
    : MPMUpdatedLagrangianUP( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : MPMUpdatedLagrangianUP( NewId, pGeometry, pProperties )
{
}
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS::MPMUpdatedLagrangianUPVMS( MPMUpdatedLagrangianUPVMS const& rOther)
    : MPMUpdatedLagrangianUP(rOther)

{
}

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

MPMUpdatedLagrangianUPVMS&  MPMUpdatedLagrangianUPVMS::operator=(MPMUpdatedLagrangianUPVMS const& rOther)
{
    MPMUpdatedLagrangianUP::operator=(rOther);


    return *this;
}
//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUPVMS::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new MPMUpdatedLagrangianUPVMS( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Element::Pointer MPMUpdatedLagrangianUPVMS::Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive< MPMUpdatedLagrangianUPVMS >(NewId, pGeom, pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer MPMUpdatedLagrangianUPVMS::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    MPMUpdatedLagrangianUPVMS NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    NewElement.m_mp_pressure = m_mp_pressure;

    NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();

    NewElement.mDeformationGradientF0 = mDeformationGradientF0;

    NewElement.mDeterminantF0 = mDeterminantF0;

    return Element::Pointer( new MPMUpdatedLagrangianUPVMS(NewElement) );
}
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
MPMUpdatedLagrangianUPVMS::~MPMUpdatedLagrangianUPVMS()
{
}


//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::InitializeGeneralVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MPMUpdatedLagrangian::InitializeGeneralVariables(rVariables, rCurrentProcessInfo);
    GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int voigt_dimension = (dimension - 1) * (dimension - 1) + 2;

    // Initialize stabilization parameters
    rVariables.tau1 = 0;
    rVariables.tau2 = 0;
    rVariables.BulkModulus = CalculateBulkModulus();
    rVariables.DynamicCoefficient = 0;
    rVariables.DynamicRHS = ZeroVector(dimension);
    rVariables.DiscreteAcceleration = ZeroVector(dimension);

    // Set Pressure and Pressure Gradient in gauss points
    rVariables.PressureGP = 0;
    rVariables.PressureGradient = ZeroVector(dimension);

    // Set identity tensor matrix
    rVariables.TensorIdentityMatrix = ZeroMatrix(voigt_dimension, voigt_dimension);

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************


void MPMUpdatedLagrangianUPVMS::CalculateElementalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    // Create and initialize element variables:
    GeneralVariables Variables;
    this->InitializeGeneralVariables(Variables, rCurrentProcessInfo);

    const Vector& r_N = row(GetGeometry().ShapeFunctionsValues(), 0);
    const bool is_explicit = (rCurrentProcessInfo.Has(IS_EXPLICIT))
        ? rCurrentProcessInfo.GetValue(IS_EXPLICIT)
        : false;

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& r_constitutive_law_options = Values.GetOptions();
    r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    if (!is_explicit)
    {
        r_constitutive_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        r_constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematics(Variables, rCurrentProcessInfo);

        // Set general variables to constitutivelaw parameters
        this->SetGeneralVariables(Variables, Values, r_N);

        // Calculate Material Response
        /* NOTE:
        The function below will call CalculateMaterialResponseCauchy() by default and then (may)
        call CalculateMaterialResponseKirchhoff() in the constitutive_law.*/
        mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);

        // The material points have constant mass, while density and volume
        // are updated every time step.
        // Update MP_Density
        mMP.density = GetProperties()[DENSITY] / Variables.detFT;


        // Compute ASGS stabilization variables
        CalculateStabilizationVariables(Variables);

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
        Vector volume_force = (mMP.volume_acceleration * mMP.mass) + (mMP.body_force * mMP.mass);
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

void MPMUpdatedLagrangianUPVMS::CalculateStabilizationVariables(GeneralVariables& rVariables)
{
    KRATOS_TRY

    CalculatePressureVariables(rVariables);
    ConvertPressureGradientInVoigt(rVariables.PressureGradient, rVariables.PressureGradientVoigt);
    CalculateTensorIdentityMatrix(rVariables.TensorIdentityMatrix);
    CalculateTaus(rVariables);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculatePressureVariables(GeneralVariables& rVariables)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    for ( unsigned int j = 0; j < number_of_nodes; j++ )
    {
        const double nodal_pressure = r_geometry[j].FastGetSolutionStepValue(PRESSURE);
        rVariables.PressureGP += r_N(0, j) * nodal_pressure;

        for ( unsigned int i = 0; i < dimension; i++ )
        {
            rVariables.PressureGradient[i] += rVariables.DN_DX(j, i) * nodal_pressure;
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateDynamicStabilizationVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    const double delta_time = rCurrentProcessInfo[DELTA_TIME];
    KRATOS_ERROR_IF(delta_time <= std::numeric_limits<double>::epsilon())
        << "DELTA_TIME must be positive to calculate dynamic ASGS stabilization variables." << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(NEWMARK_BETA))
        << "NEWMARK_BETA is required to calculate dynamic ASGS stabilization variables." << std::endl;
    const double beta = rCurrentProcessInfo[NEWMARK_BETA];
    KRATOS_ERROR_IF(beta <= std::numeric_limits<double>::epsilon())
        << "NEWMARK_BETA must be positive to calculate dynamic ASGS stabilization variables." << std::endl;

    array_1d<double, 3> previous_velocity = ZeroVector(3);
    array_1d<double, 3> previous_acceleration = ZeroVector(3);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3>& r_nodal_previous_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1);
        const array_1d<double, 3>& r_nodal_previous_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, 1);

        noalias(previous_velocity) += r_N(0, i) * r_nodal_previous_velocity;
        noalias(previous_acceleration) += r_N(0, i) * r_nodal_previous_acceleration;
    }

    const double density = GetProperties()[DENSITY] / rVariables.detFT;
    const double velocity_coefficient = 1.0 / (beta * delta_time);
    const double acceleration_coefficient = 0.5 / beta - 1.0;

    rVariables.DynamicCoefficient = density / (beta * delta_time * delta_time);
    rVariables.DynamicRHS = ZeroVector(dimension);
    rVariables.DiscreteAcceleration = ZeroVector(dimension);

    for (unsigned int i = 0; i < dimension; ++i) {
        rVariables.DynamicRHS[i] = density * (
            velocity_coefficient * previous_velocity[i]
            + acceleration_coefficient * previous_acceleration[i]);
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

double MPMUpdatedLagrangianUPVMS::CalculateElementSize() const
{

    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry().GetGeometryParent(0);
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();

    if (dimension == 2) {
        if (number_of_nodes == 3) {
            return ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 4) {
            return ElementSizeCalculator<2, 4>::AverageElementSize(r_geometry);
        }
    } else if (dimension == 3) {
        if (number_of_nodes == 4) {
            return ElementSizeCalculator<3, 4>::AverageElementSize(r_geometry);
        } else if (number_of_nodes == 8) {
            return ElementSizeCalculator<3, 8>::AverageElementSize(r_geometry);
        }
    }

    KRATOS_ERROR << "Unsupported geometry in MPMUpdatedLagrangianUPVMS element size calculation. Dimension: "
                 << dimension << ", number of nodes: " << number_of_nodes << std::endl;

    KRATOS_CATCH("")
}

// Calculate stabilization parameters

void MPMUpdatedLagrangianUPVMS::CalculateTaus(GeneralVariables& rVariables)
{
    KRATOS_TRY

    const double constant1 = 1.0;
    const double constant2 = 1.0;
    const double shear_modulus = CalculateShearModulus();
    const double characteristic_element_size = CalculateElementSize();

    rVariables.tau1 = constant1 * characteristic_element_size * characteristic_element_size / (2.0 * shear_modulus);
    rVariables.tau2 = 2 * constant2 * shear_modulus;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

double MPMUpdatedLagrangianUPVMS::CalculateShearModulus() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(GetProperties().Has(YOUNG_MODULUS))
        << "YOUNG_MODULUS is required to calculate ASGS stabilization parameters." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(POISSON_RATIO))
        << "POISSON_RATIO is required to calculate ASGS stabilization parameters." << std::endl;

    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio = GetProperties()[POISSON_RATIO];

    return young_modulus / (2.0 * (1.0 + poisson_ratio));

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

double MPMUpdatedLagrangianUPVMS::CalculateBulkModulus() const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(GetProperties().Has(YOUNG_MODULUS))
        << "YOUNG_MODULUS is required to calculate ASGS stabilization parameters." << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(POISSON_RATIO))
        << "POISSON_RATIO is required to calculate ASGS stabilization parameters." << std::endl;

    const double& young_modulus = GetProperties()[YOUNG_MODULUS];
    const double& poisson_ratio = GetProperties()[POISSON_RATIO];
    const double tolerance = 10.e-7;
    const double denominator = 3.0 * (1.0 - 2.0 * poisson_ratio);

    if (std::abs(denominator) < tolerance) {
        return 1.0e16;
    }

    const double bulk_modulus = young_modulus / denominator;

    return std::isfinite(bulk_modulus) ? bulk_modulus : 1.0e16;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************


void MPMUpdatedLagrangianUPVMS::CalculateTensorIdentityMatrix(Matrix& rTensorIdentityMatrix)
{

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    const unsigned int voigt_dimension = (dimension - 1) * (dimension - 1) + 2;
    rTensorIdentityMatrix = ZeroMatrix(voigt_dimension, voigt_dimension);

    if (dimension == 2)
    {
        for (unsigned int i = 0; i < voigt_dimension; i++)
        {
            for (unsigned int j = 0; j < voigt_dimension; j++)
            {
                rTensorIdentityMatrix(i, j) = CalculateTensorIdentityComponent(
                    this->msIndexVoigt2D3C[i][0], this->msIndexVoigt2D3C[i][1],
                    this->msIndexVoigt2D3C[j][0], this->msIndexVoigt2D3C[j][1]);
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < voigt_dimension; i++)
        {
            for (unsigned int j = 0; j < voigt_dimension; j++)
            {
                rTensorIdentityMatrix(i, j) = CalculateTensorIdentityComponent(
                    msIndexVoigt3D6C[i][0], msIndexVoigt3D6C[i][1],
                    msIndexVoigt3D6C[j][0], msIndexVoigt3D6C[j][1]);
            }
        }
    }

}

//************************************************************************************
//************************************************************************************

double MPMUpdatedLagrangianUPVMS::CalculateTensorIdentityComponent(
    const unsigned int a,
    const unsigned int b,
    const unsigned int c,
    const unsigned int d) const
{

    const auto delta = [](const unsigned int i, const unsigned int j) {
        return i == j ? 1.0 : 0.0;
    };

    const double identity_product = delta(a, b) * delta(c, d);
    const double symmetric_identity = 0.5 * (delta(a, c) * delta(b, d) + delta(a, d) * delta(b, c));

    return identity_product - 2.0 * symmetric_identity;
}


//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::ConvertPressureGradientInVoigt(
    const Vector& rPressureGradient,
    Vector& rPressureGradientVoigt)
{
    KRATOS_TRY
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int voigt_dimension = (dimension - 1) * (dimension - 1) + 2;
    rPressureGradientVoigt = ZeroVector(voigt_dimension);

    if (dimension == 2)
    {
        rPressureGradientVoigt(0) = rPressureGradient(0);
        rPressureGradientVoigt(1) = rPressureGradient(1);
        rPressureGradientVoigt(2) = 0.5 * (rPressureGradient(0) + rPressureGradient(1));
    }
    else if (dimension == 3)
    {
        rPressureGradientVoigt(0) = rPressureGradient(0);
        rPressureGradientVoigt(1) = rPressureGradient(1);
        rPressureGradientVoigt(2) = rPressureGradient(2);
        rPressureGradientVoigt(3) = 0.5 * (rPressureGradient(0) + rPressureGradient(1));
        rPressureGradientVoigt(4) = 0.5 * (rPressureGradient(2) + rPressureGradient(1));
        rPressureGradientVoigt(5) = 0.5 * (rPressureGradient(2) + rPressureGradient(0));
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddRHS(
    VectorType& rRightHandSideVector,
    GeneralVariables& rVariables,
    Vector& rVolumeForce,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{

    // Operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // Operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= Stabilized terms of the momentum equation
    CalculateAndAddStabilizedDisplacementVMS( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

    // Operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressureVMS( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight);

}

//************************************************************************************
//************************************************************************************




void MPMUpdatedLagrangianUPVMS::CalculateAndAddStabilizedDisplacementVMS(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    Vector displacement_test_pressure_gradient_projection = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector displacement_test_deviatoric_pressure_gradient_projection = prod(
        prod(trans(rVariables.B), rVariables.TensorIdentityMatrix),
        rVariables.PressureGradientVoigt);
    const double volumetric_strain = this->CalculateVolumetricStrainFunction(rVariables);
    const double volumetric_residual = rVariables.PressureGP / rVariables.BulkModulus - volumetric_strain;

    unsigned int indexi = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_up = dimension * i + i;
        for ( unsigned int jdim = 0; jdim < dimension; jdim ++ )
        {
            const double momentum_residual = -rVolumeForce[jdim] / rIntegrationWeight
                - rVariables.PressureGradient[jdim];

            rRightHandSideVector[index_up + jdim] -= rVariables.tau1
                * momentum_residual
                * displacement_test_pressure_gradient_projection(i)
                * rIntegrationWeight;

            rRightHandSideVector[index_up + jdim] -= rVariables.tau1
                * momentum_residual
                * displacement_test_deviatoric_pressure_gradient_projection(indexi)
                * rIntegrationWeight;

            rRightHandSideVector[index_up + jdim] += rVariables.tau2
                * volumetric_residual
                * rVariables.DN_DX(i, jdim)
                * rIntegrationWeight;
            indexi++;
        }
    }


    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddStabilizedPressureVMS(VectorType& rRightHandSideVector,
        GeneralVariables& rVariables,
        Vector& rVolumeForce,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int index_p = dimension;
    const double volumetric_strain = this->CalculateVolumetricStrainFunction(rVariables);
    const double volumetric_strain_linearization = this->CalculateVolumetricStrainLinearization(rVariables);
    const double volumetric_residual = rVariables.PressureGP / rVariables.BulkModulus - volumetric_strain;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idime = 0; idime < dimension; idime++ )
        {
            const double momentum_residual = -rVariables.PressureGradient[idime]
                - rVolumeForce[idime] / rIntegrationWeight;

            rRightHandSideVector[index_p] -= rVariables.tau1
                * volumetric_strain_linearization
                * rVariables.DN_DX(i, idime)
                * momentum_residual
                * rIntegrationWeight;
        }

        rRightHandSideVector[index_p] -= rVariables.tau2
            * volumetric_residual
            * r_N(0, i)
            * (1.0 / rVariables.BulkModulus)
            * rIntegrationWeight;

        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}
//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddLHS(
    MatrixType& rLeftHandSideMatrix,
    GeneralVariables& rVariables,
    const double& rIntegrationWeight,
    const ProcessInfo& rCurrentProcessInfo)
{

    // Operation performed: add Km to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kg to the rLefsHandSideMatrix
    if (!rCurrentProcessInfo.Has(IGNORE_GEOMETRIC_STIFFNESS))
    {
        CalculateAndAddKuugUP(rLeftHandSideMatrix, rVariables, rIntegrationWeight);
    }

    // Operation performed: add Kup to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKup( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operations performed: add Kpu to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp to the rLefsHandSideMatrix
    MPMUpdatedLagrangianUP::CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kuu stabilization to the rLefsHandSideMatrix
    CalculateAndAddKuuStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kup stabilization to the rLefsHandSideMatrix
    CalculateAndAddKupStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operations performed: add Kpu stabilization to the rLefsHandSideMatrix
    CalculateAndAddKpuStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

    // Operation performed: add Kpp stabilization to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, rIntegrationWeight );

}

//************************************************************************************
//************************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddKuuStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    Vector row_displacement_pressure_gradient_projection = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector row_displacement_deviatoric_pressure_gradient_projection = prod(
        Matrix(trans(prod(rVariables.TensorIdentityMatrix, rVariables.B))),
        rVariables.PressureGradientVoigt);
    Vector column_displacement_pressure_gradient_projection = row_displacement_pressure_gradient_projection;
    Vector column_displacement_deviatoric_pressure_gradient_projection = prod(
        Matrix(prod(trans(rVariables.B), rVariables.TensorIdentityMatrix)),
        rVariables.PressureGradientVoigt);
    const double volumetric_strain_linearization = this->CalculateVolumetricStrainLinearization(rVariables);


    unsigned int indexi = 0;
    unsigned int indexj = 0;

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int idim = 0; idim < dimension ; idim ++)
        {
            indexj = 0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
                for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
                {
                    rLeftHandSideMatrix(indexi + i, indexj + j) -= rVariables.tau1
                        * row_displacement_pressure_gradient_projection(i)
                        * column_displacement_pressure_gradient_projection(j)
                        * rIntegrationWeight;

                    rLeftHandSideMatrix(indexi + i, indexj + j) -= rVariables.tau1
                        * row_displacement_deviatoric_pressure_gradient_projection(indexi)
                        * column_displacement_pressure_gradient_projection(j)
                        * rIntegrationWeight;

                    rLeftHandSideMatrix(indexi + i, indexj + j) -= rVariables.tau1
                        * row_displacement_pressure_gradient_projection(i)
                        * column_displacement_deviatoric_pressure_gradient_projection(indexj)
                        * rIntegrationWeight;

                    rLeftHandSideMatrix(indexi + i, indexj + j) -= rVariables.tau1
                        * row_displacement_deviatoric_pressure_gradient_projection(indexi)
                        * column_displacement_deviatoric_pressure_gradient_projection(indexj)
                        * rIntegrationWeight;

                    rLeftHandSideMatrix(indexi + i, indexj + j) += rVariables.tau2
                        * volumetric_strain_linearization
                        * rVariables.DN_DX(i, idim)
                        * rVariables.DN_DX(j, jdim)
                        * rIntegrationWeight;
                    indexj++;
                }
            }
            indexi++;
        }
    }

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************
void MPMUpdatedLagrangianUPVMS::CalculateAndAddKupStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Vector row_displacement_pressure_gradient_projection = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector column_pressure_deviatoric_pressure_gradient_projection = prod(
        Matrix(trans(prod(rVariables.TensorIdentityMatrix, rVariables.B))),
        rVariables.PressureGradientVoigt);

    // Assemble components considering added DOF matrix system
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_p  = dimension;
        unsigned int index_up = dimension * i + i;

        unsigned int indexj = 0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            const double pressure_compressibility_shape_function = r_N(0, j) / rVariables.BulkModulus;

            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1
                    * row_displacement_pressure_gradient_projection(i)
                    * rVariables.DN_DX(j, k)
                    * rIntegrationWeight;

                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau1
                    * column_pressure_deviatoric_pressure_gradient_projection(indexj)
                    * rVariables.DN_DX(i, k)
                    * rIntegrationWeight;

                rLeftHandSideMatrix(index_up + k, index_p) -= rVariables.tau2
                    * pressure_compressibility_shape_function
                    * rVariables.DN_DX(i, k)
                    * rIntegrationWeight;
                indexj++;
            }
            index_p += (dimension + 1);
        }
    }


    KRATOS_CATCH( "" )
}



//***********************************************************************************
//***********************************************************************************
void MPMUpdatedLagrangianUPVMS::CalculateAndAddKpuStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    Vector column_displacement_pressure_gradient_projection = prod(rVariables.DN_DX, rVariables.PressureGradient);
    Vector column_displacement_deviatoric_pressure_gradient_projection = prod(
        Matrix(prod(trans(rVariables.B), rVariables.TensorIdentityMatrix)),
        rVariables.PressureGradientVoigt);
    const double volumetric_strain_linearization = this->CalculateVolumetricStrainLinearization(rVariables);

    // Assemble components considering added DOF matrix system
    unsigned int index_p = dimension;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const double pressure_compressibility_shape_function = r_N(0, i) / rVariables.BulkModulus;
        unsigned int indexj = 0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            unsigned int index_up = dimension * j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau1
                    * volumetric_strain_linearization
                    * rVariables.DN_DX(i, k)
                    * column_displacement_pressure_gradient_projection(j)
                    * rIntegrationWeight;

                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau1
                    * volumetric_strain_linearization
                    * rVariables.DN_DX(i, k)
                    * column_displacement_deviatoric_pressure_gradient_projection(indexj)
                    * rIntegrationWeight;

                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau2
                    * volumetric_strain_linearization
                    * pressure_compressibility_shape_function
                    * rVariables.DN_DX(j, k)
                    * rIntegrationWeight;

                indexj++;
            }
        }
        index_p += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();
    unsigned int indexpi = dimension;
    Matrix pressure_gradient_matrix = prod(rVariables.DN_DX, trans(rVariables.DN_DX));
    const double volumetric_strain_linearization = this->CalculateVolumetricStrainLinearization(rVariables);
    const double inverse_bulk_modulus = 1.0 / rVariables.BulkModulus;

    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        unsigned int indexpj = dimension;
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            rLeftHandSideMatrix(indexpi, indexpj) -= rVariables.tau1
                * volumetric_strain_linearization
                * pressure_gradient_matrix(i, j)
                * rIntegrationWeight;

            rLeftHandSideMatrix(indexpi, indexpj) += rVariables.tau2
                * inverse_bulk_modulus
                * inverse_bulk_modulus
                * r_N(0, i)
                * r_N(0, j)
                * rIntegrationWeight;

            indexpj += (dimension + 1);
        }

        indexpi += (dimension + 1);
    }

    KRATOS_CATCH( "" )
}



void MPMUpdatedLagrangianUPVMS::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MPMUpdatedLagrangianUP )
}

void MPMUpdatedLagrangianUPVMS::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MPMUpdatedLagrangianUP )
}
} // Namespace Kratos
