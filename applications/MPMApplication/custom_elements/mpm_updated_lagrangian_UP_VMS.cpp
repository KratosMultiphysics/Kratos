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
#include <omp.h>
#include <sstream>
#include <cmath>
#include <limits>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/mpm_updated_lagrangian.hpp"
#include "custom_elements/mpm_updated_lagrangian_UP_VMS.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"
#include "includes/checks.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_utilities/material_point_generator_utility.h"
#include "utilities/element_size_calculator.h"

namespace Kratos
{


    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt3D6C [6][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1}, {1, 2}, {0, 2} };
    const unsigned int MPMUpdatedLagrangianUPVMS::msIndexVoigt2D4C [4][2] = { {0, 0}, {1, 1}, {2, 2}, {0, 1} };
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

    // Set Identity matrices
    rVariables.TensorIdentityMatrix = ZeroMatrix(voigt_dimension,voigt_dimension);

    // Set variables for OSGS stabilization
    rVariables.ResProjDisplGP = ZeroVector(dimension);
    rVariables.ResProjPressGP = 0;

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
        CalculateStabilizationVariables(Variables, rCurrentProcessInfo);

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
        mMP.body_force = this->ComputeMaterialPointBodyForce();

        Vector volume_force = (mMP.volume_acceleration * mMP.mass ) + (mMP.body_force * mMP.mass);
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

void MPMUpdatedLagrangianUPVMS::CalculateStabilizationVariables(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    CalculatePressureVariables(rVariables);
    ConvertPressureGradientInVoigt(rVariables.PressureGradient, rVariables.PressureGradientVoigt);
    CalculateTensorIdentityMatrix(rVariables.TensorIdentityMatrix);
    CalculateTaus(rVariables, rCurrentProcessInfo);

    const bool is_dynamic = rCurrentProcessInfo.Has(IS_DYNAMIC)
        ? rCurrentProcessInfo.GetValue(IS_DYNAMIC)
        : false;
    if (is_dynamic && !GetProperties().Has(DYNAMIC_VISCOSITY)) { //
        CalculateDynamicStabilizationVariables(rVariables, rCurrentProcessInfo);
    }

    // Compute Residual Projection in integration points
    if (rCurrentProcessInfo.GetValue(STABILIZATION_TYPE) == 3)
    {
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        const unsigned int dimension = r_geometry.WorkingSpaceDimension();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        array_1d<double, 3 > nodal_resprojdispl = ZeroVector(3);
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            rVariables.ResProjPressGP += r_N(0,j) * r_geometry[j].FastGetSolutionStepValue(RESPROJ_PRESS);

            nodal_resprojdispl = r_geometry[j].FastGetSolutionStepValue(RESPROJ_DISPL,0);
            for (unsigned int k = 0; k < dimension; k++)
            {
               rVariables.ResProjDisplGP[k] += r_N(0, j) * nodal_resprojdispl[k];
            }
        }

    }
    // Compute Shear modulus and Bulk Modulus
    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO))
    {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];
        rVariables.ShearModulus = young_modulus / (2.0 * (1.0 + poisson_ratio));
        rVariables.BulkModulus  = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));

        const double tolerance = 10.e-7;
        const bool check = bool( (poisson_ratio > 0.5-tolerance ) || (poisson_ratio < (-1.0 + tolerance)) );
        if (POISSON_RATIO.Key() == 0 || check==true)  rVariables.BulkModulus= 1e16;

        if (rVariables.BulkModulus!=rVariables.BulkModulus)
            rVariables.BulkModulus= 1e16;
    }
    else if (GetProperties().Has(DYNAMIC_VISCOSITY))
    {
        rVariables.ShearModulus = GetProperties()[DYNAMIC_VISCOSITY];
        rVariables.BulkModulus = GetProperties()[BULK_MODULUS];
    }

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
    // rho * a_{n+1} = DynamicCoefficient * delta_displacement_current - known_dynamic_terms.
    // DynamicRHS stores the full dynamic residual contribution to be subtracted in
    // the stabilized momentum residual.

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
    array_1d<double, 3> current_displacement = ZeroVector(3);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3>& r_nodal_previous_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY, 1);
        const array_1d<double, 3>& r_nodal_previous_acceleration = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, 1);

        noalias(previous_velocity) += r_N(0, i) * r_nodal_previous_velocity;
        noalias(previous_acceleration) += r_N(0, i) * r_nodal_previous_acceleration;
        for (unsigned int j = 0; j < dimension; ++j) {
            current_displacement[j] += r_N(0, i) * rVariables.CurrentDisp(i, j);
        }
    }

    const double density = mMP.density;
    const double velocity_coefficient = 1.0 / (beta * delta_time);
    const double acceleration_coefficient = 0.5 / beta - 1.0;

    rVariables.DynamicCoefficient = density / (beta * delta_time * delta_time);
    rVariables.DynamicRHS = ZeroVector(dimension);
    rVariables.DiscreteAcceleration = ZeroVector(dimension);

    for (unsigned int i = 0; i < dimension; ++i) {
        const double known_dynamic_terms = density * (
            velocity_coefficient * previous_velocity[i]
            + acceleration_coefficient * previous_acceleration[i]);
        rVariables.DynamicRHS[i] = known_dynamic_terms
            - rVariables.DynamicCoefficient * current_displacement[i];
        rVariables.DiscreteAcceleration[i] = rVariables.DynamicCoefficient / density * current_displacement[i]
            - known_dynamic_terms / density;
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

void MPMUpdatedLagrangianUPVMS::CalculateTaus(
    GeneralVariables& rVariables,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double characteristic_element_size = CalculateElementSize();

    if (GetProperties().Has(DYNAMIC_VISCOSITY)) {

        const double constant1 = 16;
        const double constant2 = 4;
        const double delta_time = rCurrentProcessInfo[DELTA_TIME];
        const double& viscosity = GetProperties()[DYNAMIC_VISCOSITY]/ delta_time;
         GeometryType& r_geometry = GetGeometry();
        const unsigned int dimension = r_geometry.WorkingSpaceDimension();
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const unsigned int number_of_nodes = r_geometry.PointsNumber();
        VectorType veloc_mp = ZeroVector(3);
         double norm_veloc = 0.0;
        array_1d<double, 3> nodal_veloc = ZeroVector(3);
        for (unsigned int j = 0; j < number_of_nodes; ++j) {
            nodal_veloc = r_geometry[j].FastGetSolutionStepValue(VELOCITY, 0);
            for (unsigned int k = 0; k < dimension; ++k) {
                veloc_mp[k] += r_N(0, j) * nodal_veloc[k];
            }
        }
        for (unsigned int k = 0; k < dimension; ++k) {
            norm_veloc += std::pow(veloc_mp[k], 2);
        }
        norm_veloc = std::sqrt(norm_veloc); 
        rVariables.tau1 = constant1 * viscosity / std::pow(characteristic_element_size, 2)
            + constant2 * ((mMP.density * norm_veloc) / (characteristic_element_size * delta_time)); 
        rVariables.tau1 = std::pow(rVariables.tau1, -1);
        rVariables.tau2 = std::pow(characteristic_element_size, 2) / rVariables.tau1;
        //rVariables.tau2 = 2.0 * constant2 * viscosity;
    } else {
        const double constant1 = 1.0;
        const double constant2 = 1.0;
        const double shear_modulus = CalculateShearModulus();
        rVariables.tau1 = constant1 * characteristic_element_size * characteristic_element_size / (2.0 * shear_modulus);
        rVariables.tau2 = 2.0 * constant2 * shear_modulus;
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************


double MPMUpdatedLagrangianUPVMS::CalculateShearModulus() const
{
    KRATOS_TRY

    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO)) {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];

        return young_modulus / (2.0 * (1.0 + poisson_ratio));
    } else if (GetProperties().Has(DYNAMIC_VISCOSITY)) {
        return GetProperties()[DYNAMIC_VISCOSITY];
    }

    KRATOS_ERROR << "YOUNG_MODULUS and POISSON_RATIO or DYNAMIC_VISCOSITY are required "
        << "to calculate ASGS stabilization parameters." << std::endl;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

double MPMUpdatedLagrangianUPVMS::CalculateBulkModulus() const
{
    KRATOS_TRY

    if (GetProperties().Has(YOUNG_MODULUS) && GetProperties().Has(POISSON_RATIO)) {
        const double& young_modulus = GetProperties()[YOUNG_MODULUS];
        const double& poisson_ratio = GetProperties()[POISSON_RATIO];

        double bulk_modulus = young_modulus / (3.0 * (1.0 - 2.0 * poisson_ratio));
        const double tolerance = 10.0e-7;
        const bool is_incompressible_limit = (poisson_ratio > 0.5 - tolerance)
            || (poisson_ratio < -1.0 + tolerance);

        if (is_incompressible_limit || bulk_modulus != bulk_modulus) {
            bulk_modulus = 1.0e16;
        }

        return bulk_modulus;
    } else if (GetProperties().Has(DYNAMIC_VISCOSITY)) {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(BULK_MODULUS))
            << "BULK_MODULUS is required when DYNAMIC_VISCOSITY is used "
            << "to calculate ASGS stabilization parameters." << std::endl;
        return GetProperties()[BULK_MODULUS];
    }

    KRATOS_ERROR << "YOUNG_MODULUS and POISSON_RATIO or DYNAMIC_VISCOSITY are required "
        << "to calculate ASGS stabilization parameters." << std::endl;

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
                - rVariables.PressureGradient[jdim]
                + rVariables.DynamicRHS[jdim];

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
                - rVolumeForce[idime] / rIntegrationWeight
                + rVariables.DynamicRHS[idime];

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
    const Matrix& r_N = GetGeometry().ShapeFunctionsValues();

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
                    if (idim == jdim) {
                        rLeftHandSideMatrix(indexi + i, indexj + j) += rVariables.tau1
                            * rVariables.DynamicCoefficient
                            * r_N(0, j)
                            * row_displacement_pressure_gradient_projection(i)
                            * rIntegrationWeight;

                        rLeftHandSideMatrix(indexi + i, indexj + j) += rVariables.tau1
                            * rVariables.DynamicCoefficient
                            * r_N(0, j)
                            * row_displacement_deviatoric_pressure_gradient_projection(indexi)
                            * rIntegrationWeight;
                    }

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

                rLeftHandSideMatrix(index_p, index_up + k) -= rVariables.tau1
                    * volumetric_strain_linearization
                    * rVariables.DN_DX(i, k)
                    * rVariables.DynamicCoefficient
                    * r_N(0, j)
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

//***********************************************************************************
// Calculate Projections for OSGS stabilization
//***********************************************************************************

void MPMUpdatedLagrangianUPVMS::CalculateProjections(const ProcessInfo &rCurrentProcessInfo)
{

    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int number_of_nodes  = r_geometry.size();
    const unsigned int dimension        = r_geometry.WorkingSpaceDimension();
    const Vector& r_N       = row(GetGeometry().ShapeFunctionsValues(), 0);
    Vector MomentumRes      = ZeroVector(3);
    double ConservRes       = 0.0;
    VectorType momentum_rhs = ZeroVector(number_of_nodes*dimension);
    VectorType conserv_rhs  = ZeroVector(number_of_nodes);
    VectorType NodalArea    = ZeroVector(number_of_nodes);

     // Set the Initialize General Variables
    GeneralVariables Variables;
    this-> InitializeGeneralVariables(Variables,rCurrentProcessInfo);

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags &ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Set Variables
    this-> CalculateKinematics(Variables, rCurrentProcessInfo);
    this-> SetGeneralVariables(Variables, Values, r_N);
    this-> CalculatePressureVariables(Variables);
    const bool is_dynamic = rCurrentProcessInfo.Has(IS_DYNAMIC)
        ? rCurrentProcessInfo.GetValue(IS_DYNAMIC)
        : false;
    if (is_dynamic) {
        this-> CalculateDynamicStabilizationVariables(Variables, rCurrentProcessInfo);
    }


    mMP.volume = mMP.mass / mMP.density;
    mMP.body_force = this->ComputeMaterialPointBodyForce();
    Vector volume_force = ((mMP.volume_acceleration * mMP.mass ) + (mMP.body_force * mMP.mass)); ///mMP.volume 


    // Compute the Residual
    this->ComputeResidual(Variables,volume_force,MomentumRes,ConservRes);

    // Loop over the nodes for projecting
    for (unsigned int i = 0; i < number_of_nodes; i++)
    {
        double W = mMP.volume*r_N(i);
        unsigned int row = i*dimension;
        for (unsigned int d = 0; d < dimension; d++)
            momentum_rhs[row+d] += W*MomentumRes[d];
        NodalArea[i] += W;
        conserv_rhs[i] += W*ConservRes;
    }

     // Add to Global nodal variables
    for (SizeType i = 0; i < number_of_nodes; ++i)
      {
          r_geometry[i].SetLock(); // So it is safe to write in the node in OpenMP
          array_1d<double,3>& rMomValue = r_geometry[i].FastGetSolutionStepValue(RESPROJ_DISPL);
          unsigned int row = i*dimension;
          for (unsigned int d = 0; d < dimension; d++) rMomValue[d] += momentum_rhs[row + d];

          r_geometry[i].FastGetSolutionStepValue(RESPROJ_PRESS) += conserv_rhs[i];
          r_geometry[i].FastGetSolutionStepValue(NODAL_AREA) += NodalArea[i];
          r_geometry[i].UnSetLock(); // Free the node for other threads
      }

  
   
}


void MPMUpdatedLagrangianUPVMS::ComputeResidual(GeneralVariables& rVariables, Vector& rVolumeForce, Vector& rResidualU, double& rResidualP)
{

    GeometryType& r_geometry = this->GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    for (unsigned int k = 0; k < dimension; ++k) rResidualU[k] =  rVolumeForce[k] + rVariables.PressureGradient[k] +rVariables.DiscreteAcceleration[k]; //
    rResidualP =  -rVariables.PressureGP/rVariables.BulkModulus + (rVariables.detFT*rVariables.detF0-1);

    //std::cout<<mMP.acceleration[1]<<"\n";

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
