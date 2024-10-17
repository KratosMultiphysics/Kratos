// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_elements/timoshenko_beam_element_2D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void LinearTimoshenkoBeamElement2D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_5;
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoBeamElement2D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoBeamElement2D3N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoBeamElement2D3N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetShapeFunctionsValues(
    VectorType& rN,
    const double L,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);
    const double xi_square = std::pow(xi, 2);
    const double xi_cube   = std::pow(xi, 3);
    const double xi_quad   = std::pow(xi, 4);
    const double xi_quint  = std::pow(xi, 5);
    const double Phi_square = std::pow(Phi, 2);

    const double denom1 = 80.0 * Phi_square - 20.0 * Phi - 4.0;
    const double denom2 = 32.0 * Phi + 8.0;
    const double denom3 = 160.0 * Phi_square - 40.0 * Phi - 8.0;

    rN[0] = (-40.0 * std::pow(Phi, 2) - 10.0 * Phi) / denom1 * xi + (16.0 * Phi + 8.0) / denom2 * xi_square + (40.0 * Phi + 10.0) / denom3 * xi_cube + (-4.0) / denom2 * xi_quad + (-6.0) / denom3 * xi_quint;
    rN[1] = (-L * Phi) / denom1 * xi + L / denom2 * xi_square + L / denom3 * xi_cube + (-L) / denom2 * xi_quad + (2.0 * L * Phi - L) / denom3 * xi_quint;
    rN[2] = 1.0 + (-32.0 * Phi - 16.0) / denom2 * xi_square + 8.0 / denom2 * xi_quad;
    rN[3] = (-18.0 * L * Phi - 2.0 * L) / denom1 * xi + (40.0 * L * Phi + 8.0 * L) / denom3 * xi_cube + (-4.0 * L * Phi - 4.0 * L) / denom3 * xi_quint;
    rN[4] = (40.0 * std::pow(Phi, 2) + 10.0 * Phi) / denom1 * xi + (16.0 * Phi + 8.0) / denom2 * xi_square + (-40.0 * Phi - 10.0) / denom3 * xi_cube + (-4.0) / denom2 * xi_quad + 6.0 / denom3 * xi_quint;
    rN[5] = (-L * Phi) / denom1 * xi + (-L) / denom2 * xi_square + L / denom3 * xi_cube + L / denom2 * xi_quad + (2.0 * L * Phi - L) / denom3 * xi_quint;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double L,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);

    const double xi_square = std::pow(xi, 2);
    const double xi_cube   = std::pow(xi, 3);
    const double xi_quad   = std::pow(xi, 4);
    const double Phi_square = std::pow(Phi, 2);

    rN[0] = -30.0 * xi_quad / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 16.0 * xi_cube / (32.0 * Phi + 8.0) + 3.0 * xi_square * (40.0 * Phi + 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * xi * (16.0 * Phi + 8.0) / (32.0 * Phi + 8.0) + (-40.0 * Phi_square - 10.0 * Phi) / (80.0 * Phi_square - 20.0 * Phi - 4.0);
    rN[1] = -L * Phi / (80.0 * Phi_square - 20.0 * Phi - 4.0) - 4.0 * L * xi_cube / (32.0 * Phi + 8.0) + 3.0 * L * xi_square / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * L * xi / (32.0 * Phi + 8.0) + 5 * xi_quad * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[2] = 32.0 * xi_cube / (32.0 * Phi + 8.0) + 2.0 * xi * (-32.0 * Phi - 16.0) / (32.0 * Phi + 8.0);
    rN[3] = 5.0 * xi_quad * (-4.0 * L * Phi - 4.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 3.0 * xi_square * (40.0 * L * Phi + 8.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + (-18.0 * L * Phi - 2.0 * L) / (80.0 * Phi_square - 20.0 * Phi - 4.0);
    rN[4] = 30.0 * xi_quad / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 16.0 * xi_cube / (32.0 * Phi + 8.0) + 3.0 * xi_square * (-40.0 * Phi - 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * xi * (16.0 * Phi + 8.0) / (32.0 * Phi + 8.0) + (40.0 * Phi_square + 10.0 * Phi) / (80.0 * Phi_square - 20.0 * Phi - 4.0);
    rN[5] = -L * Phi / (80.0 * Phi_square - 20.0 * Phi - 4.0) + 4.0 * L * xi_cube / (32.0 * Phi + 8.0) + 3.0 * L * xi_square / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 2.0 * L * xi / (32.0 * Phi + 8.0) + 5.0 * xi_quad * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    
    rN *= 2.0 / L;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetSecondDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double L,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);

    const double xi_square = std::pow(xi, 2);
    const double xi_cube   = std::pow(xi, 3);
    const double Phi_square = std::pow(Phi, 2);

    rN[0] = -120.0 * xi_cube / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 48.0 * xi_square / (32.0 * Phi + 8.0) + 6.0 * xi * (40.0 * Phi + 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * (16.0 * Phi + 8.0) / (32.0 * Phi + 8.0);
    rN[1] = -12.0 * L * xi_square / (32.0 * Phi + 8.0) + 6.0 * L * xi / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * L / (32.0 * Phi + 8.0) + 20.0 * xi_cube * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[2] = 96.0 * xi_square / (32.0 * Phi + 8.0) + 2.0 * (-32.0 * Phi - 16.0) / (32.0 * Phi + 8.0);
    rN[3] = 20.0 * xi_cube * (-4.0 * L * Phi - 4.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 6.0 * xi * (40.0 * L * Phi + 8.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[4] = 120.0 * xi_cube / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 48.0 * xi_square / (32.0 * Phi + 8.0) + 6.0 * xi * (-40.0 * Phi - 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 2.0 * (16.0 * Phi + 8.0) / (32.0 * Phi + 8.0);
    rN[5] = 12.0 * L * xi_square / (32.0 * Phi + 8.0) + 6.0 * L * xi / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 2.0 * L / (32.0 * Phi + 8.0) + 20.0 * xi_cube * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);

    rN *= std::pow(2.0 / L, 2);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetThirdDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double L,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);

    const double xi_square = std::pow(xi, 2);
    const double Phi_square = std::pow(Phi, 2);

    rN[0] = -360.0 * xi_square / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 96.0 * xi / (32.0 * Phi + 8.0) + 6.0 * (40.0 * Phi + 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[1] = -24.0 * L * xi / (32.0 * Phi + 8.0) + 6.0 * L / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 60.0 * xi_square * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[2] = 192.0 * xi / (32.0 * Phi + 8.0);
    rN[3] = 60.0 * xi_square * (-4.0 * L * Phi - 4.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 6.0 * (40.0 * L * Phi + 8.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[4] = 360.0 * xi_square / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 96.0 * xi / (32.0 * Phi + 8.0) + 6.0 * (-40.0 * Phi - 10.0) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[5] = 24.0 * L * xi / (32.0 * Phi + 8.0) + 6.0 * L / (160.0 * Phi_square - 40.0 * Phi - 8.0) + 60.0 * xi_square * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN *= std::pow(2.0 / L, 3);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetFourthDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double L,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);

    const double Phi_square = std::pow(Phi, 2);

    rN[0] = -720.0 * xi / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 96.0 / (32.0 * Phi + 8.0);
    rN[1] = -24.0 * L / (32.0 * Phi + 8.0) + 120.0 * xi * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[2] = 192.0 / (32.0 * Phi + 8.0);
    rN[3] = 120.0 * xi * (-4.0 * L * Phi - 4.0 * L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);
    rN[4] = 720.0 * xi / (160.0 * Phi_square - 40.0 * Phi - 8.0) - 96.0 / (32.0 * Phi + 8.0);
    rN[5] = 24.0 * L / (32.0 * Phi + 8.0) + 120.0 * xi * (2.0 * L * Phi - L) / (160.0 * Phi_square - 40.0 * Phi - 8.0);

    rN *= std::pow(2.0 / L, 4);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);
    VectorType N_derivative(6), N_third_derivative(6);
    GetFirstDerivativesShapeFunctionsValues(N_derivative, Length, Phi, xi);
    GetThirdDerivativesShapeFunctionsValues(N_third_derivative, Length, Phi, xi);
    // v' + (Phi * L^2 / 12) * v'''
    noalias(rN) = N_derivative + Phi * std::pow(Length, 2) / 12.0 * N_third_derivative;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetFirstDerivativesNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 6)
        rN.resize(6, false);
    VectorType N_second_derivative(6), N_fourth_derivative(6);
    GetSecondDerivativesShapeFunctionsValues(N_second_derivative, Length, Phi, xi);
    GetFourthDerivativesShapeFunctionsValues(N_fourth_derivative, Length, Phi, xi);
    // v'' + (Phi * L^2 / 12) * v''''
    noalias(rN) = N_second_derivative + Phi * std::pow(Length, 2) / 12.0 * N_fourth_derivative;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 3)
        rN.resize(3, false);
    rN[0] = 0.5 * xi * (xi - 1.0);
    rN[1] = (1.0 - std::pow(xi, 2));
    rN[2] = 0.5 * xi * (xi + 1.0);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetFirstDerivativesNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    ) const
{
    if (rN.size() != 3)
        rN.resize(3, false);
    rN[0] = xi - 0.5;
    rN[1] = -2.0 * xi;
    rN[2] = xi + 0.5;
    rN *= 2.0 / Length;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::GetNodalValuesVector(VectorType& rNodalValues) const
{
    if (rNodalValues.size() != 9)
        rNodalValues.resize(9, false);
    const auto &r_geom = GetGeometry();

    const double angle = GetAngle();

    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, 3, 3> T;
        BoundedVector<double, 9> global_values;
        BoundedMatrix<double, 9, 9> global_size_T;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

        const auto& r_displ_0 = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_displ_1 = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_displ_2 = r_geom[2].FastGetSolutionStepValue(DISPLACEMENT);

        global_values[0] = r_displ_0[0];
        global_values[1] = r_displ_0[1];
        global_values[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);

        global_values[3] = r_displ_1[0];
        global_values[4] = r_displ_1[1];
        global_values[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);

        global_values[6] = r_displ_2[0];
        global_values[7] = r_displ_2[1];
        global_values[8] = r_geom[2].FastGetSolutionStepValue(ROTATION_Z);

        // We rotate to local axes
        noalias(rNodalValues) = prod(trans(global_size_T), global_values);

    } else {
        const auto& r_displ_0 = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_displ_1 = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_displ_2 = r_geom[2].FastGetSolutionStepValue(DISPLACEMENT);

        rNodalValues[0] = r_displ_0[0];
        rNodalValues[1] = r_displ_0[1];
        rNodalValues[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);
        rNodalValues[3] = r_displ_1[0];
        rNodalValues[4] = r_displ_1[1];
        rNodalValues[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);
        rNodalValues[6] = r_displ_2[0];
        rNodalValues[7] = r_displ_2[1];
        rNodalValues[8] = r_geom[2].FastGetSolutionStepValue(ROTATION_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement2D3N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_u0_derivatives(3);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_u0_derivatives, Length, Phi, xi);
    return N_u0_derivatives[0] * rNodalValues[0] + N_u0_derivatives[2] * rNodalValues[3] + N_u0_derivatives[1] * rNodalValues[6];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement2D3N::CalculateShearStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_derivatives(6), N_theta(6);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const VectorType N_s = N_derivatives - N_theta;
    return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[2] + N_s[4] * rNodalValues[4] +
           N_s[5] * rNodalValues[5] + N_s[2] * rNodalValues[7] + N_s[3] * rNodalValues[8];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement2D3N::CalculateBendingCurvature(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_theta_derivatives(6);
    GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
           N_theta_derivatives[4] * rNodalValues[4] + N_theta_derivatives[5] * rNodalValues[5] +
           N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[8];
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
)
{
    const double angle = GetAngle();

    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, 3, 3> T, Tt;
        BoundedMatrix<double, 9, 9> global_size_T, aux_product;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);
        noalias(aux_product) = prod(rLHS, trans(global_size_T));
        noalias(rLHS) = prod(global_size_T, aux_product);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{
    const double angle = GetAngle();
    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, 3, 3> T;
        BoundedMatrix<double, 9, 9> global_size_T;
        BoundedVector<double, 9> local_rhs;
        noalias(local_rhs) = rRHS;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

        noalias(rRHS) = prod(global_size_T, local_rhs);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{
    const double angle = GetAngle();
    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, 3, 3> T;
        BoundedMatrix<double, 9, 9> global_size_T, aux_product;
        BoundedVector<double, 9> local_rhs;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForBeam(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NBeam(T, global_size_T);

        noalias(local_rhs) = rRHS;
        noalias(rRHS) = prod(global_size_T, local_rhs);

        noalias(aux_product) = prod(rLHS, trans(global_size_T));
        noalias(rLHS) = prod(global_size_T, aux_product);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearTimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement2D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearTimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
