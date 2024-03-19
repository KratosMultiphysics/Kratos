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

void TimoshenkoBeamElement2D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
        BaseType::InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer TimoshenkoBeamElement2D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TimoshenkoBeamElement2D3N::Pointer p_new_elem = Kratos::make_intrusive<TimoshenkoBeamElement2D3N>
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

void TimoshenkoBeamElement2D3N::GetShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // const double one_plus_phi = 1.0 + Phi;
    // const double xi_square = xi * xi;
    // rN[0] = (xi - 1.0) * (xi + xi_square - 2.0 * one_plus_phi) / (4.0 * one_plus_phi);
    // rN[1] = (1.0 - xi_square) * (1.0 - xi + Phi) * Length / (8.0 * one_plus_phi);
    // rN[2] = (1.0 + xi) * (xi - xi_square + 2.0 * one_plus_phi) / (4.0 * one_plus_phi);
    // rN[3] = (xi_square - 1.0) * (1.0 + xi + Phi) * Length / (8.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // const double one_plus_phi = 1.0 + Phi;
    // const double xi_square = xi * xi;
    // rN[0] = (-6.0 + 6.0 * xi_square - 4.0 * Phi) / (4.0 * one_plus_phi * Length);
    // rN[1] = (-1.0 + 3.0 * xi_square - 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi);
    // rN[2] = (6.0 - 6.0 * xi_square + 4.0 * Phi) / (4.0 * one_plus_phi * Length);
    // rN[3] = (-1.0 + 3.0 * xi_square + 2.0 * xi * one_plus_phi) / (4.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetSecondDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // const double one_plus_phi = 1.0 + Phi;
    // const double L_square = std::pow(Length, 2);
    // rN[0] = 6.0 * xi / (one_plus_phi * L_square);
    // rN[1] = (-1.0 + 3.0 * xi - Phi) / (one_plus_phi * Length);
    // rN[2] = -6.0 * xi / (one_plus_phi * L_square);
    // rN[3] = (1.0 + 3.0 * xi + Phi) / (one_plus_phi * Length);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetThirdDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // const double one_plus_phi = 1.0 + Phi;
    // const double L_square = std::pow(Length, 2);
    // const double L_cube   = std::pow(Length, 3);
    // rN[0] = 12.0  / (one_plus_phi * L_cube);
    // rN[1] = 6.0   / (one_plus_phi * L_square);
    // rN[2] = -12.0 / (one_plus_phi * L_cube);
    // rN[3] = 6.0   / (one_plus_phi * L_square);
}

// /***********************************************************************************/
// /***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // const double one_plus_phi = 1.0 + Phi;
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // rN[0] = (3.0 * xi * xi - 3.0) / (2.0 * one_plus_phi * Length);
    // rN[1] = (xi - 1.0) * (1.0 + 3.0 * xi - 2.0 * Phi) / (4.0 * one_plus_phi);
    // rN[2] = (3.0 - 3.0 * xi * xi) / (2 * one_plus_phi * Length);
    // rN[3] = (1.0 + xi) * (3.0 * xi - 1.0 + 2.0 * Phi) / (4.0 * one_plus_phi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetFirstDerivativesNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // const double one_plus_phi = 1.0 + Phi;
    // if (rN.size() != 4)
    //     rN.resize(4, false);
    // rN[0] = 3.0 * xi / (one_plus_phi * Length);
    // rN[1] = (-0.5 * Phi + 1.5 * xi - 0.5) / (one_plus_phi);
    // rN[2] = (-3.0 * xi) / (Length * Phi + Length);
    // rN[3] = (0.5 * Phi + 1.5 * xi + 0.5) / (one_plus_phi);
    // rN *= (2.0 / Length);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 2)
    //     rN.resize(2, false);
    // rN[0] = 0.5 * (1.0 - xi);
    // rN[1] = 0.5 * (1.0 + xi);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetFirstDerivativesNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi
    )
{
    // if (rN.size() != 2)
    //     rN.resize(2, false);
    // const double inverse_l = 1.0 / Length;
    // rN[0] = -inverse_l;
    // rN[1] =  inverse_l;
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::GetNodalValuesVector(VectorType& rNodalValues)
{
    // if (rNodalValues.size() != 6)
    //     rNodalValues.resize(6, false);
    // const auto &r_geom = GetGeometry();

    // const double angle = StructuralMechanicsElementUtilities::
    //     GetReferenceRotationAngle2D2NBeam(GetGeometry());

    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedVector<double, 6> global_values;
    //     BoundedMatrix<double, 6, 6> global_size_T;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixFor2D2NBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(T, global_size_T);

    //     global_values[0] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    //     global_values[1] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    //     global_values[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);
    //     global_values[3] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    //     global_values[4] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    //     global_values[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);

    //     // We rotate to local axes
    //     noalias(rNodalValues) = prod(trans(global_size_T), global_values);

    // } else {
    //     rNodalValues[0] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_X);
    //     rNodalValues[1] = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
    //     rNodalValues[2] = r_geom[0].FastGetSolutionStepValue(ROTATION_Z);
    //     rNodalValues[3] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_X);
    //     rNodalValues[4] = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT_Y);
    //     rNodalValues[5] = r_geom[1].FastGetSolutionStepValue(ROTATION_Z);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

double TimoshenkoBeamElement2D3N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    // VectorType N_u0_derivatives(2);
    // GetFirstDerivativesNu0ShapeFunctionsValues(N_u0_derivatives, Length, Phi, xi);
    // return N_u0_derivatives[0] * rNodalValues[0] + N_u0_derivatives[1] * rNodalValues[3];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double TimoshenkoBeamElement2D3N::CalculateShearStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    // VectorType N_derivatives(4), N_theta(4);
    // GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    // GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    // const VectorType N_s = N_derivatives - N_theta;
    // return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[2] + N_s[2] * rNodalValues[4] + 
    //        N_s[3] * rNodalValues[5];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

double TimoshenkoBeamElement2D3N::CalculateBendingCurvature(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    )
{
    // VectorType N_theta_derivatives(4);
    // GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    // return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[2] +
    //        N_theta_derivatives[2] * rNodalValues[4] + N_theta_derivatives[3] * rNodalValues[5];
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
)
{
    // const double angle = StructuralMechanicsElementUtilities::
    //     GetReferenceRotationAngle2D2NBeam(GetGeometry());

    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T, Tt;
    //     BoundedMatrix<double, 6, 6> global_size_T, aux_product;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixFor2D2NBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(T, global_size_T);
    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{
    // const double angle = StructuralMechanicsElementUtilities::
    //     GetReferenceRotationAngle2D2NBeam(GetGeometry());
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 6, 6> global_size_T;
    //     BoundedVector<double, 6> local_rhs;
    //     noalias(local_rhs) = rRHS;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixFor2D2NBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(T, global_size_T);

    //     noalias(rRHS) = prod(global_size_T, local_rhs);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
)
{
    // const double angle = StructuralMechanicsElementUtilities::
    //     GetReferenceRotationAngle2D2NBeam(GetGeometry());
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, 3, 3> T;
    //     BoundedMatrix<double, 6, 6> global_size_T, aux_product;
    //     BoundedVector<double, 6> local_rhs;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixFor2D2NBeam(T, angle);
    //     StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(T, global_size_T);

    //     noalias(local_rhs) = rRHS;
    //     noalias(rRHS) = prod(global_size_T, local_rhs);

    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElement2D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TimoshenkoBeamElement2D2N);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
