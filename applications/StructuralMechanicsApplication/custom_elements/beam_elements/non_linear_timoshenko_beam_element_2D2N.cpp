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

// Project includes
#include "non_linear_timoshenko_beam_element_2D2N.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer NonLinearTimoshenkoBeamElement2D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    NonLinearTimoshenkoBeamElement2D2N::Pointer p_new_elem = Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[1] = simple_N[0];
    rN[4] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetFirstDerivativesNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[1] = simple_N[0];
    rN[4] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[2] = simple_N[0];
    rN[5] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesNThetaShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetFirstDerivativesNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[2] = simple_N[0];
    rN[5] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[0] = simple_N[0];
    rN[3] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::GetFirstDerivativesNu0ShapeFunctionsValues(
    VectorType& rN,
    const double Length,
    const double Phi,
    const double xi) const
{
    VectorType simple_N(2);
    BaseType::GetFirstDerivativesNu0ShapeFunctionsValues(simple_N, Length, Phi, xi);

    rN.clear();
    rN[0] = simple_N[0];
    rN[3] = simple_N[1];
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, true, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    VectorType dummy_rhs;
    CalculateAll(rLeftHandSideMatrix, dummy_rhs, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    MatrixType dummy_lhs;
    CalculateAll(dummy_lhs, rRightHandSideVector, rCurrentProcessInfo, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY

    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();


    if (ComputeLHS) {
        if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) {
            rLHS.resize(mat_size, mat_size, false);
        }
        noalias(rLHS) = ZeroMatrix(mat_size, mat_size);
    }

    if (ComputeRHS) {
        if (rRHS.size() != mat_size) {
            rRHS.resize(mat_size, false);
        }
        noalias(rRHS) = ZeroVector(mat_size);
    }

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, ComputeRHS);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, ComputeLHS);

    const double length = CalculateLength(); // Reference length
    const double Phi = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
    const double J = 0.5 * length;
    const double area = GetCrossArea();

    // Let's initialize the cl values
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

    VectorType Ntheta(mat_size), dNtheta(mat_size);
    VectorType Nu(mat_size), dNu(mat_size);
    VectorType Nv(mat_size), dNv(mat_size);

    MatrixType B(strain_size, mat_size);
    double theta, du, dv;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J;

        // Let's fill the shape functions and their derivatives
        GetNu0ShapeFunctionsValues(Nu, length, Phi, xi); // u
        GetFirstDerivativesNu0ShapeFunctionsValues(dNu, length, Phi, xi);

        GetNThetaShapeFunctionsValues(Ntheta, length, Phi, xi); // rotation
        GetFirstDerivativesNThetaShapeFunctionsValues(dNtheta, length, Phi, xi);

        GetShapeFunctionsValues(Nv, length, Phi, xi); // v
        GetFirstDerivativesShapeFunctionsValues(dNv, length, Phi, xi);

        theta = inner_prod(Ntheta, nodal_values);
        du = inner_prod(dNu, nodal_values);
        dv = inner_prod(dNv, nodal_values);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);

        strain_vector[0] = (1.0 + du) * cos_theta + dv * sin_theta - 1.0; // axial strain
        strain_vector[1] = inner_prod(dNtheta, nodal_values); // curvature
        strain_vector[2] = dv * cos_theta - sin_theta * (1.0 + du); // shear strain

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);

        const double N = stress_vector[0]; // Axial force
        const double M = stress_vector[1]; // Bending moment
        const double V = stress_vector[2]; // Shear force

        // axial strain increment
        noalias(row(B, 0)) = -sin_theta * Ntheta + cos_theta * dNu - sin_theta * du * Ntheta + sin_theta * dNv + dv * cos_theta * Ntheta;
        // curvature increment
        noalias(row(B, 1)) = dNtheta;
        // shear strain increment
        noalias(row(B, 2)) = cos_theta * dNv - sin_theta * dv * Ntheta - cos_theta * Ntheta - sin_theta * dNu - du * cos_theta * Ntheta;

        if (ComputeRHS) {
            // Internal forces
            noalias(rRHS) -= prod(trans(B), stress_vector) * jacobian_weight;
            // Body forces
            noalias(rRHS) += Nu * local_body_forces[0] * jacobian_weight * area;
            noalias(rRHS) += Nv * local_body_forces[1] * jacobian_weight * area;
        }

        if (ComputeLHS) {
            MatrixType BtD(mat_size, strain_size);
            noalias(BtD) = prod(trans(B), constitutive_matrix);
            // Material K
            noalias(rLHS) += prod(BtD, B) * jacobian_weight;

            // Geometric K, derivative of the B matrix with respect to the nodal displacements
            MatrixType dB_axial(mat_size, mat_size); // dB_axial / da(e)
            MatrixType dB_shear(mat_size, mat_size); // dB_shear / da(e)

            noalias(dB_axial) = -cos_theta * outer_prod(Ntheta, Ntheta) - sin_theta * outer_prod(dNu, Ntheta) - outer_prod(Ntheta, sin_theta * dNu + du * cos_theta * Ntheta) +
                                cos_theta * outer_prod(dNv, Ntheta) + outer_prod(Ntheta, dNv * cos_theta - dv * sin_theta * Ntheta);
            noalias(dB_shear) = -sin_theta * outer_prod(dNv, Ntheta) - outer_prod(Ntheta, sin_theta * dNv + cos_theta * dv * Ntheta) - cos_theta * outer_prod(dNu, Ntheta) -
                                outer_prod(Ntheta, -(1.0 + du) * sin_theta * Ntheta + cos_theta * dNu);
            noalias(rLHS) += (N * dB_axial + V * dB_shear) * jacobian_weight;
        }

    }
    if (ComputeLHS && ComputeRHS) {
        RotateAll(rLHS, rRHS, r_geometry);
    } else if (ComputeLHS) {
        RotateLHS(rLHS, r_geometry);
    } else if (ComputeRHS) {
        RotateRHS(rRHS, r_geometry);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void NonLinearTimoshenkoBeamElement2D2N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType mat_size = GetDoFsPerNode() * GetGeometry().size();
    rOutput.resize(r_integration_points.size());
    const auto &r_props = GetProperties();

    if (rVariable == AXIAL_FORCE || rVariable == SHEAR_FORCE || rVariable == BENDING_MOMENT || rVariable == AXIAL_STRAIN || rVariable == SHEAR_STRAIN || rVariable == BENDING_STRAIN) {

        const auto &r_geometry = GetGeometry();

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        const double length = CalculateLength();
        const double Phi    = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
        VectorType strain_vector(strain_size), stress_vector(strain_size);
        StructuralMechanicsElementUtilities::InitializeConstitutiveLawValuesForStressCalculation(cl_values, strain_vector, stress_vector);

        VectorType nodal_values(mat_size);
        GetNodalValuesVector(nodal_values);

        VectorType Ntheta(mat_size), dNtheta(mat_size);
        VectorType Nu(mat_size), dNu(mat_size);
        VectorType Nv(mat_size), dNv(mat_size);
        double theta, du, dv;

        // Loop over the integration points
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi = r_integration_points[IP].X();

            // Let's fill the shape functions and their derivatives
            GetNu0ShapeFunctionsValues(Nu, length, Phi, xi); // u
            GetFirstDerivativesNu0ShapeFunctionsValues(dNu, length, Phi, xi);

            GetNThetaShapeFunctionsValues(Ntheta, length, Phi, xi); // rotation
            GetFirstDerivativesNThetaShapeFunctionsValues(dNtheta, length, Phi, xi);

            GetShapeFunctionsValues(Nv, length, Phi, xi); // v
            GetFirstDerivativesShapeFunctionsValues(dNv, length, Phi, xi);

            theta = inner_prod(Ntheta, nodal_values);
            du = inner_prod(dNu, nodal_values);
            dv = inner_prod(dNv, nodal_values);

            const double sin_theta = std::sin(theta);
            const double cos_theta = std::cos(theta);

            strain_vector[0] = (1.0 + du) * cos_theta + dv * sin_theta - 1.0; // axial strain

            if (rVariable == AXIAL_STRAIN) {
                rOutput[IP] = strain_vector[0];
                return;
            }

            strain_vector[1] = inner_prod(dNtheta, nodal_values); // curvature


            if (rVariable == BENDING_STRAIN) {
                rOutput[IP] = strain_vector[1];
                return;
            }

            strain_vector[2] = dv * cos_theta - sin_theta * (1.0 + du); // shear strain

            if (rVariable == SHEAR_STRAIN) {
                rOutput[IP] = strain_vector[2];
                return;
            }

            mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
            const Vector &r_generalized_stresses = cl_values.GetStressVector();

            if (rVariable == AXIAL_FORCE) {
                rOutput[IP] = r_generalized_stresses[0];
            } else if (rVariable == BENDING_MOMENT) {
                rOutput[IP] = r_generalized_stresses[1];
            } else if (rVariable == SHEAR_FORCE) {
                rOutput[IP] = r_generalized_stresses[2];
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

} // namespace Kratos
