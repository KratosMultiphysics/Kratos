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
#include "timoshenko_beam_element_3D2N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoBeamElement3D2N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoBeamElement3D2N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoBeamElement3D2N>
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

array_1d<double, 3> LinearTimoshenkoBeamElement3D2N::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    BoundedMatrix<double, 3, 3> T;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(GetGeometry());
    array_1d<double, 3> local_body_force = prod(T, body_force);
    return local_body_force;
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 3, 3> LinearTimoshenkoBeamElement3D2N::GetConsistentFrenetSerretMatrix3D(
    const GeometryType& rGeometry
    ) const
{
    BoundedMatrix<double, 3, 3> T;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(rGeometry);
    return T;
}

/***********************************************************************************/
/***********************************************************************************/


void LinearTimoshenkoBeamElement3D2N::GetNodalValuesVector(
    VectorType& rNodalValues
    ) const
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const SizeType num_nodes = r_geom.size();
    const SizeType global_size = GetDoFsPerNode() * num_nodes;

    if (rNodalValues.size() != global_size)
        rNodalValues.resize(global_size, false);

    BoundedMatrix<double, 3, 3> T;
    // From global to local, need to transpose the matrix
    noalias(T) = trans(GetConsistentFrenetSerretMatrix3D(r_geom));

    const auto& r_displ_0    = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_0 = r_geom[0].FastGetSolutionStepValue(ROTATION);
    const auto& r_displ_1    = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_1 = r_geom[1].FastGetSolutionStepValue(ROTATION);

    // Here we rotate the vectors to local axes
    const VectorType& r_local_displ_0 = prod(T, r_displ_0);
    const VectorType& r_local_displ_1 = prod(T, r_displ_1);

    // Due to dextrogyr axes system, the y rot is inverted
    T(1, 0) *= -1.0;
    T(1, 1) *= -1.0;
    T(1, 2) *= -1.0;

    const VectorType& r_local_rot_0 = prod(T, r_rotation_0);
    const VectorType& r_local_rot_1 = prod(T, r_rotation_1);

    rNodalValues[0] = r_local_displ_0[0];
    rNodalValues[1] = r_local_displ_0[1];
    rNodalValues[2] = r_local_displ_0[2];

    rNodalValues[3] = r_local_rot_0[0];
    rNodalValues[4] = r_local_rot_0[1];
    rNodalValues[5] = r_local_rot_0[2];

    rNodalValues[6] = r_local_displ_1[0];
    rNodalValues[7] = r_local_displ_1[1];
    rNodalValues[8] = r_local_displ_1[2];

    rNodalValues[9]  = r_local_rot_1[0];
    rNodalValues[10] = r_local_rot_1[1];
    rNodalValues[11] = r_local_rot_1[2];

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateAxialStrain(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_u_derivatives(2);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, Length, Phi, xi);
    return N_u_derivatives[0] * rNodalValues[0] + N_u_derivatives[1] * rNodalValues[6];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_derivatives(4), N_theta(4);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const VectorType N_s = N_derivatives - N_theta;
    return N_s[0] * rNodalValues[1] + N_s[1] * rNodalValues[5] + N_s[2] * rNodalValues[7] + N_s[3] * rNodalValues[11];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateShearStrainXZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_derivatives(4), N_theta(4);
    GetFirstDerivativesShapeFunctionsValues(N_derivatives, Length, Phi, xi);
    GetNThetaShapeFunctionsValues(N_theta, Length, Phi, xi);
    const VectorType N_s = N_derivatives - N_theta;
    return N_s[0] * rNodalValues[2] + N_s[1] * rNodalValues[4] + N_s[2] * rNodalValues[8] + N_s[3] * rNodalValues[10];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureX(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_theta_x_derivatives(2);
    GetFirstDerivativesNu0ShapeFunctionsValues(N_theta_x_derivatives, Length, Phi, xi);
    return N_theta_x_derivatives[0] * rNodalValues[3] + N_theta_x_derivatives[1] * rNodalValues[9];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureY(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_theta_derivatives(4);
    GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    return N_theta_derivatives[0] * rNodalValues[2] + N_theta_derivatives[1] * rNodalValues[4] +
           N_theta_derivatives[2] * rNodalValues[8] + N_theta_derivatives[3] * rNodalValues[10];
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoBeamElement3D2N::CalculateBendingCurvatureZ(
    const double Length,
    const double Phi,
    const double xi,
    const VectorType& rNodalValues
    ) const
{
    VectorType N_theta_derivatives(4);
    GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, Length, Phi, xi);
    return N_theta_derivatives[0] * rNodalValues[1] + N_theta_derivatives[1] * rNodalValues[5] +
           N_theta_derivatives[2] * rNodalValues[7] + N_theta_derivatives[3] * rNodalValues[11];
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::AssembleGlobalRotationMatrix(
    const BoundedMatrix<double, 3, 3>& rT,
    BoundedMatrix<double, 12, 12>& rGlobalT
)
{
    const SizeType nnodes = GetGeometry().size();
    rGlobalT.clear();

    for (IndexType block = 0; block < 2*nnodes; ++block) {
        for (IndexType i = 0; i < rT.size1(); ++i) {
            for (IndexType j = 0; j < rT.size2(); ++j) {
                if (block == 1 || block == 3) { // blocks affecting the rotations
                    if (j == 1) { // the y rotation
                        rGlobalT(3 * block + i, 3 * block + j) = - rT(i, j);
                    } else {
                        rGlobalT(3 * block + i, 3 * block + j) = rT(i, j);
                    }
                } else {
                    rGlobalT(3 * block + i, 3 * block + j) = rT(i, j);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


void LinearTimoshenkoBeamElement3D2N::RotateLHS(
    MatrixType& rLHS,
    const GeometryType& rGeometry
    )
{
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T, aux_product;
    noalias(T) = GetConsistentFrenetSerretMatrix3D(rGeometry);
    AssembleGlobalRotationMatrix(T, global_size_T);

    noalias(aux_product) = prod(rLHS, trans(global_size_T));
    noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateRHS(
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T;
    BoundedVector<double, 12> local_rhs;
    noalias(local_rhs) = rRHS;
    noalias(T) = GetConsistentFrenetSerretMatrix3D(rGeometry);
    AssembleGlobalRotationMatrix(T, global_size_T);
    noalias(rRHS) = prod(global_size_T, local_rhs);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const GeometryType& rGeometry
    )
{
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T, aux_product;
    noalias(T) = GetConsistentFrenetSerretMatrix3D(rGeometry);
    AssembleGlobalRotationMatrix(T, global_size_T);

    BoundedVector<double, 12> local_rhs;
    noalias(local_rhs) = rRHS;
    noalias(rRHS) = prod(global_size_T, local_rhs);

    noalias(aux_product) = prod(rLHS, trans(global_size_T));
    noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateGeneralizedStrainsVector(
    VectorType& rStrain,
    const double Length,
    const double Phi,
    const double xi,
    const VectorType &rNodalValues
    ) const
{
    if (rStrain.size() != 6)
        rStrain.resize(6, false);

    const auto& r_props = GetProperties();
    const double Phi_rot_z  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, Length);
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, Length, 1);

    rStrain[0] = CalculateAxialStrain(Length, 0.0, xi, rNodalValues);
    rStrain[1] = CalculateBendingCurvatureX(Length, 0.0, xi, rNodalValues);
    rStrain[2] = CalculateBendingCurvatureY(Length, Phi_rot_y, xi, rNodalValues);
    rStrain[3] = CalculateBendingCurvatureZ(Length, Phi_rot_z, xi, rNodalValues);
    rStrain[4] = CalculateShearStrainXY(Length, Phi_rot_z, xi, rNodalValues);
    rStrain[5] = CalculateShearStrainXZ(Length, Phi_rot_y, xi, rNodalValues);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType local_trans_deflection_size = mat_size - 4 * number_of_nodes;

    if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) {
        rLHS.resize(mat_size, mat_size, false);
    }
    rLHS.clear();

    if (rRHS.size() != mat_size) {
        rRHS.resize(mat_size, false);
    }
    rRHS.clear();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double length = CalculateLength();
    const double Phi_rot_z  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length, 1);
    const double J      = 0.5 * length;
    const double area   = GetCrossArea();

    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

    VectorType global_size_N(mat_size);
    VectorType N_u_derivatives(number_of_nodes);
    VectorType N_theta_derivatives(local_trans_deflection_size);
    VectorType N_theta(local_trans_deflection_size);
    VectorType N_derivatives(local_trans_deflection_size);
    VectorType N_u(number_of_nodes);
    VectorType N_shape(local_trans_deflection_size);
    VectorType N_s(local_trans_deflection_size);

    // Loop over the integration points (IP)
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        global_size_N.clear();
        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J;

        CalculateGeneralizedStrainsVector(strain_vector, length, 0.0, xi, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N  = r_generalized_stresses[0];
        const double Mx = r_generalized_stresses[1];
        const double My = r_generalized_stresses[2];
        const double Mz = r_generalized_stresses[3];
        const double Vy = r_generalized_stresses[4];
        const double Vz = r_generalized_stresses[5];

        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl       = r_constitutive_matrix(0, 0);
        const double dMx_dkappa_x = r_constitutive_matrix(1, 1);
        const double dMy_dkappa_y = r_constitutive_matrix(2, 2);
        const double dMz_dkappa_z = r_constitutive_matrix(3, 3);
        const double dVy_dgamma_xy = r_constitutive_matrix(4, 4);
        const double dVz_dgamma_xz = r_constitutive_matrix(5, 5);

        // Axial DoFs shape functions (u and theta_x)
        GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, length, 0.0, xi);
        GetNu0ShapeFunctionsValues(N_u, length, 0.0, xi);

        // Axial contributions
        GlobalSizeAxialVector(global_size_N, N_u_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dN_dEl * jacobian_weight;
        noalias(rRHS) -= global_size_N * N * jacobian_weight;

        // Torsional contributions
        GlobalSizeVectorAxialRotation(global_size_N, N_u_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMx_dkappa_x * jacobian_weight;
        noalias(rRHS) -= global_size_N * Mx * jacobian_weight;

        // External body forces in u
        GlobalSizeAxialVector(global_size_N, N_u);
        noalias(rRHS) += global_size_N * local_body_forces[0] * jacobian_weight * area;

        // Transverse DoFs shape functions {v, theta_z}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_z, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_z, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_z, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_z, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in z contributions
        GlobalSizeVectorTransversalY(global_size_N, N_theta_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMz_dkappa_z * jacobian_weight;
        noalias(rRHS) -= global_size_N * Mz * jacobian_weight;

        // Shear XY contributions
        GlobalSizeVectorTransversalY(global_size_N, N_s);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dVy_dgamma_xy * jacobian_weight;
        noalias(rRHS) -= global_size_N * Vy * jacobian_weight;

		// External body forces in v
        GlobalSizeVectorTransversalY(global_size_N, N_shape);
        noalias(rRHS) += global_size_N * local_body_forces[1] * jacobian_weight * area;

        // Transverse DoFs shape functions {w, theta_y}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_y, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_y, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_y, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_y, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in z contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_theta_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMy_dkappa_y * jacobian_weight;
        noalias(rRHS) -= global_size_N * My * jacobian_weight;

        // Shear XZ contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_s);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dVz_dgamma_xz * jacobian_weight;
        noalias(rRHS) -= global_size_N * Vz * jacobian_weight;

        // External body forces in w
        GlobalSizeVectorTransversalZ(global_size_N, N_shape);
        noalias(rRHS) += global_size_N * local_body_forces[2] * jacobian_weight * area;
    }

    RotateAll(rLHS, rRHS, r_geometry);

    KRATOS_CATCH("LinearTimoshenkoBeamElement3D2N::CalculateLocalSystem")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY

    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType local_trans_deflection_size = mat_size - 4 * number_of_nodes;

    if (rRHS.size() != mat_size) {
        rRHS.resize(mat_size, false);
    }
    rRHS.clear();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const double length = CalculateLength();
    const double Phi_rot_z  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length, 1);
    const double J    = 0.5 * length;
    const double area = GetCrossArea();

    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

    VectorType global_size_N(mat_size);
    VectorType N_u_derivatives(number_of_nodes);
    VectorType N_theta_derivatives(local_trans_deflection_size);
    VectorType N_theta(local_trans_deflection_size);
    VectorType N_derivatives(local_trans_deflection_size);
    VectorType N_u(number_of_nodes);
    VectorType N_shape(local_trans_deflection_size);
    VectorType N_s(local_trans_deflection_size);

    // Loop over the integration points (IP)
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        global_size_N.clear();
        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J;

        CalculateGeneralizedStrainsVector(strain_vector, length, 0.0, xi, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N  = r_generalized_stresses[0];
        const double Mx = r_generalized_stresses[1];
        const double My = r_generalized_stresses[2];
        const double Mz = r_generalized_stresses[3];
        const double Vy = r_generalized_stresses[4];
        const double Vz = r_generalized_stresses[5];

        // Axial DoFs shape functions (u and theta_x)
        GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, length, 0.0, xi);
        GetNu0ShapeFunctionsValues(N_u, length, 0.0, xi);

        // Axial contributions
        GlobalSizeAxialVector(global_size_N, N_u_derivatives);
        noalias(rRHS) -= global_size_N * N * jacobian_weight;

        // External body forces in u
        GlobalSizeAxialVector(global_size_N, N_u);
        noalias(rRHS) += global_size_N * local_body_forces[0] * jacobian_weight * area;

        // Torsional contributions
        GlobalSizeVectorAxialRotation(global_size_N, N_u_derivatives);
        noalias(rRHS) -= global_size_N * Mx * jacobian_weight;

        // Transverse DoFs shape functions {v, theta_z}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_z, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_z, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_z, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_z, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in z contributions
        GlobalSizeVectorTransversalY(global_size_N, N_theta_derivatives);
        noalias(rRHS) -= global_size_N * Mz * jacobian_weight;

        // Shear XY contributions
        GlobalSizeVectorTransversalY(global_size_N, N_s);
        noalias(rRHS) -= global_size_N * Vy * jacobian_weight;

        // External body forces in v
        GlobalSizeVectorTransversalY(global_size_N, N_shape);
        noalias(rRHS) += global_size_N * local_body_forces[1] * jacobian_weight * area;

        // Transverse DoFs shape functions {w, theta_y}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_y, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_y, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_y, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_y, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in y contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_theta_derivatives);
        noalias(rRHS) -= global_size_N * My * jacobian_weight;

        // Shear XZ contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_s);
        noalias(rRHS) -= global_size_N * Vz * jacobian_weight;

        // External body forces in w
        GlobalSizeVectorTransversalZ(global_size_N, N_shape);
        noalias(rRHS) += global_size_N * local_body_forces[2] * jacobian_weight * area;
    }

    RotateRHS(rRHS, r_geometry);

    KRATOS_CATCH("LinearTimoshenkoBeamElement3D2N::CalculateRightHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType mat_size = GetDoFsPerNode() * GetGeometry().size();
    rOutput.resize(r_integration_points.size());

    if (rVariable == AXIAL_STRAIN     ||
        rVariable == SHEAR_STRAIN_Y   ||
        rVariable == SHEAR_STRAIN_Z   ||
        rVariable == BENDING_STRAIN_X ||
        rVariable == BENDING_STRAIN_Y ||
        rVariable == BENDING_STRAIN_Z )
    {
        IndexType component = 0;
        if (rVariable == AXIAL_STRAIN) {
            component = 0;
        } else if (rVariable == BENDING_STRAIN_X) {
            component = 1;
        } else if (rVariable == BENDING_STRAIN_Y) {
            component = 2;
        } else if (rVariable == BENDING_STRAIN_Z) {
            component = 3;
        } else if (rVariable == SHEAR_STRAIN_Y) {
            component = 4;
        } else if (rVariable == SHEAR_STRAIN_Z) {
            component = 5;
        }
        const double length = CalculateLength();
        VectorType strain_vector(strain_size);

        VectorType nodal_values(mat_size);
        GetNodalValuesVector(nodal_values);

        // Loop over the integration points (IP)
        const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi     = r_integration_points[IP].X();
            CalculateGeneralizedStrainsVector(strain_vector, length, 0.0, xi, nodal_values);
            rOutput[IP] = strain_vector[component];
        }

    } else if (rVariable == AXIAL_FORCE ||
        rVariable == SHEAR_FORCE_Y      ||
        rVariable == SHEAR_FORCE_Z      ||
        rVariable == BENDING_MOMENT_X   ||
        rVariable == BENDING_MOMENT_Y   ||
        rVariable == BENDING_MOMENT_Z )
    {
        std::vector<Vector> pk2_stress;
        BaseType::CalculateOnIntegrationPoints(PK2_STRESS_VECTOR, pk2_stress, rCurrentProcessInfo);

        IndexType component = 0;
        if (rVariable == AXIAL_FORCE) {
            component = 0;
        } else if (rVariable == BENDING_MOMENT_X) {
            component = 1;
        } else if (rVariable == BENDING_MOMENT_Y) {
            component = 2;
        } else if (rVariable == BENDING_MOMENT_Z) {
            component = 3;
        } else if (rVariable == SHEAR_FORCE_Y) {
            component = 4;
        } else if (rVariable == SHEAR_FORCE_Z) {
            component = 5;
        }

        // Loop over the integration points
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            rOutput[IP] = pk2_stress[IP][component];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType mat_size = GetDoFsPerNode() * number_of_nodes;
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType local_trans_deflection_size = mat_size - 4 * number_of_nodes;

    if (rLHS.size1() != mat_size || rLHS.size2() != mat_size) {
        rLHS.resize(mat_size, mat_size, false);
    }
    rLHS.clear();

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double length = CalculateLength();
    const double Phi_rot_z  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length);
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhi(r_props, length, 1);
    const double J     = 0.5 * length;

    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    VectorType nodal_values(mat_size);
    GetNodalValuesVector(nodal_values);

    VectorType global_size_N(mat_size);
    VectorType N_u_derivatives(number_of_nodes);
    VectorType N_theta_derivatives(local_trans_deflection_size);
    VectorType N_theta(local_trans_deflection_size);
    VectorType N_derivatives(local_trans_deflection_size);
    VectorType N_u(number_of_nodes);
    VectorType N_shape(local_trans_deflection_size);
    VectorType N_s(local_trans_deflection_size);

    // Loop over the integration points (IP)
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        global_size_N.clear();
        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J;

        CalculateGeneralizedStrainsVector(strain_vector, length, 0.0, xi, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);

        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl       = r_constitutive_matrix(0, 0);
        const double dMx_dkappa_x = r_constitutive_matrix(1, 1);
        const double dMy_dkappa_y = r_constitutive_matrix(2, 2);
        const double dMz_dkappa_z = r_constitutive_matrix(3, 3);
        const double dVy_dgamma_xy = r_constitutive_matrix(4, 4);
        const double dVz_dgamma_xz = r_constitutive_matrix(5, 5);

        // Axial DoFs shape functions (u and theta_x)
        GetFirstDerivativesNu0ShapeFunctionsValues(N_u_derivatives, length, 0.0, xi);
        GetNu0ShapeFunctionsValues(N_u, length, 0.0, xi);

        // Axial contributions
        GlobalSizeAxialVector(global_size_N, N_u_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dN_dEl * jacobian_weight;

        // Torsional contributions
        GlobalSizeVectorAxialRotation(global_size_N, N_u_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMx_dkappa_x * jacobian_weight;

        // Transverse DoFs shape functions {v, theta_z}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_z, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_z, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_z, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_z, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in z contributions
        GlobalSizeVectorTransversalY(global_size_N, N_theta_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMz_dkappa_z * jacobian_weight;

        // Shear XY contributions
        GlobalSizeVectorTransversalY(global_size_N, N_s);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dVy_dgamma_xy * jacobian_weight;

        // Transverse DoFs shape functions {w, theta_y}
        GetNThetaShapeFunctionsValues(N_theta, length, Phi_rot_y, xi);
        GetFirstDerivativesNThetaShapeFunctionsValues(N_theta_derivatives, length, Phi_rot_y, xi);
        GetShapeFunctionsValues(N_shape, length, Phi_rot_y, xi);
        GetFirstDerivativesShapeFunctionsValues(N_derivatives, length, Phi_rot_y, xi);
        noalias(N_s) = N_derivatives - N_theta;

        // Bending in z contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_theta_derivatives);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dMy_dkappa_y * jacobian_weight;

        // Shear XZ contributions
        GlobalSizeVectorTransversalZ(global_size_N, N_s);
        noalias(rLHS) += outer_prod(global_size_N, global_size_N) * dVz_dgamma_xz * jacobian_weight;
    }

    RotateLHS(rLHS, r_geometry);

    KRATOS_CATCH("LinearTimoshenkoBeamElement3D2N::CalculateLeftHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(r_integration_points.size());

    if (rVariable == LOCAL_AXIS_1 || rVariable == LOCAL_AXIS_2) {
        const auto& r_local_axis = this->Has(rVariable) ? this->GetValue(rVariable) : ZeroVector(3);
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            noalias(rOutput[IP]) = r_local_axis;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::save(
    Serializer& rSerializer
    ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::load(
    Serializer& rSerializer
    )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
