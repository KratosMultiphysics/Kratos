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

void LinearTimoshenkoBeamElement3D2N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = GetDoFsPerNode();

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X    , rot_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y    , rot_pos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos + 3).EquationId();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode();
    rElementalDofList.resize(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * dofs_per_node;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index + 3] = r_geom[i].pGetDof(ROTATION_X    );
        rElementalDofList[index + 4] = r_geom[i].pGetDof(ROTATION_Y    );
        rElementalDofList[index + 5] = r_geom[i].pGetDof(ROTATION_Z    );
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> LinearTimoshenkoBeamElement3D2N::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return array_1d<double, 3>();
    // const double angle = GetAngle();
    // const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);

    // const double c = std::cos(angle);
    // const double s = std::sin(angle);
    // array_1d<double, 3> local_body_force = ZeroVector(3);
    // local_body_force[0] = c * body_force[0] + s * body_force[1];
    // local_body_force[1] = -s * body_force[0] + c * body_force[1];
    // return local_body_force;
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
    
    BoundedVector<double, 12> global_values;
    BoundedMatrix<double, 3, 3> T;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(r_geom);

    const auto& r_displ_0    = r_geom[0].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_0 = r_geom[0].FastGetSolutionStepValue(ROTATION);

    // Here we rotate the vectors to local axes
    const VectorType& r_local_displ_0 = prod(T, r_displ_0);
    const VectorType& r_local_rot_0   = prod(T, r_rotation_0);

    global_values[0] = r_local_displ_0[0];
    global_values[1] = r_local_displ_0[1];
    global_values[2] = r_local_displ_0[2];

    global_values[3] = r_local_rot_0[0];
    global_values[4] = r_local_rot_0[1];
    global_values[5] = r_local_rot_0[2];

    const auto& r_displ_1    = r_geom[1].FastGetSolutionStepValue(DISPLACEMENT);
    const auto& r_rotation_1 = r_geom[1].FastGetSolutionStepValue(ROTATION);

    const VectorType& r_local_displ_1 = prod(T, r_displ_1);
    const VectorType& r_local_rot_1   = prod(T, r_rotation_1);

    global_values[6] = r_rotation_1[0];
    global_values[7] = r_rotation_1[1];
    global_values[8] = r_rotation_1[2];

    global_values[9]  = r_local_rot_1[0];
    global_values[10] = r_local_rot_1[1];
    global_values[11] = r_local_rot_1[2];

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
    const VectorType N_s = N_derivatives + N_theta;
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
    for (IndexType block = 0; block < 4; ++block) {
        for (IndexType i = 0; i < rT.size1(); ++i) {
            for (IndexType j = 0; j < rT.size2(); ++j) {
                rGlobalT(3 * block + i, 3 * block + j) = rT(i, j);
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
    const SizeType num_nodes = rGeometry.size();
    const SizeType global_size = GetDoFsPerNode() * num_nodes;
    
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T, aux_product;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(rGeometry);
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
    const SizeType num_nodes = rGeometry.size();
    const SizeType global_size = GetDoFsPerNode() * num_nodes;
    
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T;
    BoundedVector<double, 12> local_rhs;
    noalias(local_rhs) = rRHS;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(rGeometry);
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
    const SizeType num_nodes = rGeometry.size();
    const SizeType global_size = GetDoFsPerNode() * num_nodes;
    
    BoundedMatrix<double, 3, 3> T;
    BoundedMatrix<double, 12, 12> global_size_T, aux_product;
    noalias(T) = StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(rGeometry);
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
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhiY(r_props, Length);

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
    const double Phi_rot_y  = StructuralMechanicsElementUtilities::CalculatePhiY(r_props, length);
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
        const double dVz_dgamma_xz = r_constitutive_matrix(5, 6);
        
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













    }

    KRATOS_CATCH("LinearTimoshenkoBeamElement3D2N::CalculateLocalSystem")
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{

}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoBeamElement3D2N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{

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
