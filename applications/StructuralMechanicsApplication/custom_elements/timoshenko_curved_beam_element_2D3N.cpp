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
#include "timoshenko_curved_beam_element_2D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void LinearTimoshenkoCurvedBeamElement2D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
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

void LinearTimoshenkoCurvedBeamElement2D3N::InitializeMaterial()
{
    KRATOS_TRY

    if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        const auto& r_properties = GetProperties();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer LinearTimoshenkoCurvedBeamElement2D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoCurvedBeamElement2D3N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement2D3N>
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

void LinearTimoshenkoCurvedBeamElement2D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = this->GetGeometry()[0].GetDofPosition(ROTATION_Z);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos ).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta
    rElementalDofList.resize(dofs_per_node * number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * dofs_per_node;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(ROTATION_Z    );
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetShapeFunctionsValues(
    const double xi
    ) const
{
    array_3 N;
    N[0] = 0.5 * xi * (xi - 1.0);
    N[1] = 0.5 * xi * (xi + 1.0);
    N[2] = 1.0 - std::pow(xi, 2);
    return N;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetFirstDerivativesShapeFunctionsValues(
    const double xi,
    const double J
    ) const
{
    return GetLocalFirstDerivativesShapeFunctionsValues(xi) / J;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetLocalFirstDerivativesShapeFunctionsValues(
    const double xi
    ) const
{
    array_3 dN;
    dN[0] = xi - 0.5;
    dN[1] = xi + 0.5;
    dN[2] = -2.0 * xi;
    return dN;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetSecondDerivativesShapeFunctionsValues(
    const double xi,
    const double J
    ) const
{
    return GetLocalSecondDerivativesShapeFunctionsValues(xi) / std::pow(J, 2);
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetLocalSecondDerivativesShapeFunctionsValues(
    const double xi
    ) const
{
    array_3 d2N;
    d2N[0] = 1.0;
    d2N[1] = 1.0;
    d2N[2] = -2.0;
    return d2N;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetTangentandTransverseUnitVectors(
    const double xi,
    array_3& rt,
    array_3& rn
    ) const
{
    const auto& r_geom = GetGeometry();

    array_3 dN_dxi   = GetLocalFirstDerivativesShapeFunctionsValues(xi);
    array_3 d2N_dxi2 = GetLocalSecondDerivativesShapeFunctionsValues(xi);

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    double d2x_dxi2 = 0.0;
    double d2y_dxi2 = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[i];
        dy_dxi += r_coords_node[1] * dN_dxi[i];

        d2x_dxi2 += r_coords_node[0] * d2N_dxi2[i];
        d2y_dxi2 += r_coords_node[1] * d2N_dxi2[i];
    }

    array_3 x_prime, x_2prime, b, aux;
    x_prime.clear();
    x_2prime.clear();
    rt.clear();
    rn.clear();
    b.clear();

    x_prime[0] = dx_dxi;
    x_prime[1] = dy_dxi;

    x_2prime[0] = d2x_dxi2;
    x_2prime[1] = d2y_dxi2;

    noalias(rt) = x_prime / norm_2(x_prime);

    noalias(aux) = MathUtils<double>::CrossProduct(x_prime, x_2prime);
    const double norm = norm_2(MathUtils<double>::CrossProduct(x_prime, x_2prime));
    if (norm != 0.0) // if the beam is not curved
        noalias(b) = aux / norm;
    else
        b[2] = 1.0;

    noalias(rn) = MathUtils<double>::CrossProduct(rt, b);
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 2, 2> LinearTimoshenkoCurvedBeamElement2D3N::GetFrenetSerretMatrix(
    const double xi,
    const array_3& rt,
    const array_3& rn
    ) const
{
    BoundedMatrix<double, 2, 2> T;
    T.clear();

    T(0, 0) = rt[0];
    T(0, 1) = rt[1];

    T(1, 0) = rn[0];
    T(1, 1) = rn[1];
    return T;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement2D3N::GetJacobian(
    const double xi
    ) const
{
    const auto& r_geom = GetGeometry();
    const array_3 dN_dxi = GetLocalFirstDerivativesShapeFunctionsValues(xi);

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[i];
        dy_dxi += r_coords_node[1] * dN_dxi[i];
    }
    return std::sqrt(std::pow(dx_dxi, 2) + std::pow(dy_dxi, 2));
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetNodalValuesVector(
    GlobalSizeVector& rNodalValues
    ) const
{
    const auto &r_geom = GetGeometry();

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

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement2D3N::array_3 LinearTimoshenkoCurvedBeamElement2D3N::GetBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::GetShapeFunctionsValuesGlobalVectors(
    const array_3 &rShapeFunctions,
    GlobalSizeVector &rNshape,
    GlobalSizeVector &rNu,
    GlobalSizeVector &rNtheta
    ) const
{
    // deflection v
    rNshape.clear();
    rNshape[1] = rShapeFunctions[0];
    rNshape[4] = rShapeFunctions[1];
    rNshape[7] = rShapeFunctions[2];

    // axial u
    rNu.clear();
    rNu[0] = rShapeFunctions[0];
    rNu[3] = rShapeFunctions[1];
    rNu[6] = rShapeFunctions[2];

    // rotation
    rNtheta.clear();
    rNtheta[2] = rShapeFunctions[0];
    rNtheta[5] = rShapeFunctions[1];
    rNtheta[8] = rShapeFunctions[2];
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rLHS.size1() != SystemSize || rLHS.size2() != SystemSize) {
        rLHS.resize(SystemSize, SystemSize, false);
    }
    noalias(rLHS) = ZeroMatrix(SystemSize, SystemSize);

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double area = r_props[CROSS_AREA];

    // Let's initialize the constitutive law values
    VectorType strain_vector(StrainSize), stress_vector(StrainSize);
    MatrixType constitutive_matrix(StrainSize, StrainSize);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values, B_b, dNu, dN_theta, N_shape, Nu, N_theta, dN_shape;
    BoundedMatrix<double, 2, 2> C_gamma, frenet_serret;
    BoundedVector<double, 2> N_forces, Gamma; // axial ans shear forces, strains
    BoundedMatrix<double, 2, 9> B_s, aux_B_s;
    C_gamma.clear();

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        const array_3 shape_functions   = GetShapeFunctionsValues(xi);
        const array_3 d_shape_functions = GetFirstDerivativesShapeFunctionsValues(xi, J);

        // Get shape functions for deflection, axial displ and rotations
        GetShapeFunctionsValuesGlobalVectors(shape_functions, N_shape, Nu, N_theta);
        GetShapeFunctionsValuesGlobalVectors(d_shape_functions, dN_shape, dNu, dN_theta);

        array_3 t, n;
        GetTangentandTransverseUnitVectors(xi, t, n);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
        noalias(B_b) =  dN_theta;

        // we fill aux_B_s
        for (IndexType i = 0; i < SystemSize; ++i) {
            aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
            aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
        }
        noalias(B_s) = prod(frenet_serret, aux_B_s);

        noalias(Gamma) = prod(B_s, nodal_values);
        strain_vector[0] = Gamma[0]; // axial strain
        strain_vector[1] = inner_prod(B_b, nodal_values); // curvature
        strain_vector[2] = Gamma[1]; // shear strain

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N = r_generalized_stresses[0];
        const double M = r_generalized_stresses[1];
        const double V = r_generalized_stresses[2];

        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl    = r_constitutive_matrix(0, 0);
        const double dM_dkappa = r_constitutive_matrix(1, 1);
        const double dV_dgamma = r_constitutive_matrix(2, 2);

        N_forces[0] = N;
        N_forces[1] = V;

        C_gamma(0, 0) = dN_dEl;
        C_gamma(1, 1) = dV_dgamma;

        noalias(rRHS) -= jacobian_weight * (prod(trans(B_s), N_forces) + M * B_b);
        noalias(rLHS) += jacobian_weight * (prod(trans(B_s), Matrix(prod(C_gamma, B_s))) + dM_dkappa * outer_prod(B_b, B_b));

        auto body_forces = GetBodyForce(*this, r_integration_points, IP);
        noalias(rRHS) += Nu      * body_forces[0] * jacobian_weight * area;
        noalias(rRHS) += N_shape * body_forces[1] * jacobian_weight * area;

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rLHS.size1() != SystemSize || rLHS.size2() != SystemSize) {
        rLHS.resize(SystemSize, SystemSize, false);
    }
    noalias(rLHS) = ZeroMatrix(SystemSize, SystemSize);

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Let's initialize the constitutive law values
    VectorType strain_vector(StrainSize), stress_vector(StrainSize);
    MatrixType constitutive_matrix(StrainSize, StrainSize);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values, B_b, dNu, dN_theta, N_shape, Nu, N_theta, dN_shape;
    BoundedMatrix<double, 2, 2> C_gamma, frenet_serret;
    BoundedVector<double, 2> N_forces, Gamma; // axial ans shear forces, strains
    BoundedMatrix<double, 2, 9> B_s, aux_B_s;
    C_gamma.clear();

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        const array_3 shape_functions   = GetShapeFunctionsValues(xi);
        const array_3 d_shape_functions = GetFirstDerivativesShapeFunctionsValues(xi, J);

        GetShapeFunctionsValuesGlobalVectors(shape_functions, N_shape, Nu, N_theta);
        GetShapeFunctionsValuesGlobalVectors(d_shape_functions, dN_shape, dNu, dN_theta);

        array_3 t, n;
        GetTangentandTransverseUnitVectors(xi, t, n);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
        noalias(B_b) =  dN_theta;

        // we fill aux_B_s
        for (IndexType i = 0; i < SystemSize; ++i) {
            aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
            aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
        }
        noalias(B_s) = prod(frenet_serret, aux_B_s);

        noalias(Gamma) = prod(B_s, nodal_values);
        strain_vector[0] = Gamma[0]; // axial strain
        strain_vector[2] = Gamma[1]; // shear strain
        strain_vector[1] = inner_prod(B_b, nodal_values); // curvature

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);

        const MatrixType& r_constitutive_matrix = cl_values.GetConstitutiveMatrix();
        const double dN_dEl    = r_constitutive_matrix(0, 0);
        const double dM_dkappa = r_constitutive_matrix(1, 1);
        const double dV_dgamma = r_constitutive_matrix(2, 2);

        C_gamma(0, 0) = dN_dEl;
        C_gamma(1, 1) = dV_dgamma;

        noalias(rLHS) += jacobian_weight * (prod(trans(B_s), Matrix(prod(C_gamma, B_s))) + dM_dkappa * outer_prod(B_b, B_b));

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const double area = r_props[CROSS_AREA];

    // Let's initialize the constitutive law values
    VectorType strain_vector(StrainSize), stress_vector(StrainSize);
    MatrixType constitutive_matrix(StrainSize, StrainSize);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values, B_b, dNu, dN_theta, N_shape, Nu, N_theta, dN_shape;
    BoundedMatrix<double, 2, 2> frenet_serret;
    BoundedVector<double, 2> N_forces, Gamma; // axial ans shear forces, strains
    BoundedMatrix<double, 2, 9> B_s, aux_B_s;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        const array_3 shape_functions   = GetShapeFunctionsValues(xi);
        const array_3 d_shape_functions = GetFirstDerivativesShapeFunctionsValues(xi, J);

        GetShapeFunctionsValuesGlobalVectors(shape_functions, N_shape, Nu, N_theta);
        GetShapeFunctionsValuesGlobalVectors(d_shape_functions, dN_shape, dNu, dN_theta);

        array_3 t, n;
        GetTangentandTransverseUnitVectors(xi, t, n);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
        noalias(B_b) =  dN_theta;

        // we fill aux_B_s
        for (IndexType i = 0; i < SystemSize; ++i) {
            aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
            aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
        }
        noalias(B_s) = prod(frenet_serret, aux_B_s);

        noalias(Gamma) = prod(B_s, nodal_values);
        strain_vector[0] = Gamma[0]; // axial strain
        strain_vector[2] = Gamma[1]; // shear strain
        strain_vector[1] = inner_prod(B_b, nodal_values); // curvature

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = cl_values.GetStressVector();
        const double N = r_generalized_stresses[0];
        const double M = r_generalized_stresses[1];
        const double V = r_generalized_stresses[2];

        N_forces[0] = N;
        N_forces[1] = V;

        noalias(rRHS) -= jacobian_weight * (prod(trans(B_s), N_forces) + M * B_b);

        auto body_forces = GetBodyForce(*this, r_integration_points, IP);
        noalias(rRHS) += Nu      * body_forces[0] * jacobian_weight * area;
        noalias(rRHS) += N_shape * body_forces[1] * jacobian_weight * area;

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(r_integration_points.size());
    const auto &r_props = GetProperties();

    if (rVariable == AXIAL_FORCE || rVariable == BENDING_MOMENT || rVariable == SHEAR_FORCE) {
        const auto &r_geometry = GetGeometry();
        const int strain_component = (rVariable == AXIAL_STRAIN) ? 0 : (rVariable == SHEAR_STRAIN) ? 2 : 1;

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Let's initialize the cl values
        VectorType strain_vector(StrainSize), stress_vector(StrainSize);
        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);

        GlobalSizeVector nodal_values, B_b, dNu, dN_theta, N_shape, Nu, N_theta, dN_shape;
        BoundedMatrix<double, 2, 2> frenet_serret;
        BoundedVector<double, 2> N_forces, Gamma; // axial and shear forces, strains
        BoundedMatrix<double, 2, 9> B_s, aux_B_s;

        // Loop over the integration points
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi     = r_integration_points[IP].X();
            const double J      = GetJacobian(xi);

            GetNodalValuesVector(nodal_values);

            const array_3 shape_functions   = GetShapeFunctionsValues(xi);
            const array_3 d_shape_functions = GetFirstDerivativesShapeFunctionsValues(xi, J);

            GetShapeFunctionsValuesGlobalVectors(shape_functions, N_shape, Nu, N_theta);
            GetShapeFunctionsValuesGlobalVectors(d_shape_functions, dN_shape, dNu, dN_theta);

            array_3 t, n;
            GetTangentandTransverseUnitVectors(xi, t, n);
            noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
            noalias(B_b) =  dN_theta;

            // we fill aux_B_s
            for (IndexType i = 0; i < SystemSize; ++i) {
                aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
                aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
            }
            noalias(B_s) = prod(frenet_serret, aux_B_s);

            noalias(Gamma) = prod(B_s, nodal_values);
            strain_vector[0] = Gamma[0]; // axial strain
            strain_vector[2] = Gamma[1]; // shear strain
            strain_vector[1] = inner_prod(B_b, nodal_values); // curvature

            mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
            const Vector &r_generalized_stresses = cl_values.GetStressVector();

            rOutput[IP] = r_generalized_stresses[strain_component];
        }
    } else if (rVariable == AXIAL_STRAIN || rVariable == BENDING_STRAIN || rVariable == SHEAR_STRAIN) {
        const int strain_component = (rVariable == AXIAL_STRAIN) ? 0 : (rVariable == SHEAR_STRAIN) ? 2 : 1;

        // Let's initialize the cl values
        VectorType strain_vector(StrainSize);

        GlobalSizeVector nodal_values, B_b, dNu, dN_theta, N_shape, Nu, N_theta, dN_shape;
        BoundedMatrix<double, 2, 2> frenet_serret;
        BoundedVector<double, 2> N_forces, Gamma; // axial and shear forces, strains
        BoundedMatrix<double, 2, 9> B_s, aux_B_s;

        // Loop over the integration points
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi     = r_integration_points[IP].X();
            const double J      = GetJacobian(xi);

            GetNodalValuesVector(nodal_values);

            const array_3 shape_functions   = GetShapeFunctionsValues(xi);
            const array_3 d_shape_functions = GetFirstDerivativesShapeFunctionsValues(xi, J);

            GetShapeFunctionsValuesGlobalVectors(shape_functions, N_shape, Nu, N_theta);
            GetShapeFunctionsValuesGlobalVectors(d_shape_functions, dN_shape, dNu, dN_theta);

            array_3 t, n;
            GetTangentandTransverseUnitVectors(xi, t, n);
            noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n);
            noalias(B_b) =  dN_theta;

            // we fill aux_B_s
            for (IndexType i = 0; i < SystemSize; ++i) {
                aux_B_s(0, i) = dNu[i] + t[1] * N_theta[i];
                aux_B_s(1, i) = dN_shape[i] - t[0] * N_theta[i];
            }
            noalias(B_s) = prod(frenet_serret, aux_B_s);

            noalias(Gamma) = prod(B_s, nodal_values);
            strain_vector[0] = Gamma[0]; // axial strain
            strain_vector[2] = Gamma[1]; // shear strain
            strain_vector[1] = inner_prod(B_b, nodal_values); // curvature

            rOutput[IP] = strain_vector[strain_component];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType r_integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != r_integration_points_number) {
            rValues.resize(r_integration_points_number);
        }
        for (IndexType point_number = 0; point_number < r_integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

int LinearTimoshenkoCurvedBeamElement2D3N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement2D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
