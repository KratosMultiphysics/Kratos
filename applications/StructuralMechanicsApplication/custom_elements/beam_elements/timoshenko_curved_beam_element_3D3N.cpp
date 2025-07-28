// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "timoshenko_curved_beam_element_3D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void LinearTimoshenkoCurvedBeamElement3D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
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

void LinearTimoshenkoCurvedBeamElement3D3N::InitializeMaterial()
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

Element::Pointer LinearTimoshenkoCurvedBeamElement3D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTimoshenkoCurvedBeamElement3D3N::Pointer p_new_elem = Kratos::make_intrusive<LinearTimoshenkoCurvedBeamElement3D3N>
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

void LinearTimoshenkoCurvedBeamElement3D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    IndexType local_index = 0;

    if (rResult.size() != SystemSize)
        rResult.resize(SystemSize, false);

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos       ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1   ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2   ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X    , rot_pos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y    , rot_pos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z    , rot_pos + 2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = DoFperNode; // u, v, theta
    rElementalDofList.resize(SystemSize);

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

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetShapeFunctionsValues(
    const double xi
    ) const
{
    array_3 N;
    const double half_xi = 0.5 * xi;
    N[0] = half_xi * (xi - 1.0);
    N[1] = half_xi * (xi + 1.0);
    N[2] = 1.0 - xi * xi;
    return N;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetFirstDerivativesShapeFunctionsValues(
    const double xi,
    const double J
    ) const
{
    return GetLocalFirstDerivativesShapeFunctionsValues(xi) / J;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetLocalFirstDerivativesShapeFunctionsValues(
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

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetSecondDerivativesShapeFunctionsValues(
    const double xi,
    const double J
    ) const
{
    return GetLocalSecondDerivativesShapeFunctionsValues(xi) / std::pow(J, 2);
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetLocalSecondDerivativesShapeFunctionsValues(
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

void LinearTimoshenkoCurvedBeamElement3D3N::GetTangentandTransverseUnitVectors(
    const double xi,
    array_3& rt,
    array_3& rn,
    array_3& rb
    ) const
{
    const auto& r_geom = GetGeometry();

    array_3 dN_dxi   = GetLocalFirstDerivativesShapeFunctionsValues(xi);
    array_3 d2N_dxi2 = GetLocalSecondDerivativesShapeFunctionsValues(xi);

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;
    double dz_dxi = 0.0;

    double d2x_dxi2 = 0.0;
    double d2y_dxi2 = 0.0;
    double d2z_dxi2 = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[i];
        dy_dxi += r_coords_node[1] * dN_dxi[i];
        dz_dxi += r_coords_node[2] * dN_dxi[i];

        d2x_dxi2 += r_coords_node[0] * d2N_dxi2[i];
        d2y_dxi2 += r_coords_node[1] * d2N_dxi2[i];
        d2z_dxi2 += r_coords_node[2] * d2N_dxi2[i];
    }

    array_3 x_prime, x_2prime, b, aux;
    x_prime.clear();
    x_2prime.clear();
    rt.clear();
    rn.clear();
    rb.clear();

    x_prime[0] = dx_dxi;
    x_prime[1] = dy_dxi;
    x_prime[2] = dz_dxi;

    x_2prime[0] = d2x_dxi2;
    x_2prime[1] = d2y_dxi2;
    x_2prime[2] = d2z_dxi2;

    noalias(rt) = x_prime / norm_2(x_prime);

    noalias(aux) = MathUtils<double>::CrossProduct(x_prime, x_2prime);
    const double norm = norm_2(aux);

    // Compute the binormal vector
    if (norm > Tolerance) { // if the beam is curved
        noalias(rb) = aux / norm;
    } else {
        rb[2] = 1.0;
        // Ortogonalize b
        rb -= MathUtils<double>::Dot(rb, rt) * rt;
        double norm_b = norm_2(rb);

        if (norm_b > Tolerance) {
            rb /= norm_b;
        } else {
            rb[0] = 0.0;
            rb[1] = 1.0;
            rb[2] = 0.0;
            rb -= MathUtils<double>::Dot(rb, rt) * rt;
            rb /= norm_2(rb);
        }
    }

    noalias(rn) = MathUtils<double>::CrossProduct(rt, rb);
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 3, 3> LinearTimoshenkoCurvedBeamElement3D3N::GetFrenetSerretMatrix(
    const double xi,
    const array_3& rt,
    const array_3& rn,
    const array_3& rb
    ) const
{
    BoundedMatrix<double, Dimension, Dimension> T;
    T.clear();

    for (IndexType i = 0; i < Dimension; ++i) {
        T(0, i) = rt[i];
        T(1, i) = rb[i];
        T(2, i) = rn[i];
    }

    return T;
}

/***********************************************************************************/
/***********************************************************************************/

double LinearTimoshenkoCurvedBeamElement3D3N::GetJacobian(
    const double xi
    ) const
{
    const auto& r_geom = GetGeometry();
    const array_3 dN_dxi = GetLocalFirstDerivativesShapeFunctionsValues(xi);

    double dx_dxi = 0.0;
    double dy_dxi = 0.0;
    double dz_dxi = 0.0;

    for (IndexType i = 0; i < NumberOfNodes; ++i) {
        const auto &r_coords_node = r_geom[i].GetInitialPosition();
        dx_dxi += r_coords_node[0] * dN_dxi[i];
        dy_dxi += r_coords_node[1] * dN_dxi[i];
        dz_dxi += r_coords_node[2] * dN_dxi[i];
    }
    return std::sqrt(std::pow(dx_dxi, 2) + std::pow(dy_dxi, 2) + std::pow(dz_dxi, 2));
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::GetNodalValuesVector(
    GlobalSizeVector& rNodalValues
    ) const
{
    const auto &r_geom = GetGeometry();

    for (IndexType i_node = 0; i_node < NumberOfNodes; ++i_node) {
        const auto& r_node = r_geom[i_node];
        const auto& r_displ    = r_node.FastGetSolutionStepValue(DISPLACEMENT);
        const auto& r_rotation = r_node.FastGetSolutionStepValue(ROTATION);

        rNodalValues[DoFperNode * i_node]     = r_displ[0];
        rNodalValues[DoFperNode * i_node + 1] = r_displ[1];
        rNodalValues[DoFperNode * i_node + 2] = r_displ[2];
        rNodalValues[DoFperNode * i_node + 3] = r_rotation[0];
        rNodalValues[DoFperNode * i_node + 4] = r_rotation[1];
        rNodalValues[DoFperNode * i_node + 5] = r_rotation[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

BoundedMatrix<double, 6, 18> LinearTimoshenkoCurvedBeamElement3D3N::CalculateB(
    const array_3 & rN,
    const array_3 & rdN,
    const array_3 & rt
) const
{
    BoundedMatrix<double, 6, 18> B;
    B.clear();

    for (IndexType i_node = 0; i_node < NumberOfNodes; ++i_node) { // loop over node blocks
        const IndexType local_col_0 = i_node * DoFperNode;

        // We fill the derivative terms (Romero's notation)
        B(0, local_col_0)     = rdN[i_node]; // axial
        B(1, local_col_0 + 1) = rdN[i_node]; // shear xy
        B(2, local_col_0 + 2) = rdN[i_node]; // shear xz
        B(3, local_col_0 + 3) = rdN[i_node]; // kappa x
        B(4, local_col_0 + 4) = rdN[i_node]; // kappa y
        B(5, local_col_0 + 5) = rdN[i_node]; // kappa z

        // We fill the remaining terms
        // axial
        B(0, local_col_0 + 4) = -rt[2] * rN[i_node];  
        B(0, local_col_0 + 5) =  rt[1] * rN[i_node];  

        // shear xy
        B(1, local_col_0 + 3) =  rt[2] * rN[i_node];
        B(1, local_col_0 + 5) = -rt[0] * rN[i_node];

        // shear xz
        B(2, local_col_0 + 3) = -rt[1] * rN[i_node];
        B(2, local_col_0 + 4) =  rt[0] * rN[i_node];
    }
    return B;
}

/***********************************************************************************/
/***********************************************************************************/

LinearTimoshenkoCurvedBeamElement3D3N::array_3 LinearTimoshenkoCurvedBeamElement3D3N::GetBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    return StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

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

    const double area = GetCrossArea();

    // Let's initialize the constitutive law values
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values;
    BoundedMatrix<double, 3, 3> frenet_serret;
    BoundedMatrix<double, 6, 6> element_frenet_serret;
    BoundedMatrix<double, 6, 18> B;
    GlobalSizeVector Nu, Nv, Nw;
    array_3 t, n, b, shape_functions, d_shape_functions;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        noalias(shape_functions)   = GetShapeFunctionsValues(xi);
        noalias(d_shape_functions) = GetFirstDerivativesShapeFunctionsValues(xi, J);

        GetTangentandTransverseUnitVectors(xi, t, n, b);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n, b);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(frenet_serret, element_frenet_serret);
        noalias(B) = CalculateB(shape_functions, d_shape_functions, t);
        B = prod(element_frenet_serret, B);

        noalias(strain_vector) = CalculateStrainVector(B, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = ConvertGeneralizedVectorComponents(cl_values.GetStressVector());
        const MatrixType& r_constitutive_matrix = ConvertConstitutiveMatrixComponents(cl_values.GetConstitutiveMatrix());

        noalias(rRHS) -= jacobian_weight * prod(trans(B), r_generalized_stresses);
        noalias(rLHS) += jacobian_weight * prod(trans(B), Matrix(prod(r_constitutive_matrix, B)));

        auto body_forces = GetBodyForce(*this, r_integration_points, IP);
        CalculateDisplacementInterpolationVectors(Nu, Nv, Nw, shape_functions);
        const double area_weight = area * weight * J;
        noalias(rRHS) += Nu * body_forces[0] * area_weight;
        noalias(rRHS) += Nv * body_forces[1] * area_weight;
        noalias(rRHS) += Nw * body_forces[2] * area_weight;

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

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
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values;
    BoundedMatrix<double, 3, 3> frenet_serret;
    BoundedMatrix<double, 6, 6> element_frenet_serret;
    BoundedMatrix<double, 6, 18> B;
    array_3 t, n, b, shape_functions, d_shape_functions;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        noalias(shape_functions)   = GetShapeFunctionsValues(xi);
        noalias(d_shape_functions) = GetFirstDerivativesShapeFunctionsValues(xi, J);

        GetTangentandTransverseUnitVectors(xi, t, n, b);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n, b);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(frenet_serret, element_frenet_serret);
        noalias(B) = CalculateB(shape_functions, d_shape_functions, t);
        B = prod(element_frenet_serret, B);

        noalias(strain_vector) = CalculateStrainVector(B, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const MatrixType& r_constitutive_matrix = ConvertConstitutiveMatrixComponents(cl_values.GetConstitutiveMatrix());

        noalias(rLHS) += jacobian_weight * (prod(trans(B), Matrix(prod(r_constitutive_matrix, B))));

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;
    const auto &r_props    = GetProperties();
    const auto &r_geometry = GetGeometry();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const double area = GetCrossArea();

    // Let's initialize the constitutive law values
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    // Initialize required matrices/vectors...
    GlobalSizeVector nodal_values, Nu, Nv, Nw;
    BoundedMatrix<double, 3, 3> frenet_serret;
    BoundedMatrix<double, 6, 6> element_frenet_serret;
    BoundedMatrix<double, 6, 18> B;
    array_3 t, n, b, shape_functions, d_shape_functions;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double J      = GetJacobian(xi);
        const double jacobian_weight = weight * J;

        GetNodalValuesVector(nodal_values);

        noalias(shape_functions)   = GetShapeFunctionsValues(xi);
        noalias(d_shape_functions) = GetFirstDerivativesShapeFunctionsValues(xi, J);

        GetTangentandTransverseUnitVectors(xi, t, n, b);
        noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n, b);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(frenet_serret, element_frenet_serret);
        noalias(B) = CalculateB(shape_functions, d_shape_functions, t);
        B = prod(element_frenet_serret, B);

        noalias(strain_vector) = CalculateStrainVector(B, nodal_values);

        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        const Vector &r_generalized_stresses = ConvertGeneralizedVectorComponents(cl_values.GetStressVector());

        noalias(rRHS) -= jacobian_weight * prod(trans(B), r_generalized_stresses);

        auto body_forces = GetBodyForce(*this, r_integration_points, IP);
        CalculateDisplacementInterpolationVectors(Nu, Nv, Nw, shape_functions);
        const double area_weight = area * weight * J;
        noalias(rRHS) += Nu * body_forces[0] * area_weight;
        noalias(rRHS) += Nv * body_forces[1] * area_weight;
        noalias(rRHS) += Nw * body_forces[2] * area_weight;

    } // IP loop
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
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

        VectorType strain_vector(StrainSize);
        GlobalSizeVector nodal_values;
        GetNodalValuesVector(nodal_values);
        BoundedMatrix<double, 6, 18> B;
        BoundedMatrix<double, 3, 3> frenet_serret;
        BoundedMatrix<double, 6, 6> element_frenet_serret;
        array_3 t, n, b, shape_functions, d_shape_functions;

        // Loop over the integration points (IP)
        const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi = r_integration_points[IP].X();
            const double J  = GetJacobian(xi);

            GetNodalValuesVector(nodal_values);

            noalias(shape_functions)   = GetShapeFunctionsValues(xi);
            noalias(d_shape_functions) = GetFirstDerivativesShapeFunctionsValues(xi, J);

            GetTangentandTransverseUnitVectors(xi, t, n, b);
            noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n, b);
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(frenet_serret, element_frenet_serret);
            noalias(B) = CalculateB(shape_functions, d_shape_functions, t);
            B = prod(element_frenet_serret, B);
            noalias(strain_vector) = CalculateStrainVector(B, nodal_values);
            rOutput[IP] = strain_vector[component];
        }

    } else if (rVariable == AXIAL_FORCE ||
        rVariable == SHEAR_FORCE_Y      ||
        rVariable == SHEAR_FORCE_Z      ||
        rVariable == BENDING_MOMENT_X   ||
        rVariable == BENDING_MOMENT_Y   ||
        rVariable == BENDING_MOMENT_Z )
    {

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

        ConstitutiveLaw::Parameters cl_values(GetGeometry(), GetProperties(), rProcessInfo);
        VectorType strain_vector(StrainSize), stress_vector(StrainSize);
        StructuralMechanicsElementUtilities::InitializeConstitutiveLawValuesForStressCalculation(cl_values, strain_vector, stress_vector);
        GlobalSizeVector nodal_values;
        GetNodalValuesVector(nodal_values);
        BoundedMatrix<double, 6, 18> B;
        BoundedMatrix<double, 3, 3> frenet_serret;
        BoundedMatrix<double, 6, 6> element_frenet_serret;
        array_3 t, n, b, shape_functions, d_shape_functions;

        // Loop over the integration points (IP)
        const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());
        for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
            const double xi = r_integration_points[IP].X();
            const double J  = GetJacobian(xi);

            GetNodalValuesVector(nodal_values);

            noalias(shape_functions)   = GetShapeFunctionsValues(xi);
            noalias(d_shape_functions) = GetFirstDerivativesShapeFunctionsValues(xi, J);

            GetTangentandTransverseUnitVectors(xi, t, n, b);
            noalias(frenet_serret) = GetFrenetSerretMatrix(xi, t, n, b);
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NBeam(frenet_serret, element_frenet_serret);
            noalias(B) = CalculateB(shape_functions, d_shape_functions, t);
            B = prod(element_frenet_serret, B);
            noalias(strain_vector) = CalculateStrainVector(B, nodal_values);

            mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
            rOutput[IP] = cl_values.GetStressVector()[component];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Vector LinearTimoshenkoCurvedBeamElement3D3N::CalculateStrainVector(
    const BoundedMatrix<double, 6, 18>& rB,
    const GlobalSizeVector& rNodalValues
    ) const
{
    const Vector aux_strain = prod(rB, rNodalValues);

    Vector strain_vector(StrainSize);

    // We reorder from Romero et al. to the one in CL
    strain_vector[0] = aux_strain[0];
    strain_vector[1] = aux_strain[3];
    strain_vector[2] = aux_strain[4];
    strain_vector[3] = aux_strain[5];
    strain_vector[4] = aux_strain[1];
    strain_vector[5] = aux_strain[2];

    return strain_vector;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::CalculateOnIntegrationPoints(
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

int LinearTimoshenkoCurvedBeamElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void LinearTimoshenkoCurvedBeamElement3D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/


double LinearTimoshenkoCurvedBeamElement3D3N::GetCrossArea()
{
    return GetProperties()[CROSS_AREA];
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
