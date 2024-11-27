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
#include "linear_truss_element_3D.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
Element::Pointer LinearTrussElement3D<TNNodes, TDimension>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTrussElement3D<TNNodes, TDimension>::Pointer p_new_elem = Kratos::make_intrusive<LinearTrussElement3D<TNNodes, TDimension>>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    IndexType local_index = 0;

    if (rResult.size() != SystemSize)
        rResult.resize(SystemSize, false);

    const IndexType xpos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    rElementalDofList.resize(SystemSize);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * DofsPerNode;
        rElementalDofList[index]     = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index + 1] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index + 2] = r_geom[i].pGetDof(DISPLACEMENT_Z);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetShapeFunctionsValues(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();
    array_1d<double, NNodes> base_N;
    noalias(base_N) = GetBaseShapeFunctions(xi);

    if constexpr (NNodes == 2) {
        rN[0] = base_N[0];
        rN[3] = base_N[1];
    } else { // 3N
        rN[0] = base_N[0];
        rN[3] = base_N[1];
        rN[6] = base_N[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetShapeFunctionsValuesY(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();
    array_1d<double, NNodes> base_N;
    noalias(base_N) = GetBaseShapeFunctions(xi);

    if constexpr (NNodes == 2) {
        rN[1] = base_N[0];
        rN[4] = base_N[1];
    } else { // 3N
        rN[1] = base_N[0];
        rN[4] = base_N[1];
        rN[7] = base_N[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetShapeFunctionsValuesZ(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();
    array_1d<double, NNodes> base_N;
    noalias(base_N) = GetBaseShapeFunctions(xi);

    if constexpr (NNodes == 2) {
        rN[2] = base_N[0];
        rN[5] = base_N[1];
    } else { // 3N
        rN[2] = base_N[0];
        rN[5] = base_N[1];
        rN[8] = base_N[2];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetFirstDerivativesShapeFunctionsValues(
    SystemSizeBoundedArrayType& rdN_dX,
    const double Length,
    const double xi
    ) const
{
    if (rdN_dX.size() != SystemSize)
        rdN_dX.resize(SystemSize, false);

    rdN_dX.clear();

    if constexpr (NNodes == 2) {
        const double inverse_l = 1.0 / Length;
        rdN_dX[0] = -inverse_l;
        rdN_dX[3] = inverse_l;
    } else { // 3N
        rdN_dX[0] = xi - 0.5;
        rdN_dX[3] = xi + 0.5;
        rdN_dX[6] = -2.0 * xi;
        rdN_dX *= 2.0 / Length; // The Jacobian
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
BoundedMatrix<double, 3, 3> LinearTrussElement3D<TNNodes, TDimension>::GetFrenetSerretMatrix() const
{
    const auto &r_geom = GetGeometry();
    BoundedMatrix<double, 3, 3> T;
    T.clear(); // global to local

    array_1d<double, 3> t;
    array_1d<double, 3> n;
    array_1d<double, 3> m;

    // t is the axis of the truss
    noalias(t) = r_geom[1].GetInitialPosition() - r_geom[0].GetInitialPosition();
    t /= norm_2(t);

    n.clear();
    n[1] = 1.0;

    if (norm_2(t-n) <= 1.0e-8) { // colineal, hence we use another aux vector
        n.clear();
        n[2] = 1.0;
    }

    // Gram-Schmidt ortogonalization
    n = n - inner_prod(t, n) / inner_prod(t, t) * t;
    n /= norm_2(n);

    noalias(m) = MathUtils<double>::CrossProduct(t, n);

    T(0, 0) = t[0];
    T(0, 1) = t[1];
    T(0, 2) = t[2];

    T(1, 0) = n[0];
    T(1, 1) = n[1];
    T(1, 2) = n[2];

    T(2, 0) = m[0];
    T(2, 1) = m[1];
    T(2, 2) = m[2];

    return T;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::GetNodalValuesVector(SystemSizeBoundedArrayType& rNodalValues) const
{
    if (rNodalValues.size() != SystemSize)
        rNodalValues.resize(SystemSize, false);
    const auto &r_geom = GetGeometry();

    BoundedVector<double, SystemSize> global_values;

    // We fill the vector with global values
    for (SizeType i = 0; i < NNodes; ++i) {
        const auto& r_displ = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT);
        global_values[i * DofsPerNode]     = r_displ[0];
        global_values[i * DofsPerNode + 1] = r_displ[1];
        global_values[i * DofsPerNode + 2] = r_displ[2];
    }

        BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
        noalias(T) = GetFrenetSerretMatrix(); // global to local
        BoundedMatrix<double, SystemSize, SystemSize> global_size_T;

        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D3NTruss(T, global_size_T);
        }

        noalias(rNodalValues) = prod(global_size_T, global_values);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
array_1d<double, 3> LinearTrussElement3D<TNNodes, TDimension>::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    // const double angle = GetAngle();
    // const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);

    // const double c = std::cos(angle);
    // const double s = std::sin(angle);
    array_1d<double, 3> local_body_force = ZeroVector(3);
    // local_body_force[0] = c * body_force[0] + s * body_force[1];
    // local_body_force[1] = -s * body_force[0] + c * body_force[1];
    return local_body_force;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;

    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rLHS.size1() != SystemSize || rLHS.size2() != SystemSize) {
        rLHS.resize(SystemSize, SystemSize, false);
    }
    noalias(rLHS) = ZeroMatrix(SystemSize, SystemSize);

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double length = CalculateLength();
    const double J      = 0.5 * length;
    const double area   = r_props[CROSS_AREA];

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B, N_shape, N_shapeY;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, integration_points, IP);

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double jacobian_weight = weight * J * area;
        GetShapeFunctionsValues(N_shape, length, xi);
        GetShapeFunctionsValuesY(N_shapeY, length, xi);
        GetFirstDerivativesShapeFunctionsValues(B, length, xi);


        strain_vector[0] = inner_prod(B, nodal_values);
        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rLHS) += outer_prod(B, B) * constitutive_matrix(0, 0) * jacobian_weight;
        noalias(rRHS) -= B * stress_vector[0] * jacobian_weight;

        noalias(rRHS) += N_shape * local_body_forces[0] * jacobian_weight;
        noalias(rRHS) += N_shapeY * local_body_forces[1] * jacobian_weight;
    }
    RotateAll(rLHS, rRHS); // rotate to global

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
   KRATOS_TRY;
    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rLHS.size1() != SystemSize || rLHS.size2() != SystemSize) {
        rLHS.resize(SystemSize, SystemSize, false);
    }
    noalias(rLHS) = ZeroMatrix(SystemSize, SystemSize);

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double length = CalculateLength();
    const double J      = 0.5 * length;
    const double area   = r_props[CROSS_AREA];

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, integration_points, IP);

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double jacobian_weight = weight * J * area;
        GetFirstDerivativesShapeFunctionsValues(B, length, xi);

        strain_vector[0] = inner_prod(B, nodal_values);
        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rLHS) += outer_prod(B, B) * constitutive_matrix(0, 0) * jacobian_weight;
    }
    RotateLHS(rLHS); // rotate to global

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
   KRATOS_TRY;
    const auto &r_props = GetProperties();
    const auto &r_geometry = GetGeometry();

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    noalias(rRHS) = ZeroVector(SystemSize);

    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double length = CalculateLength();
    const double J      = 0.5 * length;
    const double area   = r_props[CROSS_AREA];

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B, N_shape, N_shapeY;

    // Loop over the integration points
    for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, integration_points, IP);

        const double xi     = integration_points[IP].X();
        const double weight = integration_points[IP].Weight();
        const double jacobian_weight = weight * J * area;
        GetShapeFunctionsValues(N_shape, length, xi);
        GetShapeFunctionsValuesY(N_shapeY, length, xi);
        GetFirstDerivativesShapeFunctionsValues(B, length, xi);

        strain_vector[0] = inner_prod(B, nodal_values);
        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rRHS) -= B * stress_vector[0] * jacobian_weight;

        noalias(rRHS) += N_shape * local_body_forces[0] * jacobian_weight;
        noalias(rRHS) += N_shapeY * local_body_forces[1] * jacobian_weight;
    }
    RotateRHS(rRHS); // rotate to global

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::RotateLHS(
    MatrixType& rLHS
)
{
    // const double angle = GetAngle();

    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, DofsPerNode, DofsPerNode> T, Tt;
    //     BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
    //     if constexpr (NNodes == 2) {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
    //     } else {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
    //     }

    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::RotateRHS(
    VectorType& rRHS
)
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
    //     BoundedMatrix<double, SystemSize, SystemSize> global_size_T;
    //     BoundedVector<double, SystemSize> local_rhs;
    //     noalias(local_rhs) = rRHS;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
    //     if constexpr (NNodes == 2) {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
    //     } else {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
    //     }

    //     noalias(rRHS) = prod(global_size_T, local_rhs);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS
)
{
    // const double angle = GetAngle();
    // if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
    //     BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
    //     BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
    //     BoundedVector<double, SystemSize> local_rhs;
    //     StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);

    //     if constexpr (NNodes == 2) {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
    //     } else {
    //         StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
    //     }

    //     noalias(local_rhs) = rRHS;
    //     noalias(rRHS) = prod(global_size_T, local_rhs);

    //     noalias(aux_product) = prod(rLHS, trans(global_size_T));
    //     noalias(rLHS) = prod(global_size_T, aux_product);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes, SizeType TDimension>
void LinearTrussElement3D<TNNodes, TDimension>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

template class LinearTrussElement3D<2, 3>;
template class LinearTrussElement3D<3, 3>;

} // Namespace Kratos
