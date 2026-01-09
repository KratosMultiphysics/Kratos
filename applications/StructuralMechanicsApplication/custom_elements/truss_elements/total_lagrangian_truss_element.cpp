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
#include "total_lagrangian_truss_element.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (this->UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = CalculateIntegrationMethod();
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

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::InitializeMaterial()
{
    KRATOS_TRY
    const auto &r_props = GetProperties();

    if (r_props[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry   = GetGeometry();
        auto N_values            = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = r_props[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_props, r_geometry, N_values);
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
Element::Pointer TotalLagrangianTrussElement<TDimension>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    TotalLagrangianTrussElement<TDimension>::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianTrussElement<TDimension>>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::EquationIdVector(
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
        if constexpr (Dimension == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetDofList(
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
        if constexpr (Dimension == 3) {
            rElementalDofList[index + 2] = r_geom[i].pGetDof(DISPLACEMENT_Z);
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetShapeFunctionsValues(
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

    if constexpr (Dimension == 2) {
        rN[0] = base_N[0];
        rN[2] = base_N[1];
    } else {
        rN[0] = base_N[0];
        rN[3] = base_N[1];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetShapeFunctionsValuesY(
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

    if constexpr (Dimension == 2) {
        rN[1] = base_N[0];
        rN[3] = base_N[1];
    } else {
        rN[1] = base_N[0];
        rN[4] = base_N[1];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetShapeFunctionsValuesZ(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();

    if constexpr (Dimension == 3) {
        array_1d<double, NNodes> base_N;
        noalias(base_N) = GetBaseShapeFunctions(xi);
        rN[2] = base_N[0];
        rN[5] = base_N[1];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetFirstDerivativesShapeFunctionsValues(
    SystemSizeBoundedArrayType& rdN_dX,
    const double RefLength,
    const double xi
    ) const
{
    if (rdN_dX.size() != SystemSize)
        rdN_dX.resize(SystemSize, false);

    rdN_dX.clear();

    Vector coord(3);
    Matrix dN_de(NNodes, 1);
    coord.clear();
    coord[0] = xi;
    GetGeometry().ShapeFunctionsLocalGradients(dN_de, coord);

    if constexpr (Dimension == 2) {
        rdN_dX[0] = dN_de(0, 0);
        rdN_dX[2] = dN_de(1, 0);
    } else {
        rdN_dX[0] = dN_de(0, 0);
        rdN_dX[3] = dN_de(1, 0);
    }
    rdN_dX *= 2.0 / RefLength; // The Jacobian
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
BoundedMatrix<double, 3, 3> TotalLagrangianTrussElement<TDimension>::GetFrenetSerretMatrix() const
{
    // The Frenet-Serret in the CURRENT configuration
    return StructuralMechanicsElementUtilities::GetFrenetSerretMatrix3D(GetGeometry(), true);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetNodalValuesVector(
    SystemSizeBoundedArrayType& rNodalValues
) const
{
    if (rNodalValues.size() != SystemSize)
        rNodalValues.resize(SystemSize, false);

    const auto &r_geom = GetGeometry();
    BoundedVector<double, SystemSize> global_values;
    BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
    BoundedMatrix<double, SystemSize, SystemSize> global_size_T;

    if constexpr (Dimension == 2) {
        const double angle = GetCurrentAngle();
        // We fill the vector with global values
        for (SizeType i = 0; i < NNodes; ++i) {
            const auto& r_displ = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            global_values[i * DofsPerNode]     = r_displ[0];
            global_values[i * DofsPerNode + 1] = r_displ[1];
        }

        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
        noalias(rNodalValues) = prod(trans(global_size_T), global_values);
    } else {
        // We fill the vector with global values
        for (SizeType i = 0; i < NNodes; ++i) {
            const auto& r_displ = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT);
            global_values[i * DofsPerNode]     = r_displ[0];
            global_values[i * DofsPerNode + 1] = r_displ[1];
            global_values[i * DofsPerNode + 2] = r_displ[2];
        }

        noalias(T) = GetFrenetSerretMatrix(); // global to local
        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D3NTruss(T, global_size_T);
        }
        noalias(rNodalValues) = prod(global_size_T, global_values);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
array_1d<double, 3> TotalLagrangianTrussElement<TDimension>::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    array_1d<double, 3> local_body_force;

    if constexpr (Dimension == 2) {
        const double angle = GetCurrentAngle();

        const double c = std::cos(angle);
        const double s = std::sin(angle);
        local_body_force[0] = c * body_force[0] + s * body_force[1];
        local_body_force[1] = -s * body_force[0] + c * body_force[1];
        local_body_force[2] = 0.0;
        return local_body_force;
    } else {
        BoundedMatrix<double, Dimension, Dimension> T;
        noalias(T) = GetFrenetSerretMatrix(); // global to local
        noalias(local_body_force) = prod(T, body_force);
        return local_body_force;
    }

    return array_1d<double, 3>();
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
typename TotalLagrangianTrussElement<TDimension>::MatrixType TotalLagrangianTrussElement<TDimension>::
CalculateGeometricStiffnessMatrix(
    const double Stress_x,
    const double Ref_Length
)
{
    KRATOS_TRY
    MatrixType K_geometric(SystemSize, SystemSize);
    K_geometric.clear();

    for (IndexType i = 0; i < SystemSize; ++i) {
        for (IndexType j = 0; j < SystemSize; ++j) {
            if (i == j) {
                K_geometric(i, i) = 1.0;
            } else if (i == j + Dimension) {
                K_geometric(i, j) = -1.0;
            } else if (i + Dimension == j) {
                K_geometric(i, j) = -1.0;
            }
        }
    }
    return K_geometric * Stress_x / std::pow(Ref_Length, 2);

    KRATOS_CATCH("CalculateGeometricStiffnessMatrix")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateLocalSystem(
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
    rLHS.clear();

    if (rRHS.size() != SystemSize) {
        rRHS.resize(SystemSize, false);
    }
    rRHS.clear();

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double ref_length = CalculateReferenceLength();
    const double curr_length = CalculateCurrentLength();
    const double J = 0.5 * ref_length; // Jacobian to the natural coordinates
    const double ref_area = r_props[CROSS_AREA];
    const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
    const auto pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B, N_shape, N_shapeY, N_shapeZ;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J * ref_area;
        GetShapeFunctionsValues(N_shape, ref_length, xi);
        GetShapeFunctionsValuesY(N_shapeY, ref_length, xi);
        GetShapeFunctionsValuesZ(N_shapeZ, ref_length, xi);
        GetFirstDerivativesShapeFunctionsValues(B, ref_length, xi);

        strain_vector[0] = CalculateGreenLagrangeStrain(curr_length, ref_length);

        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rRHS) -= Fxx * B * (stress_vector[0] + pre_stress) * jacobian_weight;

        noalias(rLHS) += Fxx * Fxx * outer_prod(B, B) * constitutive_matrix(0, 0) * jacobian_weight; // Material stiffness
        noalias(rLHS) += CalculateGeometricStiffnessMatrix(stress_vector[0] + pre_stress, ref_length) * jacobian_weight; // Geometric stiffness

        noalias(rRHS) += N_shape  * local_body_forces[0] * jacobian_weight;
        noalias(rRHS) += N_shapeY * local_body_forces[1] * jacobian_weight;
        noalias(rRHS) += N_shapeZ * local_body_forces[2] * jacobian_weight;
    }
    RotateAll(rLHS, rRHS); // rotate to global

    KRATOS_CATCH("CalculateLocalSystem")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateLeftHandSide(
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
    rLHS.clear();

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , false);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const double ref_length = CalculateReferenceLength();
    const double curr_length = CalculateCurrentLength();
    const double J = 0.5 * ref_length; // Jacobian to the natural coordinates
    const double ref_area = r_props[CROSS_AREA];
    const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
    const auto pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B; // Derivatives of the shape functions with respect to reference coordinates

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J * ref_area;
        GetFirstDerivativesShapeFunctionsValues(B, ref_length, xi);

        strain_vector[0] = CalculateGreenLagrangeStrain(curr_length, ref_length);

        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rLHS) += Fxx * Fxx * outer_prod(B, B) * constitutive_matrix(0, 0) * jacobian_weight; // Material stiffness
        noalias(rLHS) += CalculateGeometricStiffnessMatrix(stress_vector[0] + pre_stress, ref_length) * jacobian_weight; // Geometric stiffness

    }
    RotateLHS(rLHS); // rotate to global

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateRightHandSide(
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
    rRHS.clear();

    const auto& r_integration_points = IntegrationPoints(GetIntegrationMethod());

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const double ref_length = CalculateReferenceLength();
    const double curr_length = CalculateCurrentLength();
    const double J = 0.5 * ref_length;
    const double ref_area = r_props[CROSS_AREA];
    const double Fxx = curr_length / ref_length;
    const auto pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] : 0.0;

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);
    SystemSizeBoundedArrayType nodal_values(SystemSize);
    GetNodalValuesVector(nodal_values); // In local axes

    SystemSizeBoundedArrayType B, N_shape, N_shapeY, N_shapeZ;

    // Loop over the integration points
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        const auto local_body_forces = GetLocalAxesBodyForce(*this, r_integration_points, IP);

        const double xi     = r_integration_points[IP].X();
        const double weight = r_integration_points[IP].Weight();
        const double jacobian_weight = weight * J * ref_area;
        GetShapeFunctionsValues(N_shape, ref_length, xi);
        GetShapeFunctionsValuesY(N_shapeY, ref_length, xi);
        GetShapeFunctionsValuesZ(N_shapeZ, ref_length, xi);
        GetFirstDerivativesShapeFunctionsValues(B, ref_length, xi);

        strain_vector[0] = CalculateGreenLagrangeStrain(curr_length, ref_length);

        mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix

        noalias(rRHS) -= Fxx * B * (stress_vector[0] + pre_stress) * jacobian_weight;

        noalias(rRHS) += N_shape  * local_body_forces[0] * jacobian_weight;
        noalias(rRHS) += N_shapeY * local_body_forces[1] * jacobian_weight;
        noalias(rRHS) += N_shapeZ * local_body_forces[2] * jacobian_weight;
    }
    RotateRHS(rRHS); // rotate to global

    KRATOS_CATCH("CalculateRightHandSide")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::RotateLHS(
    MatrixType& rLHS
)
{
    BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
    BoundedMatrix<double, DofsPerNode, DofsPerNode> T;

    if constexpr (Dimension == 2) {
        const double angle = GetCurrentAngle();

        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);

    } else {
        noalias(T) = trans(GetFrenetSerretMatrix()); // global to local
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D2NTruss(T, global_size_T);
    }

    noalias(aux_product) = prod(rLHS, trans(global_size_T));
    noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::RotateRHS(
    VectorType& rRHS
)
{
    BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
    BoundedMatrix<double, SystemSize, SystemSize> global_size_T;
    BoundedVector<double, SystemSize> local_rhs;

    if constexpr (Dimension == 2) {
        const double angle = GetCurrentAngle();
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);

    } else {
        noalias(T) = trans(GetFrenetSerretMatrix()); // global to local
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D2NTruss(T, global_size_T);
    }
    noalias(local_rhs) = rRHS;
    noalias(rRHS) = prod(global_size_T, local_rhs);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS
)
{
    BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
    BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
    BoundedVector<double, SystemSize> local_rhs;

    if constexpr (Dimension == 2) {
        const double angle = GetCurrentAngle();
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);

    } else {
        noalias(T) = trans(GetFrenetSerretMatrix()); // global to local
        StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor3D2NTruss(T, global_size_T);
    }

    noalias(local_rhs) = rRHS;
    noalias(rRHS) = prod(global_size_T, local_rhs);

    noalias(aux_product) = prod(rLHS, trans(global_size_T));
    noalias(rLHS) = prod(global_size_T, aux_product);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateMassMatrix(
    TotalLagrangianTrussElement<TDimension>::MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY

    const auto &r_props = GetProperties();

    if (rMassMatrix.size1() != SystemSize || rMassMatrix.size2() != SystemSize) {
        rMassMatrix.resize(SystemSize, SystemSize, false);
    }
    rMassMatrix.clear();

    const double A = r_props[CROSS_AREA];
    const double L0 = CalculateReferenceLength();
    const double rho = r_props[DENSITY];
    const double total_mass = rho * A * L0;

    if (StructuralMechanicsElementUtilities::ComputeLumpedMassMatrix(r_props, rCurrentProcessInfo)) {
        // Compute lumped mass matrix
        noalias(rMassMatrix) = 0.5 * total_mass * IdentityMatrix(SystemSize, SystemSize);
    } else {
        // Compute consistent mass matrix
        for (IndexType i = 0; i < SystemSize; ++i) {
            for (IndexType j = 0; j < SystemSize; ++j) {
                if (i == j) {
                    rMassMatrix(i, i) = 2.0;
                } else if (i == j + Dimension) {
                    rMassMatrix(i, j) = 1.0;
                } else if (i + Dimension == j) {
                    rMassMatrix(i, j) = 1.0;
                }
            }
        }
        rMassMatrix *= total_mass / 6.0;
    }

    KRATOS_CATCH("CalculateMassMatrix")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateDampingMatrix(
    TotalLagrangianTrussElement<TDimension>::MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
)
{
    KRATOS_TRY
    StructuralMechanicsElementUtilities::CalculateRayleighDampingMatrix(
        *this,
        rDampingMatrix,
        rCurrentProcessInfo,
        SystemSize);
    KRATOS_CATCH("CalculateDampingMatrix")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(integration_points.size());

    if (rVariable == AXIAL_FORCE) {
        ConstitutiveLaw::Parameters cl_values(GetGeometry(), GetProperties(), rProcessInfo);
        VectorType strain_vector(1), stress_vector(1);
        MatrixType C(1,1);
        StructuralMechanicsElementUtilities::InitializeConstitutiveLawValuesForStressCalculation(cl_values, strain_vector, stress_vector, C);

        const double ref_length = CalculateReferenceLength();
        const double curr_length = CalculateCurrentLength();
        SystemSizeBoundedArrayType nodal_values(SystemSize);
        GetNodalValuesVector(nodal_values);

        SystemSizeBoundedArrayType B;
        const double ref_area = GetProperties()[CROSS_AREA];

        double pre_stress = 0.0;
        if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            pre_stress = GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
            GetFirstDerivativesShapeFunctionsValues(B, ref_length, xi);

            strain_vector[0] = CalculateGreenLagrangeStrain(curr_length, ref_length);

            mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values);

            rOutput[IP] = (cl_values.GetStressVector()[0] + pre_stress) * ref_area;
        }
    } else if (rVariable == AXIAL_STRAIN) {
        ConstitutiveLaw::Parameters cl_values(GetGeometry(), GetProperties(), rProcessInfo);
        VectorType strain_vector(1), stress_vector(1);
        MatrixType C(1,1);
        StructuralMechanicsElementUtilities::InitializeConstitutiveLawValuesForStressCalculation(cl_values, strain_vector, stress_vector, C);

        SystemSizeBoundedArrayType nodal_values(SystemSize);
        GetNodalValuesVector(nodal_values);

        const double ref_length = CalculateReferenceLength();
        const double curr_length = CalculateCurrentLength();

        SystemSizeBoundedArrayType B;

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
            GetFirstDerivativesShapeFunctionsValues(B, ref_length, xi);
            rOutput[IP] = CalculateGreenLagrangeStrain(curr_length, ref_length);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(integration_points.size());

    if (rVariable == PK2_STRESS_VECTOR) {
        ConstitutiveLaw::Parameters cl_values(GetGeometry(), GetProperties(), rProcessInfo);
        VectorType strain_vector(1), stress_vector(1);
        MatrixType constitutive_matrix(1,1);
        StructuralMechanicsElementUtilities::InitializeConstitutiveLawValuesForStressCalculation(cl_values, strain_vector, stress_vector, constitutive_matrix);

        const double ref_length = CalculateReferenceLength();
        const double curr_length = CalculateCurrentLength();

        SystemSizeBoundedArrayType nodal_values(SystemSize);
        GetNodalValuesVector(nodal_values);

        SystemSizeBoundedArrayType dN_dX;

        for (SizeType integration_point = 0; integration_point < integration_points.size(); ++integration_point) {
            GetFirstDerivativesShapeFunctionsValues(dN_dX, ref_length, integration_points[integration_point].X());
            strain_vector[0] = CalculateGreenLagrangeStrain(curr_length, ref_length);
            mConstitutiveLawVector[integration_point]->CalculateMaterialResponsePK2(cl_values);

            auto& r_stress = cl_values.GetStressVector()[0];
            if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
                r_stress += GetProperties()[TRUSS_PRESTRESS_PK2];
            }
            rOutput[integration_point] = ScalarVector(1, r_stress);
        }
    }
}
/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
int TotalLagrangianTrussElement<TDimension>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CROSS_AREA)) << "CROSS_AREA not defined in the properties" << std::endl;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
Vector TotalLagrangianTrussElement<TDimension>::GetBaseShapeFunctions(const double xi) const
{
    Vector coord(3), N(NNodes);
    coord.clear();
    coord[0] = xi;
    GetGeometry().ShapeFunctionsValues(N, coord);
    return N;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
typename TotalLagrangianTrussElement<TDimension>::MatrixType TotalLagrangianTrussElement<TDimension>::
CalculateClosedFormK(
    const double Sxx
)
{
    MatrixType K = ZeroMatrix(SystemSize, SystemSize);
    const auto &r_props = GetProperties();
    const double E = r_props[YOUNG_MODULUS];
    const double A = r_props[CROSS_AREA];
    const double L0 = CalculateReferenceLength();
    const double L = CalculateCurrentLength();
    const double Fxx = L / L0;
    const double theta = GetCurrentAngle();

    const double k = Fxx * Fxx * E * A / L0;

    // Material stiffness matrix
    K(0,0) =  k * std::cos(theta) * std::cos(theta);
    K(0,1) =  k * std::cos(theta) * std::sin(theta);
    K(1,0) =  K(0,1);
    K(1,1) =  k * std::sin(theta) * std::sin(theta);

    K(0,2) = -k * std::cos(theta) * std::cos(theta);
    K(1,2) = -k * std::sin(theta) * std::cos(theta);
    K(2,0) =  K(0,2);
    K(2,1) =  K(1,2);
    K(2,2) =  k * std::cos(theta) * std::cos(theta);

    K(3,3) =  k * std::sin(theta) * std::sin(theta);
    K(2,3) = k * std::cos(theta) * std::sin(theta);
    K(3,2) = K(2,3);
    K(1,3) = -k * std::sin(theta) * std::sin(theta);
    K(3,1) = K(1,3);
    K(0,3) = -k * std::cos(theta) * std::sin(theta);
    K(3,0) = K(0,3);

    // Geometric stiffness matrix
    const double k_fact_g = Sxx * A / L0;

    for (IndexType i = 0; i < SystemSize; ++i) {
        for (IndexType j = 0; j < SystemSize; ++j) {
            if (i == j) {
                K(i, i) += k_fact_g;
            } else if (i == j + Dimension) {
                K(i, j) += -k_fact_g;
            } else if (i + Dimension == j) {
                K(i, j) += -k_fact_g;
            }
        }
    }

    return K;
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetValuesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY
    if (rValues.size() != SystemSize) {
        rValues.resize(SystemSize, false);
    }

    for (IndexType i = 0; i < NNodes; ++i) {
        int index = i * Dimension;
        const auto& r_disp =
            GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);

        rValues[index] = r_disp[0];
        rValues[index + 1] = r_disp[1];
        if constexpr (Dimension == 3)
            rValues[index + 2] = r_disp[2];
    }
    KRATOS_CATCH("GetValuesVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY
    if (rValues.size() != SystemSize) {
        rValues.resize(SystemSize, false);
    }

    for (IndexType i = 0; i < NNodes; ++i) {
        int index = i * Dimension;
        const auto& r_vel =
            GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);

        rValues[index] = r_vel[0];
        rValues[index + 1] = r_vel[1];
        if constexpr (Dimension == 3)
            rValues[index + 2] = r_vel[2];
    }
    KRATOS_CATCH("GetFirstDerivativesVector")
}

/***********************************************************************************/
/***********************************************************************************/

template <SizeType TDimension>
void TotalLagrangianTrussElement<TDimension>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{

    KRATOS_TRY
    if (rValues.size() != SystemSize) {
        rValues.resize(SystemSize, false);
    }

    for (IndexType i = 0; i < NNodes; ++i) {
        int index = i * Dimension;
        const auto& r_acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);

        rValues[index] = r_acc[0];
        rValues[index + 1] = r_acc[1];
        if constexpr (Dimension == 3)
            rValues[index + 2] = r_acc[2];
    }

    KRATOS_CATCH("GetSecondDerivativesVector")
}

/***********************************************************************************/
/***********************************************************************************/

template class TotalLagrangianTrussElement<2>;
template class TotalLagrangianTrussElement<3>;

} // Namespace Kratos
