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
#include "linear_truss_element_2D.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
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

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::InitializeMaterial()
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

template<SizeType TNNodes>
Element::Pointer LinearTrussElement2D<TNNodes>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    LinearTrussElement2D<TNNodes>::Pointer p_new_elem = Kratos::make_intrusive<LinearTrussElement2D<TNNodes>>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::EquationIdVector(
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
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetDofList(
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
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetShapeFunctionsValues(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();
    if constexpr (NNodes == 2) {
        rN[0] = 0.5 * (1.0 - xi);
        rN[2] = 0.5 * (1.0 + xi);
    } else { // 3N
        rN[0] = 0.5 * xi * (xi - 1.0);
        rN[2] = 0.5 * xi * (xi + 1.0);
        rN[4] = (1.0 - std::pow(xi, 2));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetShapeFunctionsValuesY(
    SystemSizeBoundedArrayType& rN,
    const double Length,
    const double xi
    ) const
{
    if (rN.size() != SystemSize)
        rN.resize(SystemSize, false);

    rN.clear();
    if constexpr (NNodes == 2) {
        rN[1] = 0.5 * (1.0 - xi);
        rN[3] = 0.5 * (1.0 + xi);
    } else { // 3N
        rN[1] = 0.5 * xi * (xi - 1.0);
        rN[3] = 0.5 * xi * (xi + 1.0);
        rN[5] = (1.0 - std::pow(xi, 2));
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetFirstDerivativesShapeFunctionsValues(
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
        rdN_dX[2] = inverse_l;
    } else { // 3N
        rdN_dX[0] = xi - 0.5;
        rdN_dX[2] = xi + 0.5;
        rdN_dX[4] = -2.0 * xi;
        rdN_dX *= 2.0 / Length; // The Jacobian
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::GetNodalValuesVector(SystemSizeBoundedArrayType& rNodalValues) const
{
    if (rNodalValues.size() != SystemSize)
        rNodalValues.resize(SystemSize, false);
    const auto &r_geom = GetGeometry();

    const double angle = GetAngle();

    BoundedVector<double, SystemSize> global_values;

    // We fill the vector with global values
    for (SizeType i = 0; i < NNodes; ++i) {
        const auto& r_displ = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT);
        global_values[i * DofsPerNode]     = r_displ[0];
        global_values[i * DofsPerNode + 1] = r_displ[1];
    }

    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
        BoundedMatrix<double, SystemSize, SystemSize> global_size_T;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
        }
        // We rotate to local axes
        noalias(rNodalValues) = prod(trans(global_size_T), global_values);
    } else {
        noalias(rNodalValues) = global_values;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
array_1d<double, 3> LinearTrussElement2D<TNNodes>::GetLocalAxesBodyForce(
    const Element &rElement,
    const GeometryType::IntegrationPointsArrayType &rIntegrationPoints,
    const IndexType PointNumber
    ) const
{
    const double angle = GetAngle();
    const auto body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);

    const double c = std::cos(angle);
    const double s = std::sin(angle);
    array_1d<double, 3> local_body_force = ZeroVector(3);
    local_body_force[0] = c * body_force[0] + s * body_force[1];
    local_body_force[1] = -s * body_force[0] + c * body_force[1];
    return local_body_force;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateLocalSystem(
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

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateLeftHandSide(
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

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateRightHandSide(
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

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateLHS(
    MatrixType& rLHS
)
{
    const double angle = GetAngle();

    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, DofsPerNode, DofsPerNode> T, Tt;
        BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
        }

        noalias(aux_product) = prod(rLHS, trans(global_size_T));
        noalias(rLHS) = prod(global_size_T, aux_product);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateRHS(
    VectorType& rRHS
)
{
    const double angle = GetAngle();
    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
        BoundedMatrix<double, SystemSize, SystemSize> global_size_T;
        BoundedVector<double, SystemSize> local_rhs;
        noalias(local_rhs) = rRHS;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);
        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
        }

        noalias(rRHS) = prod(global_size_T, local_rhs);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::RotateAll(
    MatrixType& rLHS,
    VectorType& rRHS
)
{
    const double angle = GetAngle();
    if (std::abs(angle) > std::numeric_limits<double>::epsilon()) {
        BoundedMatrix<double, DofsPerNode, DofsPerNode> T;
        BoundedMatrix<double, SystemSize, SystemSize> global_size_T, aux_product;
        BoundedVector<double, SystemSize> local_rhs;
        StructuralMechanicsElementUtilities::BuildRotationMatrixForTruss(T, angle);

        if constexpr (NNodes == 2) {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D2NTruss(T, global_size_T);
        } else {
            StructuralMechanicsElementUtilities::BuildElementSizeRotationMatrixFor2D3NTruss(T, global_size_T);
        }

        noalias(local_rhs) = rRHS;
        noalias(rRHS) = prod(global_size_T, local_rhs);

        noalias(aux_product) = prod(rLHS, trans(global_size_T));
        noalias(rLHS) = prod(global_size_T, aux_product);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo
    )
{
    const auto& integration_points = IntegrationPoints(GetIntegrationMethod());
    rOutput.resize(integration_points.size());

    if (rVariable == AXIAL_FORCE) {
        const auto &r_props = GetProperties();
        const auto &r_geometry = GetGeometry();
        const double area = r_props[CROSS_AREA];

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        const double length = CalculateLength();

        // Let's initialize the cl values
        VectorType strain_vector(1), stress_vector(1);
        MatrixType C(1,1);
        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);
        cl_values.SetConstitutiveMatrix(C);
        SystemSizeBoundedArrayType nodal_values(SystemSize);
        GetNodalValuesVector(nodal_values);

        SystemSizeBoundedArrayType B;

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
             GetFirstDerivativesShapeFunctionsValues(B, length, xi);

            strain_vector[0] = inner_prod(B, nodal_values);

            mConstitutiveLawVector[IP]->CalculateMaterialResponsePK2(cl_values);
            rOutput[IP] = cl_values.GetStressVector()[0] * area;
        }
    } else if (rVariable == AXIAL_STRAIN) {
        const auto &r_props = GetProperties();
        const auto &r_geometry = GetGeometry();

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        const double length = CalculateLength();

        // Let's initialize the cl values
        VectorType strain_vector(1), stress_vector(1);
        MatrixType C(1,1);
        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);
        cl_values.SetConstitutiveMatrix(C);
        SystemSizeBoundedArrayType nodal_values(SystemSize);
        GetNodalValuesVector(nodal_values);

        SystemSizeBoundedArrayType B;

        // Loop over the integration points
        for (SizeType IP = 0; IP < integration_points.size(); ++IP) {
            const double xi = integration_points[IP].X();
             GetFirstDerivativesShapeFunctionsValues(B, length, xi);
            rOutput[IP] = inner_prod(B, nodal_values);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::CalculateOnIntegrationPoints(
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

template<SizeType TNNodes>
int LinearTrussElement2D<TNNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CROSS_AREA)) << "CROSS_AREA not defined in the properties" << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
void LinearTrussElement2D<TNNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TNNodes>
array_1d<double, TNNodes> LinearTrussElement2D<TNNodes>GetBaseShapeFunctions()
{
    array_1d<double, NNodes> N;

    if constexpr (NNodes == 2) {
        N[0] = 0.5 * (1.0 - xi);
        N[1] = 0.5 * (1.0 + xi);
    } else {
        N[0] = 0.5 * xi * (xi - 1.0);
        N[1] = 0.5 * xi * (xi + 1.0);
        N[2] = (1.0 - std::pow(xi, 2));
    }
    return N;
}

/***********************************************************************************/
/***********************************************************************************/

template class LinearTrussElement2D<2>;
template class LinearTrussElement2D<3>;

} // Namespace Kratos
