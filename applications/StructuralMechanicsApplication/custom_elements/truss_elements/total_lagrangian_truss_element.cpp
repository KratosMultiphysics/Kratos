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
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_1;

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

    const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
    const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
    const double ref_length = norm_2(ref_axis);
    const double curr_length = norm_2(curr_axis);
    const array_3 current_unit_dir = curr_axis / curr_length;
    const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
    const double ref_area = r_props[CROSS_AREA];
    const double green_lagrange_strain = 0.5 * (Fxx * Fxx - 1.0);
    const double pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    strain_vector[0] = green_lagrange_strain;

    mConstitutiveLawVector[0]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix
    const double gen_axial_force = ref_area * (stress_vector[0] + pre_stress);
    const array_3 axial_force_vector = gen_axial_force * current_unit_dir;

    for (IndexType i = 0; i < Dimension; ++i) {
        rRHS[i]             = axial_force_vector[i];
        rRHS[i + Dimension] = -axial_force_vector[i];
    }

    // Body forces
    const array_3 body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, IntegrationPoints(mThisIntegrationMethod), 0);
    const array_3 body_force_contribution = 0.5 * ref_area * ref_length * body_force; // divided by 2 for each node
    for (IndexType i = 0; i < Dimension; ++i) {
        rRHS[i]             += body_force_contribution[i];
        rRHS[i + Dimension] += body_force_contribution[i];
    }

    bounded_matrix_3 outer_unit, K_m, K_geo, K_t; // the 3x3 matrices
    noalias(outer_unit) = outer_prod(current_unit_dir, current_unit_dir);
    noalias(K_m) = (ref_area * constitutive_matrix(0, 0) / ref_length) * outer_unit;
    noalias(K_geo) = gen_axial_force / curr_length * (IdentityMatrix(3) - outer_unit);

    noalias(K_t) = K_m + K_geo;

    // Assemble K_t into rLHS
    for (IndexType i = 0; i < Dimension; ++i) {
        for (IndexType j = 0; j < Dimension; ++j) {
            rLHS(i, j)                         =  K_t(i, j);
            rLHS(i, j + Dimension)             = -K_t(i, j);
            rLHS(i + Dimension, j)             = -K_t(i, j);
            rLHS(i + Dimension, j + Dimension) =  K_t(i, j);
        }
    }

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

    const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
    const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
    const double ref_length = norm_2(ref_axis);
    const double curr_length = norm_2(curr_axis);
    const array_3 current_unit_dir = curr_axis / curr_length;
    const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
    const double ref_area = r_props[CROSS_AREA];
    const double green_lagrange_strain = 0.5 * (Fxx * Fxx - 1.0);
    const double pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    strain_vector[0] = green_lagrange_strain;

    mConstitutiveLawVector[0]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix
    const double gen_axial_force = ref_area * (stress_vector[0] + pre_stress);
    const array_3 axial_force_vector = gen_axial_force * current_unit_dir;

    bounded_matrix_3 outer_unit, K_m, K_geo, K_t; // the 3x3 matrices
    noalias(outer_unit) = outer_prod(current_unit_dir, current_unit_dir);
    noalias(K_m) = (ref_area * constitutive_matrix(0, 0) / ref_length) * outer_unit;
    noalias(K_geo) = gen_axial_force / curr_length * (IdentityMatrix(3) - outer_unit);

    noalias(K_t) = K_m + K_geo;

    // Assemble K_t into rLHS
    for (IndexType i = 0; i < Dimension; ++i) {
        for (IndexType j = 0; j < Dimension; ++j) {
            rLHS(i, j)                         =  K_t(i, j);
            rLHS(i, j + Dimension)             = -K_t(i, j);
            rLHS(i + Dimension, j)             = -K_t(i, j);
            rLHS(i + Dimension, j + Dimension) =  K_t(i, j);
        }
    }

    KRATOS_CATCH("CalculateLeftHandSide")
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

    const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
    const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
    const double ref_length = norm_2(ref_axis);
    const double curr_length = norm_2(curr_axis);
    const array_3 current_unit_dir = curr_axis / curr_length;
    const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
    const double ref_area = r_props[CROSS_AREA];
    const double green_lagrange_strain = 0.5 * (Fxx * Fxx - 1.0);
    const double pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    // Let's initialize the cl values
    VectorType strain_vector(1), stress_vector(1);
    MatrixType constitutive_matrix(1, 1); // Young modulus

    strain_vector.clear();
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    strain_vector[0] = green_lagrange_strain;

    mConstitutiveLawVector[0]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix
    const double gen_axial_force = ref_area * (stress_vector[0] + pre_stress);
    const array_3 axial_force_vector = gen_axial_force * current_unit_dir;

    for (IndexType i = 0; i < Dimension; ++i) {
        rRHS[i]             = axial_force_vector[i];
        rRHS[i + Dimension] = -axial_force_vector[i];
    }

    // Body forces
    const array_3 body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, IntegrationPoints(mThisIntegrationMethod), 0);
    const array_3 body_force_contribution = 0.5 * ref_area * ref_length * body_force; // divided by 2 for each node
    for (IndexType i = 0; i < Dimension; ++i) {
        rRHS[i]             += body_force_contribution[i];
        rRHS[i + Dimension] += body_force_contribution[i];
    }

    KRATOS_CATCH("CalculateRightHandSide")
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
    const auto& r_geometry = GetGeometry();

    if (rMassMatrix.size1() != SystemSize || rMassMatrix.size2() != SystemSize) {
        rMassMatrix.resize(SystemSize, SystemSize, false);
    }
    rMassMatrix.clear();

    const double A = r_props[CROSS_AREA];
    const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
    const double L0 = norm_2(ref_axis);

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
    rOutput.resize(1);

    if (rVariable == AXIAL_FORCE) {
        const auto& r_props = GetProperties();
        const auto& r_geometry = GetGeometry();

        const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
        const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
        const double ref_length = norm_2(ref_axis);
        const double curr_length = norm_2(curr_axis);
        const array_3 current_unit_dir = curr_axis / curr_length;
        const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
        const double ref_area = r_props[CROSS_AREA];
        const double green_lagrange_strain = 0.5 * (Fxx * Fxx - 1.0);
        const double pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Let's initialize the cl values
        VectorType strain_vector(1), stress_vector(1);
        MatrixType constitutive_matrix(1, 1); // Young modulus

        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);
        cl_values.SetConstitutiveMatrix(constitutive_matrix);

        strain_vector[0] = green_lagrange_strain;

        mConstitutiveLawVector[0]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix
        rOutput[0] = ref_area * (stress_vector[0] + pre_stress);

    } else if (rVariable == AXIAL_STRAIN) {

        const auto& r_geometry = GetGeometry();
        const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
        const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
        const double ref_length = norm_2(ref_axis);
        const double curr_length = norm_2(curr_axis);
        const array_3 current_unit_dir = curr_axis / curr_length;
        const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
        rOutput[0] = 0.5 * (Fxx * Fxx - 1.0);
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
    rOutput.resize(1);

    if (rVariable == PK2_STRESS_VECTOR) {

        const auto& r_props = GetProperties();
        const auto& r_geometry = GetGeometry();

        const array_3 ref_axis = r_geometry[1].GetInitialPosition().Coordinates() - r_geometry[0].GetInitialPosition().Coordinates();
        const array_3 curr_axis = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
        const double ref_length = norm_2(ref_axis);
        const double curr_length = norm_2(curr_axis);
        const array_3 current_unit_dir = curr_axis / curr_length;
        const double Fxx = curr_length / ref_length; // Deformation gradient in the axial direction
        const double green_lagrange_strain = 0.5 * (Fxx * Fxx - 1.0);
        const double pre_stress = r_props.Has(TRUSS_PRESTRESS_PK2) ? r_props[TRUSS_PRESTRESS_PK2] :  0.0;

        ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
        auto &r_cl_options = cl_values.GetOptions();
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS             , true);
        r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Let's initialize the cl values
        VectorType strain_vector(1), stress_vector(1);
        MatrixType constitutive_matrix(1, 1); // Young modulus

        strain_vector.clear();
        cl_values.SetStrainVector(strain_vector);
        cl_values.SetStressVector(stress_vector);
        cl_values.SetConstitutiveMatrix(constitutive_matrix);

        strain_vector[0] = green_lagrange_strain;

        mConstitutiveLawVector[0]->CalculateMaterialResponsePK2(cl_values); // fills stress and const. matrix
        rOutput[0] = ScalarVector(1, stress_vector[0] + pre_stress);
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
