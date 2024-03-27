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
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/small_displacement_mixed_strain_displacement_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedStrainDisplacementElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementMixedStrainDisplacementElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry     = GetGeometry();
    const SizeType n_nodes     = r_geometry.PointsNumber();
    const SizeType dim         = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType dof_size    = n_nodes * (dim + strain_size);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    const IndexType displ_pos  = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType strain_pos = r_geometry[0].GetDofPosition(NODAL_STRAIN_VECTOR_XX);

    IndexType aux_index = 0;
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X,          displ_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y,          displ_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XX, strain_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YY, strain_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XY, strain_pos + 2).EquationId();
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, displ_pos             ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, displ_pos + 1         ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, displ_pos + 2         ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XX, strain_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YY, strain_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_ZZ, strain_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XY, strain_pos + 3).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YZ, strain_pos + 4).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XZ, strain_pos + 5).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry     = GetGeometry();
    const SizeType n_nodes     = r_geometry.PointsNumber();
    const SizeType dim         = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType dof_size    = n_nodes * block_size;

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * block_size    ] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XX);
            rElementalDofList[i * block_size + 3] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YY);
            rElementalDofList[i * block_size + 4] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XY);
        }
    } else if (dim == 3) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * block_size    ] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * block_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * block_size + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * block_size + 3] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XX);
            rElementalDofList[i * block_size + 4] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YY);
            rElementalDofList[i * block_size + 5] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_ZZ);
            rElementalDofList[i * block_size + 6] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XY);
            rElementalDofList[i * block_size + 7] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YZ);
            rElementalDofList[i * block_size + 8] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XZ);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Integration method initialization
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Constitutive Law Vector initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size()) {
            mConstitutiveLawVector.resize(r_integration_points.size());
        }

        // Initialize material
        InitializeMaterial();

    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto &r_props    = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;
    const double tau = r_props.Has(STABILIZATION_FACTOR) ? r_props[STABILIZATION_FACTOR] : default_stabilization_factor;

    // Check RHS size
    if (rRHS.size() != matrix_size) {
        rRHS.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLHS.size1() != matrix_size || rLHS.size2() != matrix_size) {
        rLHS.resize(matrix_size, matrix_size, false);
    }
    noalias(rLHS) = ZeroMatrix(matrix_size, matrix_size);
    noalias(rRHS) = ZeroVector(matrix_size);

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.NodalDisplacements, kinematic_variables.NodalStrains);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    Vector RHSu(dim * n_nodes), RHSe(strain_size * n_nodes);
    noalias(RHSu) = ZeroVector(dim * n_nodes);
    noalias(RHSe) = ZeroVector(strain_size * n_nodes);

    Matrix K(dim * n_nodes, dim * n_nodes), G(n_nodes * strain_size, dim * n_nodes),
        M(n_nodes * strain_size, n_nodes * strain_size);
    noalias(K) = ZeroMatrix(dim * n_nodes, dim * n_nodes);
    noalias(G) = ZeroMatrix(n_nodes * strain_size, dim * n_nodes);
    noalias(M) = ZeroMatrix(n_nodes * strain_size, n_nodes * strain_size);

    // IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);

        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        // Calculate the constitutive response with  Symmetric gradient of u
        CalculateConstitutiveVariables(kinematic_variables, kinematic_variables.SymmGradientDispl, constitutive_variables,
            cons_law_values, i_gauss, r_geometry.IntegrationPoints(GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);
        const Vector stress_u = constitutive_variables.StressVector;

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables, kinematic_variables.EquivalentStrain, constitutive_variables,
            cons_law_values, i_gauss, r_geometry.IntegrationPoints(GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);
        const Vector stress_epsilon = constitutive_variables.StressVector;

        // Contributions to the RHS
        noalias(RHSe) -= w_gauss * (prod(trans(kinematic_variables.N_epsilon), stress_u) - prod(trans(kinematic_variables.N_epsilon), stress_epsilon));
        noalias(RHSu) -= w_gauss * prod(trans(kinematic_variables.B), stress_epsilon);

        // Now we add the body forces
        for (IndexType i = 0; i < n_nodes; ++i) {
            const SizeType index = dim * i;
            for (IndexType j = 0; j < dim; ++j)
                RHSu[index + j] += w_gauss * kinematic_variables.N[i] * body_force[j];
        }

        // Contributions to the LHS
        noalias(K) += tau * w_gauss * prod(trans(kinematic_variables.B),         Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
        noalias(M) += (tau - 1.0)   * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(constitutive_variables.D, kinematic_variables.N_epsilon)));
        noalias(G) += (1.0 - tau)   * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
    }
    AssembleRHS(rRHS, RHSu, RHSe);
    AssembleLHS(rLHS, K, G, M);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto &r_props    = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;
    const double tau = r_props.Has(STABILIZATION_FACTOR) ? r_props[STABILIZATION_FACTOR] : default_stabilization_factor;

    // Check LHS size
    if (rLHS.size1() != matrix_size || rLHS.size2() != matrix_size) {
        rLHS.resize(matrix_size, matrix_size, false);
    }
    noalias(rLHS) = ZeroMatrix(matrix_size, matrix_size);

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.NodalDisplacements, kinematic_variables.NodalStrains);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    Matrix K(dim * n_nodes, dim * n_nodes), G(n_nodes * strain_size, dim * n_nodes),
        M(n_nodes * strain_size, n_nodes * strain_size);
    noalias(K) = ZeroMatrix(dim * n_nodes, dim * n_nodes);
    noalias(G) = ZeroMatrix(n_nodes * strain_size, dim * n_nodes);
    noalias(M) = ZeroMatrix(n_nodes * strain_size, n_nodes * strain_size);

    // IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables, kinematic_variables.EquivalentStrain, constitutive_variables,
            cons_law_values, i_gauss, r_geometry.IntegrationPoints(GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);

        // Contributions to the LHS
        noalias(K) += tau * w_gauss * prod(trans(kinematic_variables.B),         Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
        noalias(M) += (tau - 1.0)   * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(constitutive_variables.D, kinematic_variables.N_epsilon)));
        noalias(G) += (1.0 - tau)   * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
    }
    AssembleLHS(rLHS, K, G, M);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto &r_props    = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;
    const double tau = r_props.Has(STABILIZATION_FACTOR) ? r_props[STABILIZATION_FACTOR] : default_stabilization_factor;

    // Check RHS size
    if (rRHS.size() != matrix_size) {
        rRHS.resize(matrix_size, false);
    }

    noalias(rRHS) = ZeroVector(matrix_size);

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.NodalDisplacements, kinematic_variables.NodalStrains);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    Vector RHSu(dim * n_nodes), RHSe(strain_size * n_nodes);
    noalias(RHSu) = ZeroVector(dim * n_nodes);
    noalias(RHSe) = ZeroVector(strain_size * n_nodes);

    // IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(GetIntegrationMethod()), i_gauss);

        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        // Calculate the constitutive response with  Symmetric gradient of u
        CalculateConstitutiveVariables(kinematic_variables, kinematic_variables.SymmGradientDispl, constitutive_variables,
            cons_law_values, i_gauss, r_geometry.IntegrationPoints(GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);
        const Vector stress_u = constitutive_variables.StressVector;

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables, kinematic_variables.EquivalentStrain, constitutive_variables,
            cons_law_values, i_gauss, r_geometry.IntegrationPoints(GetIntegrationMethod()), ConstitutiveLaw::StressMeasure_Cauchy);
        const Vector stress_epsilon = constitutive_variables.StressVector;

        // Contributions to the RHS
        noalias(RHSe) -= w_gauss * (prod(trans(kinematic_variables.N_epsilon), stress_u) - prod(trans(kinematic_variables.N_epsilon), stress_epsilon));
        noalias(RHSu) -= w_gauss * prod(trans(kinematic_variables.B), stress_epsilon);

        // Now we add the body forces
        for (IndexType i = 0; i < n_nodes; ++i) {
            const SizeType index = dim * i;
            for (IndexType j = 0; j < dim; ++j)
                RHSu[index + j] += w_gauss * kinematic_variables.N[i] * body_force[j];
        }
    }
    AssembleRHS(rRHS, RHSu, RHSe);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry = GetGeometry();
        const auto& r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector) {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedStrainDisplacementElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateN_EpsilonMatrix(
    Matrix &rN_Epsilon,
    const Vector &rN) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    rN_Epsilon.clear();

    if(dim == 2) {
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            rN_Epsilon(0, initial_index    ) = rN[i];
            rN_Epsilon(1, initial_index + 1) = rN[i];
            rN_Epsilon(2, initial_index + 2) = rN[i];
        }
    } else { // dim = 3
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            rN_Epsilon(0, initial_index    ) = rN[i];
            rN_Epsilon(1, initial_index + 1) = rN[i];
            rN_Epsilon(2, initial_index + 2) = rN[i];
            rN_Epsilon(3, initial_index + 3) = rN[i];
            rN_Epsilon(4, initial_index + 4) = rN[i];
            rN_Epsilon(5, initial_index + 5) = rN[i];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::GetNodalDoFsVectors(
    Vector &rU,
    Vector &rE)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    if (rU.size() != dim*n_nodes)
        rU.resize(dim * n_nodes, false);
    if (rE.size() != strain_size*n_nodes)
        rE.resize(strain_size * n_nodes, false);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const auto& r_displ = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        SizeType index = i * dim;
        for(IndexType k = 0; k < dim; ++k) {
            rU[index + k] = r_displ[k];
        }

        const auto& r_nodal_strain = r_geometry[i].FastGetSolutionStepValue(NODAL_STRAIN_VECTOR);
        index = i * strain_size;
        for(IndexType k = 0; k < strain_size; ++k) {
            rE[index + k] = r_nodal_strain[k];
        }
    }
}

void SmallDisplacementMixedStrainDisplacementElement::AssembleRHS(
    Vector &rRHS,
    const Vector &rRHSu,
    const Vector &rRHSe)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType system_size = n_nodes * (strain_size + dim);
    const SizeType displ_size = n_nodes * dim;

    if (rRHS.size() != system_size)
        rRHS.resize(system_size, false);

    for (IndexType i = 0; i < rRHSu.size(); ++i) {
        rRHS[i] = rRHSu[i];
    }

    for (IndexType i = 0; i < rRHSe.size(); ++i) {
        rRHS[i + displ_size] = rRHSe[i];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::AssembleLHS(
    Matrix &rLHS,
    const Matrix &rK,
    const Matrix &rG,
    const Matrix &rM)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType system_size = n_nodes * (strain_size + dim);
    const SizeType displ_size = n_nodes * dim;

    // Assemble K
    for (IndexType i = 0; i < rK.size1(); ++i)
        for (IndexType j = 0; j < rK.size2(); ++j)
            rLHS(i, j) = rK(i, j);

    // Assemble M
    for (IndexType i = 0; i < rM.size1(); ++i)
        for (IndexType j = 0; j < rM.size2(); ++j)
            rLHS(i + displ_size, j + displ_size) = rM(i, j);

    // Assemble G
    for (IndexType i = 0; i < rG.size1(); ++i)
        for (IndexType j = 0; j < rG.size2(); ++j)
            rLHS(i + displ_size, j) = rG(i, j);

    // Assemble G transposed
    for (IndexType i = 0; i < rG.size2(); ++i)
        for (IndexType j = 0; j < rG.size1(); ++j)
            rLHS(i, j + displ_size) = rG(j, i);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    Vector& rStrainVector,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    ) const
{
    // Here we essentially set the input parameters
    rValues.SetStrainVector(rStrainVector);
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N);

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    Vector& rStrainVector,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rStrainVector, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); // Here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

Vector SmallDisplacementMixedStrainDisplacementElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const auto aux_body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    Vector body_force(dim);
    for (IndexType d = 0; d < dim; ++d) {
        body_force(d) = aux_body_force(d);
    }
    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateKinematicVariables(
    KinematicVariables& rKinVariables,
    const IndexType IP,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // Shape functions
    rKinVariables.N = r_geometry.ShapeFunctionsValues(rKinVariables.N, r_integration_points[IP].Coordinates());

    // Calculate the inverse Jacobian
    GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_integration_points[IP], rKinVariables.J0);
    MathUtils<double>::InvertMatrix(rKinVariables.J0, rKinVariables.InvJ0, rKinVariables.detJ0);

    KRATOS_ERROR_IF(rKinVariables.detJ0 < 0.0) << "Element ID: " << this->Id() << " is inverted. det(J0) = " 
        << rKinVariables.detJ0 << std::endl;

    // Calculate the shape functions gradients
    GeometryUtils::ShapeFunctionsGradients(r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[IP],
        rKinVariables.InvJ0, rKinVariables.DN_DX);

    // Compute B
    CalculateB(rKinVariables.B, rKinVariables.DN_DX);

    // Compute N_epsilon
    CalculateN_EpsilonMatrix(rKinVariables.N_epsilon, rKinVariables.N);

    // Compute the symmetric gradient of the displacements
    noalias(rKinVariables.SymmGradientDispl) = prod(rKinVariables.B, rKinVariables.NodalDisplacements);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rKinVariables);

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateEquivalentStrain(
    KinematicVariables& rKinVars) const
{
    const auto &r_props = GetProperties();
    const double tau = r_props.Has(STABILIZATION_FACTOR) ? r_props[STABILIZATION_FACTOR] : default_stabilization_factor;
    noalias(rKinVars.EquivalentStrain) = (1.0 - tau) * prod(rKinVars.N_epsilon, rKinVars.NodalStrains) +
        tau * prod(rKinVars.B, rKinVars.NodalDisplacements);
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedStrainDisplacementElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = SmallDisplacementMixedStrainDisplacementElement::BaseType::Check(rCurrentProcessInfo);

    // Base check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < r_geometry.size(); i++) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_STRAIN_VECTOR, r_node)
        KRATOS_CHECK_DOF_IN_NODE(NODAL_STRAIN_VECTOR_XX, r_node)

        // to improve................
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable))
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    else
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // const auto& r_geometry = GetGeometry();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // const SizeType n_gauss = r_integration_points.size();
    // if (rOutput.size() != n_gauss) {
    //     rOutput.resize(n_gauss);
    // }

    // if (mConstitutiveLawVector[0]->Has( rVariable)) {
    //     GetValueOnConstitutiveLaw(rVariable, rOutput);
    // } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
    //     // Create and initialize element variables:
    //     const SizeType n_nodes = r_geometry.PointsNumber();
    //     const SizeType dim = r_geometry.WorkingSpaceDimension();
    //     const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Create the constitutive variables and values containers
    //     ConstitutiveVariables constitutive_variables(strain_size);
    //     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    //     auto& r_cons_law_options = cons_law_values.GetOptions();
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


    //     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //         // Calculate kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Call the constitutive law to update material variables
    //         if( rVariable == CAUCHY_STRESS_VECTOR) {
    //             // Compute material response
    //             CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
    //         } else {
    //             // Compute material response
    //             CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
    //         }

    //         // Check sizes and save the output stress
    //         if (rOutput[i_gauss].size() != strain_size) {
    //             rOutput[i_gauss].resize(strain_size, false);
    //         }
    //         rOutput[i_gauss] = constitutive_variables.StressVector;
    //     }
    // } else {
    //     CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedStrainDisplacementElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Quadrilateral2D4", "Tetrahedra3D4","Hexahedra3D8"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","PlaneStress","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization. This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy."
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainDisplacementElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainDisplacementElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
