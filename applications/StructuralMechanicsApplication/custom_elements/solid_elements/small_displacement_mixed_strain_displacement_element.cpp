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
#include "includes/checks.h"

// Application includes
#include "small_displacement_mixed_strain_displacement_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

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
    const SizeType u_size      = dim * n_nodes;
    const SizeType dof_size    = n_nodes * (dim + strain_size);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto &r_node = r_geometry[i_node];
            const SizeType u_block = i_node * dim;
            const SizeType E_block = i_node * strain_size + u_size;

            const IndexType displ_dof = r_node.GetDofPosition(DISPLACEMENT_X);
            const IndexType E_dof     = r_node.GetDofPosition(STRAIN_VECTOR_XX);
            const IndexType E_dof_xy  = r_node.GetDofPosition(STRAIN_VECTOR_XY);

            rResult[u_block    ] = r_node.GetDof(DISPLACEMENT_X,          displ_dof    ).EquationId();
            rResult[u_block + 1] = r_node.GetDof(DISPLACEMENT_Y,          displ_dof + 1).EquationId();
            rResult[E_block    ] = r_node.GetDof(STRAIN_VECTOR_XX,  E_dof        ).EquationId();
            rResult[E_block + 1] = r_node.GetDof(STRAIN_VECTOR_YY,  E_dof + 1    ).EquationId();
            rResult[E_block + 2] = r_node.GetDof(STRAIN_VECTOR_XY,  E_dof_xy     ).EquationId();
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto &r_node = r_geometry[i_node];
            const SizeType u_block = i_node * dim;
            const SizeType E_block = i_node * strain_size + u_size;

            const IndexType displ_dof = r_node.GetDofPosition(DISPLACEMENT_X);
            const IndexType E_dof     = r_node.GetDofPosition(STRAIN_VECTOR_XX);

            rResult[u_block    ] = r_node.GetDof(DISPLACEMENT_X,          displ_dof    ).EquationId();
            rResult[u_block + 1] = r_node.GetDof(DISPLACEMENT_Y,          displ_dof + 1).EquationId();
            rResult[u_block + 2] = r_node.GetDof(DISPLACEMENT_Z,          displ_dof + 2).EquationId();
            rResult[E_block    ] = r_node.GetDof(STRAIN_VECTOR_XX,  E_dof    ).EquationId();
            rResult[E_block + 1] = r_node.GetDof(STRAIN_VECTOR_YY,  E_dof + 1).EquationId();
            rResult[E_block + 2] = r_node.GetDof(STRAIN_VECTOR_ZZ,  E_dof + 2).EquationId();
            rResult[E_block + 3] = r_node.GetDof(STRAIN_VECTOR_XY,  E_dof + 3).EquationId();
            rResult[E_block + 4] = r_node.GetDof(STRAIN_VECTOR_YZ,  E_dof + 4).EquationId();
            rResult[E_block + 5] = r_node.GetDof(STRAIN_VECTOR_XZ,  E_dof + 5).EquationId();
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
    const SizeType u_size      = dim * n_nodes;

    if (rElementalDofList.size() != dof_size) {
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto &r_node = r_geometry[i_node];
            const SizeType u_block = i_node * dim;
            const SizeType E_block = i_node * strain_size + u_size;
            rElementalDofList[u_block    ] = r_node.pGetDof(DISPLACEMENT_X);
            rElementalDofList[u_block + 1] = r_node.pGetDof(DISPLACEMENT_Y);
            rElementalDofList[E_block    ] = r_node.pGetDof(STRAIN_VECTOR_XX);
            rElementalDofList[E_block + 1] = r_node.pGetDof(STRAIN_VECTOR_YY);
            rElementalDofList[E_block + 2] = r_node.pGetDof(STRAIN_VECTOR_XY);
        }
    } else if (dim == 3) {
        for(IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto &r_node = r_geometry[i_node];
            const SizeType u_block = i_node * dim;
            const SizeType E_block = i_node * strain_size + u_size;
            rElementalDofList[u_block    ] = r_node.pGetDof(DISPLACEMENT_X);
            rElementalDofList[u_block + 1] = r_node.pGetDof(DISPLACEMENT_Y);
            rElementalDofList[u_block + 2] = r_node.pGetDof(DISPLACEMENT_Z);
            rElementalDofList[E_block    ] = r_node.pGetDof(STRAIN_VECTOR_XX);
            rElementalDofList[E_block + 1] = r_node.pGetDof(STRAIN_VECTOR_YY);
            rElementalDofList[E_block + 2] = r_node.pGetDof(STRAIN_VECTOR_ZZ);
            rElementalDofList[E_block + 3] = r_node.pGetDof(STRAIN_VECTOR_XY);
            rElementalDofList[E_block + 4] = r_node.pGetDof(STRAIN_VECTOR_YZ);
            rElementalDofList[E_block + 5] = r_node.pGetDof(STRAIN_VECTOR_XZ);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::Initialize(
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Constitutive Law Vector initialisation
        const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
        if (mConstitutiveLawVector.size() != r_integration_points.size()) {
            mConstitutiveLawVector.resize(r_integration_points.size());
        }

        // Initialize material
        InitializeMaterial();

    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::InitializeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
        if (mConstitutiveLawVector[point_number]->RequiresInitializeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const auto &r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

        // Compute U and E
        GetNodalDoFsVectors(this_kinematic_variables.ElementSizeNodalDisplacementsVector, this_kinematic_variables.ElementSizeStrainVector);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(r_geometry, r_properties, rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &cl_options = cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);
        cl_values.SetStressVector(this_constitutive_variables.StressVector);
        cl_values.SetConstitutiveMatrix(this_constitutive_variables.D);

        const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

            noalias(this_constitutive_variables.StrainVector) = this_kinematic_variables.EquivalentStrain;

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, cl_values, point_number, r_integration_points);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->InitializeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_Cauchy);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::FinalizeSolutionStep(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    bool required = false;
    for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        if (mConstitutiveLawVector[point_number]->RequiresFinalizeMaterialResponse()) {
            required = true;
            break;
        }
    }
    if (required) {
        const auto &r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType dimension = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);

        // Compute U and E
        GetNodalDoFsVectors(this_kinematic_variables.ElementSizeNodalDisplacementsVector, this_kinematic_variables.ElementSizeStrainVector);

        ConstitutiveVariables this_constitutive_variables(strain_size);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters cl_values(r_geometry, r_properties, rCurrentProcessInfo);

        // Set constitutive law flags:
        Flags &cl_options = cl_values.GetOptions();
        cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        cl_values.SetStrainVector(this_constitutive_variables.StrainVector);
        cl_values.SetStressVector(this_constitutive_variables.StressVector);
        cl_values.SetConstitutiveMatrix(this_constitutive_variables.D);

        const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, mThisIntegrationMethod);

            noalias(this_constitutive_variables.StrainVector) = this_kinematic_variables.EquivalentStrain;

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables,  this_constitutive_variables, cl_values, point_number, r_integration_points);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(cl_values, ConstitutiveLaw::StressMeasure_Cauchy);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    auto &r_props          = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;
    const double tau = GetStabilizationFactor();

    // Check RHS size
    if (rRHS.size() != matrix_size) {
        rRHS.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLHS.size1() != matrix_size || rLHS.size2() != matrix_size) {
        rLHS.resize(matrix_size, matrix_size, false);
    }
    rLHS.clear();
    rRHS.clear();

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.ElementSizeNodalDisplacementsVector, kinematic_variables.ElementSizeStrainVector);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
    cons_law_values.SetStressVector(constitutive_variables.StrainVector);
    cons_law_values.SetConstitutiveMatrix(constitutive_variables.D);

    Vector RHSu(dim * n_nodes);
    Vector RHSe(strain_size * n_nodes);
    RHSu.clear();
    RHSe.clear();

    Matrix K(dim * n_nodes, dim * n_nodes);
    Matrix G(n_nodes * strain_size, dim * n_nodes);
    Matrix M(n_nodes * strain_size, n_nodes * strain_size);
    Matrix Q(dim * n_nodes, n_nodes * strain_size);
    K.clear();
    G.clear();
    Q.clear();
    M.clear();

    const auto& r_integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);
    SizeType n_gauss = r_integration_points.size();

    Matrix D0(strain_size, strain_size);
    mConstitutiveLawVector[0]->CalculateValue(cons_law_values, CONSTITUTIVE_MATRIX, D0);

    // Gauss IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(mThisIntegrationMethod), i_gauss);

        CalculateKinematicVariables(kinematic_variables, i_gauss, mThisIntegrationMethod);

        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        noalias(constitutive_variables.StrainVector) = kinematic_variables.EquivalentStrain;

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables,  constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

        // Now we add the body forces
        for (IndexType i = 0; i < n_nodes; ++i) {
            const SizeType index = dim * i;
            for (IndexType j = 0; j < dim; ++j) {
                RHSu[index + j] += w_gauss * kinematic_variables.N[i] * body_force[j];
            }
        }

        // Contributions to the RHS
        noalias(RHSu) -= w_gauss * prod(trans(kinematic_variables.B), constitutive_variables.StressVector);
        const Vector strain_diff = prod(kinematic_variables.B, kinematic_variables.ElementSizeNodalDisplacementsVector) - kinematic_variables.EquivalentStrain;
        noalias(RHSe) -= w_gauss * prod(trans(kinematic_variables.N_epsilon), Vector(prod(D0, strain_diff)));

        // Contributions to the LHS
        noalias(K) += w_gauss * prod(trans(kinematic_variables.B), Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
        noalias(Q) += w_gauss * prod(trans(kinematic_variables.B), Matrix(prod(constitutive_variables.D, kinematic_variables.N_epsilon)));

        noalias(M) += w_gauss * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(D0, kinematic_variables.N_epsilon)));
        noalias(G) += w_gauss * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(D0, kinematic_variables.B)));
    }

    K *= tau;
    Q *= (1.0 - tau);
    M *= (tau - 1.0);
    G *= (1.0 - tau);

    AssembleRHS(rRHS, RHSu, RHSe);
    AssembleLHS(rLHS, K, Q, M, G);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    auto &r_props          = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;
    const double tau = GetStabilizationFactor();

    // Check LHS size
    if (rLHS.size1() != matrix_size || rLHS.size2() != matrix_size) {
        rLHS.resize(matrix_size, matrix_size, false);
    }
    rLHS.clear();

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.ElementSizeNodalDisplacementsVector, kinematic_variables.ElementSizeStrainVector);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
    cons_law_values.SetStressVector(constitutive_variables.StrainVector);
    cons_law_values.SetConstitutiveMatrix(constitutive_variables.D);

    Matrix K(dim * n_nodes, dim * n_nodes);
    Matrix G(n_nodes * strain_size, dim * n_nodes);
    Matrix M(n_nodes * strain_size, n_nodes * strain_size);
    Matrix Q(dim * n_nodes, n_nodes * strain_size);
    K.clear();
    G.clear();
    Q.clear();
    M.clear();

    const auto& r_integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);
    SizeType n_gauss = r_integration_points.size();

    Matrix D0(strain_size, strain_size);
    mConstitutiveLawVector[0]->CalculateValue(cons_law_values, CONSTITUTIVE_MATRIX, D0);

    // Gauss IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(mThisIntegrationMethod), i_gauss);

        CalculateKinematicVariables(kinematic_variables, i_gauss, mThisIntegrationMethod);

        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        noalias(constitutive_variables.StrainVector) = kinematic_variables.EquivalentStrain;

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables,  constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

        // Contributions to the LHS
        noalias(K) += w_gauss * prod(trans(kinematic_variables.B), Matrix(prod(constitutive_variables.D, kinematic_variables.B)));
        noalias(Q) += w_gauss * prod(trans(kinematic_variables.B), Matrix(prod(constitutive_variables.D, kinematic_variables.N_epsilon)));

        noalias(M) += w_gauss * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(D0, kinematic_variables.N_epsilon)));
        noalias(G) += w_gauss * prod(trans(kinematic_variables.N_epsilon), Matrix(prod(D0, kinematic_variables.B)));
    }

    K *= tau;
    Q *= (1.0 - tau);
    M *= (tau - 1.0);
    G *= (1.0 - tau);
    AssembleLHS(rLHS, K, Q, M, G);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    auto &r_props          = GetProperties();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;

    // Check RHS size
    if (rRHS.size() != matrix_size) {
        rRHS.resize(matrix_size, false);
    }
    rRHS.clear();

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables(strain_size, dim, n_nodes);

    // Compute U and E
    GetNodalDoFsVectors(kinematic_variables.ElementSizeNodalDisplacementsVector, kinematic_variables.ElementSizeStrainVector);

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables(strain_size);
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, r_props, rProcessInfo);
    auto& r_cl_options = cons_law_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    cons_law_values.SetStrainVector(constitutive_variables.StrainVector);
    cons_law_values.SetStressVector(constitutive_variables.StrainVector);
    cons_law_values.SetConstitutiveMatrix(constitutive_variables.D);

    const auto& r_integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);
    SizeType n_gauss = r_integration_points.size();

    Vector RHSu(dim * n_nodes);
    Vector RHSe(strain_size * n_nodes);
    RHSu.clear();
    RHSe.clear();

    Matrix D0(strain_size, strain_size);
    mConstitutiveLawVector[0]->CalculateValue(cons_law_values, CONSTITUTIVE_MATRIX, D0);


    // IP loop
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {

        const auto body_force = GetBodyForce(r_geometry.IntegrationPoints(mThisIntegrationMethod), i_gauss);

        CalculateKinematicVariables(kinematic_variables, i_gauss, mThisIntegrationMethod);

        double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();
        if (dim == 2 && r_props.Has(THICKNESS))
            w_gauss *= r_props[THICKNESS];

        noalias(constitutive_variables.StrainVector) = kinematic_variables.EquivalentStrain;

        // Calculate the constitutive response with the equivalent stabilized strain
        CalculateConstitutiveVariables(kinematic_variables,  constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);

        // Now we add the body forces
        for (IndexType i = 0; i < n_nodes; ++i) {
            const SizeType index = dim * i;
            for (IndexType j = 0; j < dim; ++j) {
                RHSu[index + j] += w_gauss * kinematic_variables.N[i] * body_force[j];
            }
        }

        // Contributions to the RHS
        noalias(RHSu) -= w_gauss * prod(trans(kinematic_variables.B), constitutive_variables.StressVector);
        const Vector strain_diff = prod(kinematic_variables.B, kinematic_variables.ElementSizeNodalDisplacementsVector) - kinematic_variables.EquivalentStrain;
        noalias(RHSe) -= w_gauss * prod(trans(kinematic_variables.N_epsilon), Vector(prod(D0, strain_diff)));
    }

    AssembleRHS(rRHS, RHSu, RHSe);

    KRATOS_CATCH("")
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

    KRATOS_CATCH("")
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
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    rN_Epsilon.clear();

    if (dim == 2) {
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            const double N_i = rN[i];
            rN_Epsilon(0, initial_index    ) = N_i;
            rN_Epsilon(1, initial_index + 1) = N_i;
            rN_Epsilon(2, initial_index + 2) = N_i;
        }
    } else { // dim = 3
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            const double N_i = rN[i];
            rN_Epsilon(0, initial_index    ) = N_i;
            rN_Epsilon(1, initial_index + 1) = N_i;
            rN_Epsilon(2, initial_index + 2) = N_i;
            rN_Epsilon(3, initial_index + 3) = N_i;
            rN_Epsilon(4, initial_index + 4) = N_i;
            rN_Epsilon(5, initial_index + 5) = N_i;
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::GetNodalDoFsVectors(
    Vector &rU,
    Vector &rE)
{
    KRATOS_TRY
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

        const auto& r_nodal_strain = r_geometry[i].FastGetSolutionStepValue(STRAIN_VECTOR);
        index = i * strain_size;
        if (dim == 3) {
            for(IndexType k = 0; k < strain_size; ++k) {
                rE[index + k] = r_nodal_strain[k];
            }
        } else { // 2D
            rE[index    ] = r_nodal_strain[0];
            rE[index + 1] = r_nodal_strain[1];
            rE[index + 2] = r_nodal_strain[3]; // XY of the 6 components vector
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::AssembleRHS(
    Vector &rRHS,
    const Vector &rRHSu,
    const Vector &rRHSe)
{
    KRATOS_TRY

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

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::AssembleLHS(
    Matrix &rLHS,
    const Matrix &rK,
    const Matrix &rQ,
    const Matrix &rM,
    const Matrix &rG)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType displ_size = n_nodes * dim;

    // Assemble K
    for (IndexType i = 0; i < rK.size1(); ++i)
        for (IndexType j = 0; j < rK.size2(); ++j)
            rLHS(i, j) = rK(i, j);

    // Assemble M
    for (IndexType i = 0; i < rM.size1(); ++i)
        for (IndexType j = 0; j < rM.size2(); ++j)
            rLHS(i + displ_size, j + displ_size) = rM(i, j);

    // Assemble Q
    for (IndexType i = 0; i < rQ.size1(); ++i)
        for (IndexType j = 0; j < rQ.size2(); ++j)
            rLHS(i, j + displ_size) = rQ(i, j);

    // Assemble G
    for (IndexType i = 0; i < rG.size1(); ++i)
        for (IndexType j = 0; j < rG.size2(); ++j)
            rLHS(i + displ_size, j) = rG(i, j);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const IntegrationPointsArrayType& IntegrationPoints
    ) const
{
    KRATOS_TRY

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N);

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);
    rValues.SetStrainVector(rThisConstitutiveVariables.StrainVector);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
{
    KRATOS_TRY

    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); // Here the calculations are actually done

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

Vector SmallDisplacementMixedStrainDisplacementElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    KRATOS_TRY

    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const auto aux_body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    Vector body_force(dim);
    for (IndexType d = 0; d < dim; ++d) {
        body_force(d) = aux_body_force(d);
    }
    return body_force;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateKinematicVariables(
    KinematicVariables& rKinVariables,
    const IndexType IP,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = GetGeometry().IntegrationPoints(rIntegrationMethod);
    const double tau = GetStabilizationFactor();

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
    StructuralMechanicsElementUtilities::CalculateB(*this, rKinVariables.DN_DX, rKinVariables.B);

    // Compute N_epsilon
    CalculateN_EpsilonMatrix(rKinVariables.N_epsilon, rKinVariables.N);

    noalias(rKinVariables.EquivalentStrain) = (1.0 - tau) * prod(rKinVariables.N_epsilon, rKinVariables.ElementSizeStrainVector) + tau * prod(rKinVariables.B, rKinVariables.ElementSizeNodalDisplacementsVector);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

int SmallDisplacementMixedStrainDisplacementElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = BaseType::Check(rCurrentProcessInfo);

    // Base check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < r_geometry.size(); i++) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(STRAIN_VECTOR, r_node)
    }

    return check;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    const SizeType n_gauss = r_integration_points.size();

    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        rOutput[i_gauss] = 0.0;
    }

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(mThisIntegrationMethod);

    const SizeType number_of_integration_points = r_integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points, false);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number <number_of_integration_points; ++point_number ) {
            bool value;
            mConstitutiveLawVector[point_number]->GetValue( rVariable, value);
            rOutput[point_number] = value;
        }
    } else {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ++ii ) {
            bool solution;
            solution = mConstitutiveLawVector[ii]->CalculateValue( Values, rVariable, solution);
            rOutput[ii] = solution;
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const auto &r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    const SizeType number_of_integration_points = r_integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points, false);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        GetNodalDoFsVectors(kinematic_variables.ElementSizeNodalDisplacementsVector, kinematic_variables.ElementSizeStrainVector);

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables(strain_size);
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
            noalias(constitutive_variables.StrainVector) = kinematic_variables.EquivalentStrain;

            // Call the constitutive law to update material variables
            if (rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material response
                CalculateConstitutiveVariables(kinematic_variables,  constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material response
                CalculateConstitutiveVariables(kinematic_variables,  constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            noalias(rOutput[i_gauss]) = constitutive_variables.StressVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR || rVariable == STRAIN) {
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
        GetNodalDoFsVectors(kinematic_variables.ElementSizeNodalDisplacementsVector, kinematic_variables.ElementSizeStrainVector);

        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            noalias(rOutput[i_gauss]) = kinematic_variables.EquivalentStrain;
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedStrainDisplacementElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : false,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","STRAIN_VECTOR"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","STRAIN_VECTOR"],
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
            "This element implements a mixed displacement - strain formulation with Variational MultiScales (VMS) stabilization."
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","STRAIN_VECTOR"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","STRAIN_VECTOR"});
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

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::SetValuesOnIntegrationPoints(
    const Variable<int>& rVariable,
    const std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("SmallDisplacementMixedStrainDisplacementElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    const std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("SmallDisplacementMixedStrainDisplacementElement") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
