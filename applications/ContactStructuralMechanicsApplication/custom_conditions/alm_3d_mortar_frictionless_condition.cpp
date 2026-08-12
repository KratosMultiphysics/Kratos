// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// External includes

// Project includes

// #include "includes/global_variables.h"
#include "custom_conditions/alm_3d_mortar_frictionless_condition.h"

/* Utilities */
#include "utilities/geometrical_projection_utilities.h"
// #include "utilities/math_utils.h"
#include "custom_utilities/mortar_explicit_contribution_utilities.h"

namespace Kratos
{

/************************************* OPERATIONS **********************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, this->GetParentGeometry().Create(rThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer ALM3dMortarFrictionlessCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties,
    GeometryType::Pointer pMasterGeom) const
{
    return Kratos::make_intrusive<ALM3dMortarFrictionlessCondition>(NewId, pGeom, pProperties, pMasterGeom);
}

/************************************* DESTRUCTOR **********************************/
/***********************************************************************************/

ALM3dMortarFrictionlessCondition::~ALM3dMortarFrictionlessCondition() = default;

//************************** STARTING - ENDING  METHODS ***************************//
/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    // We reset the ISOLATED flag
    this->Set(ISOLATED, false);

    mSlaveNormal.resize(3, false);
    noalias(mSlaveNormal) = GetParentGeometry().UnitNormal(0);

    KRATOS_CATCH("Initialize");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("InitializeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::InitializeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("InitializeNonLinearIteration");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);

    noalias(mSlaveNormal) = GetParentGeometry().UnitNormal(0);

    KRATOS_CATCH("FinalizeSolutionStep");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);

    KRATOS_CATCH("FinalizeNonLinearIteration");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    // Calculate condition system
    CalculateConditionSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("CalculateLocalSystem");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculate condition system
    VectorType aux_right_hand_side_vector = Vector();
    CalculateConditionSystem(rLeftHandSideMatrix, aux_right_hand_side_vector, rCurrentProcessInfo, true, false);
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // Calculate condition system
    MatrixType aux_left_hand_side_matrix = Matrix();
    CalculateConditionSystem(aux_left_hand_side_matrix, rRightHandSideVector, rCurrentProcessInfo, false, true);
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateInterpolationMatrices(
    const Vector &rN_slave,
    const Vector &rN_master,
    const Vector &rN_LMgeom,
    Matrix &rNs,
    Matrix &rNm,
    Vector &rN_LM)
{
    if (rNs.size1() != 3 || rNs.size2() != 3*rN_slave.size()) {
        rNs.resize(3, 3*rN_slave.size(), false);
    }
    if (rNm.size1() != 3 || rNm.size2() != 3*rN_master.size()) {
        rNm.resize(3, 3*rN_master.size(), false);
    }
    if (rN_LM.size() != rN_LMgeom.size()) {
        rN_LM.resize(rN_LMgeom.size(), false);
    }
    rNs.clear();
    rNm.clear();

    for (IndexType i = 0; i < rN_slave.size(); ++i) {
        rNs(0, 3*i)     = rN_slave[i];
        rNs(1, 3*i + 1) = rN_slave[i];
        rNs(2, 3*i + 2) = rN_slave[i];
    }

    for (IndexType i = 0; i < rN_master.size(); ++i) {
        rNm(0, 3*i)     = rN_master[i];
        rNm(1, 3*i + 1) = rN_master[i];
        rNm(2, 3*i + 2) = rN_master[i];
    }

    noalias(rN_LM) = rN_LMgeom;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::GetValuesVector(
    Vector &rValues,
    int Step) const
{
    const SizeType system_size = GetSystemSize();
    const auto& r_slave_geometry = GetParentGeometry();
    const auto& r_master_geometry = GetPairedGeometry();

    if (rValues.size() != system_size) {
        rValues.resize(system_size, false);
    }

    rValues[0] = r_slave_geometry[0].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, Step);
    rValues[1] = r_slave_geometry[1].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, Step);
    rValues[2] = r_slave_geometry[2].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, Step);

    // todo....
}

/***********************************************************************************/
/***********************************************************************************/
bool ALM3dMortarFrictionlessCondition::FastProjectDirection(
    const GeometryType& rGeometryToProject,
    const PointType& rPointToProject,
    PointType& rPointProjected,
    const array_1d<double,3>& rNormal
    )
{
    KRATOS_TRY

    const double zero_tolerance = CheckThresholdCoefficient;
    const double slave_normal_norm = norm_2(rNormal);

    KRATOS_ERROR_IF(slave_normal_norm <= zero_tolerance)
        << "The projection normal is null in ALM3dMortarFrictionlessCondition.";

    // 1. Compute master triangle face normal
    const array_1d<double, 3> edge_01 = rGeometryToProject[1].Coordinates() - rGeometryToProject[0].Coordinates();
    const array_1d<double, 3> edge_02 = rGeometryToProject[2].Coordinates() - rGeometryToProject[0].Coordinates();
    array_1d<double, 3> master_normal;
    MathUtils<double>::CrossProduct(master_normal, edge_01, edge_02);

    const double master_normal_norm = norm_2(master_normal);
    const double edge_scale = std::max({norm_2(edge_01), norm_2(edge_02),
        norm_2(rGeometryToProject[2].Coordinates() - rGeometryToProject[1].Coordinates())});

    KRATOS_ERROR_IF(edge_scale <= zero_tolerance || master_normal_norm <= zero_tolerance * edge_scale * edge_scale)
        << "The projected geometry is a degenerate triangle in ALM3dMortarFrictionlessCondition.";

    // 2. Normalized projection direction (along slave normal)
    const array_1d<double, 3> proj_dir = rNormal / slave_normal_norm;

    // 3. Ray-plane intersection denominator check
    const double denom = inner_prod(proj_dir, master_normal);

    // Fail projection if ray and master face are parallel
    if (std::abs(denom) <= zero_tolerance * master_normal_norm) {
        noalias(rPointProjected.Coordinates()) = rPointToProject.Coordinates();
        return false;
    }

    // 4. Compute exact ray distance alpha to master plane
    const array_1d<double, 3> point_to_plane = rGeometryToProject[0].Coordinates() - rPointToProject.Coordinates();
    const double alpha = inner_prod(point_to_plane, master_normal) / denom;

    // Set exact projected point on master surface
    noalias(rPointProjected.Coordinates()) = rPointToProject.Coordinates() + alpha * proj_dir;

    // 5. Verify if projected point lies inside master triangle boundaries
    array_1d<double, 3> aux_local_coords;
    const double tol = rGeometryToProject.Area() * 1.0e-8;
    const bool is_inside_master = rGeometryToProject.IsInside(rPointProjected.Coordinates(), aux_local_coords, tol);

    return is_inside_master;

    KRATOS_CATCH("ALM3dMortarFrictionlessCondition::FastProjectDirection")
}


/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY

    const auto& r_slave_geometry = GetParentGeometry(); // where we perform the integration
    const auto& r_master_geometry = GetPairedGeometry();
    const SizeType system_size = GetSystemSize(); // 3 displacements per node (slave and master) + 1 lagrange multiplier per slave node

    if (rLeftHandSideMatrix.size1() != system_size || rLeftHandSideMatrix.size2() != system_size) {
        rLeftHandSideMatrix.resize(system_size, system_size, false);
    }
    rLeftHandSideMatrix.clear();


    if (rRightHandSideVector.size() != system_size) {
        rRightHandSideVector.resize(system_size, false);
    }
    rRightHandSideVector.clear();

    // const double scale_factor = rCurrentProcessInfo.Has(SCALE_FACTOR) ? rCurrentProcessInfo[SCALE_FACTOR] : 1.0e12;
    // const double penalty_factor = rCurrentProcessInfo.Has(INITIAL_PENALTY) ? rCurrentProcessInfo[INITIAL_PENALTY] : 1.0e12;
    const double scale_factor = 1E11;
    const double penalty_factor = 1E12;

    VectorType slave_normal(3);
    VectorType master_normal(3);
    noalias(slave_normal) = r_slave_geometry.UnitNormal(0);
    noalias(master_normal) = r_master_geometry.UnitNormal(0);

    const auto &r_properties = this->GetProperties();

    const IndexType integration_order = 5;
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e-12;
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;

    // The utility that performs the local clipping of the projected master to the slave plane
    IntegrationUtility integration_utility = IntegrationUtility(integration_order, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    ConditionArrayListType conditions_points_slave;
    // This fills the conditions_points_slave, which are LOCAL coordinates on the slave geometry
    const bool is_segmentation_inside = integration_utility.GetExactIntegration(r_slave_geometry, slave_normal, r_master_geometry, master_normal, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    if (((integration_area / r_slave_geometry.Area()) > MinIntegrationAreaRatioTolerance) && is_segmentation_inside) {

        PointType global_point; // aux variable to store the global coordinates of the integration points
        PointerVector<PointType> points_array(3); // aux variable to store the global coordinates of the clipping points
        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) { // Segmented surfaces
            for (IndexType i_node = 0; i_node < conditions_points_slave[i_geom].size(); ++i_node) {
                // Here we tranform the local coordinates to global coordinates
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            }

            DecompositionType segmented_geom(points_array);

            bool bad_shape = MortarUtilities::HeronCheck(segmented_geom);

            if (!bad_shape) { // The segmented geometry is valid
                // We perform the integration over the segmented geometry
                const auto &r_integration_points_slave = segmented_geom.IntegrationPoints(GetIntegrationMethod());

                VectorType shape_functions_slave; // Slave element shape functions evaluated at the integration point
                VectorType shape_functions_master; // Master element shape functions evaluated at the integration point
                VectorType dual_shape_functions; // LM slave shape functions evaluated at the integration point

                MatrixType slave_local_gradients; // Slave element local gradients evaluated at the integration point
                MatrixType slave_jacobian; // Slave element jacobian evaluated at the integration point
                MatrixType d_n_dslave_displacements; // Derivative of the slave normal with respect to the slave displacements

                PointType segmented_GP_global_point; // Global coordinates of the integration point of the segmented geometry
                PointType segmented_GP_local_point; // Local coordinates of the integration point of the segmented geometry
                PointType slave_GP_local_point; // Local coordinates of the integration point of the segmented geometry on the slave element
                PointType master_GP_global_point; // Global coordinates of the integration point of the segmented geometry on the master element
                PointType master_GP_local_point; // Local coordinates of the integration point of the segmented geometry on the master element

                Matrix Ns(3, system_size); // Interpolation matrix for the slave element
                Matrix Nm(3, system_size); // Interpolation matrix for the master element
                Vector N_LM(system_size); // Interpolation matrix for the Lagrange multipliers in the slave element

                Vector lm_values(3); // Lagrange multiplier values at the nodes
                lm_values[0] = r_slave_geometry[0].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                lm_values[1] = r_slave_geometry[1].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                lm_values[2] = r_slave_geometry[2].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);


                // Loop over the integration points of the segmented geometry
                for (IndexType point_number = 0; point_number < r_integration_points_slave.size(); ++point_number) {
                    segmented_GP_local_point = PointType(r_integration_points_slave[point_number].Coordinates());
                    segmented_geom.GlobalCoordinates(segmented_GP_global_point, segmented_GP_local_point); // from local to global coordinates in the segmented
                    r_slave_geometry.PointLocalCoordinates(slave_GP_local_point, segmented_GP_global_point); // from global to local coordinates in the slave
                    const bool successful_projection = FastProjectDirection(r_master_geometry, segmented_GP_global_point, master_GP_global_point, slave_normal);
                    r_master_geometry.PointLocalCoordinates(master_GP_local_point, master_GP_global_point); // from global to local coordinates in the master

                    r_slave_geometry.ShapeFunctionsLocalGradients(slave_local_gradients, slave_GP_local_point);
                    r_slave_geometry.Jacobian(slave_jacobian, slave_GP_local_point);
                    d_n_dslave_displacements = CalculateSlaveNormalDerivatives(slave_jacobian, slave_local_gradients);

                    const double integration_weight = r_integration_points_slave[point_number].Weight() * segmented_geom.DeterminantOfJacobian(segmented_GP_local_point);

                    const double gap_n = -inner_prod(segmented_GP_global_point - master_GP_global_point, slave_normal); // open contact if gap_n > 0, closed contact if gap_n < 0

                    r_slave_geometry.ShapeFunctionsValues(shape_functions_slave, slave_GP_local_point);
                    // r_slave_geometry.ShapeFunctionsValues(dual_shape_functions, slave_GP_local_point); // NOTE: for now we do not use DUAL shape functions
                    dual_shape_functions.resize(3);
                    dual_shape_functions[0] = 3.0 - 4.0 * slave_GP_local_point[0] - 4.0 * slave_GP_local_point[1];
                    dual_shape_functions[1] = 4.0 * slave_GP_local_point[0] - 1.0;
                    dual_shape_functions[2] = 4.0 * slave_GP_local_point[1] - 1.0;

                    r_master_geometry.ShapeFunctionsValues(shape_functions_master, master_GP_local_point);

                    CalculateInterpolationMatrices(shape_functions_slave, shape_functions_master, dual_shape_functions, Ns, Nm, N_LM);
                    const double interpolated_LM = inner_prod(N_LM, lm_values); // interpolated Lagrange multiplier at the integration point
                    const double augmented_lm = scale_factor * interpolated_LM + penalty_factor * gap_n; // contact pressure at the integration point

                    const bool active_contact = (augmented_lm < 0.0); // active contact if augmented_lm < 0, inactive contact if augmented_lm > 0
                    
                    if (gap_n < 0.0) {
                    //     KRATOS_WATCH("*******************************")
                        // KRATOS_WATCH(Id())
                        // KRATOS_WATCH(r_slave_geometry[0].Coordinates())
                        // KRATOS_WATCH(r_slave_geometry[0].Id())
                        // KRATOS_WATCH(r_slave_geometry[1].Id())
                        // KRATOS_WATCH(r_slave_geometry[2].Id())
                        // KRATOS_WATCH(r_master_geometry[0].Id())
                        // KRATOS_WATCH(r_master_geometry[1].Id())
                        // KRATOS_WATCH(r_master_geometry[2].Id())
                        // KRATOS_WATCH(r_slave_geometry[0].Coordinates())
                        // KRATOS_WATCH(r_slave_geometry[1].Coordinates())
                        // KRATOS_WATCH(r_slave_geometry[2].Coordinates())
                        // KRATOS_WATCH(r_master_geometry[0].Coordinates())
                        // KRATOS_WATCH(r_master_geometry[1].Coordinates())
                        // KRATOS_WATCH(r_master_geometry[2].Coordinates())
                        // KRATOS_WATCH(segmented_GP_global_point)
                        // KRATOS_WATCH(master_GP_global_point)
                        // KRATOS_WATCH(is_segmentation_inside)
                        // KRATOS_WATCH(successful_projection)
                        // KRATOS_WATCH(integration_area / r_slave_geometry.Area())
                        // KRATOS_WATCH(slave_normal)
                        // KRATOS_WATCH(master_normal)
                        // KRATOS_WATCH(gap_n)
                        // KRATOS_WATCH(active_contact)
                        // KRATOS_WATCH(master_GP_local_point)
                        // KRATOS_WATCH(segmented_GP_global_point)
                        // KRATOS_WATCH(augmented_lm)
                        // KRATOS_WATCH("*******************************")
                        // KRATOS_ERROR << "" << std::endl;
                    }
                    else
                    {
                        // KRATOS_WATCH(is_segmentation_inside)
                        // KRATOS_WATCH(successful_projection)
                        // KRATOS_WATCH(integration_area)
                    }

                    if (ComputeRHS) {
                        AddRightHandSideContribution(rRightHandSideVector, scale_factor, penalty_factor, gap_n, interpolated_LM,
                            integration_weight, Ns, Nm, N_LM, slave_normal, active_contact);
                    }
                    if (ComputeLHS) {
                        AddLeftHandSideContribution(rLeftHandSideMatrix, scale_factor, penalty_factor, integration_weight,
                            augmented_lm, Ns, Nm, N_LM, slave_normal,  segmented_GP_global_point - master_GP_global_point, active_contact, d_n_dslave_displacements);
                    }
                } // loop over the integration points of the segmented geometry
            } // if valid segmented geometry
        } // loop over segmented surfaces
    }

    KRATOS_CATCH("CalculateConditionSystem");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::AddRightHandSideContribution(
    VectorType& rRightHandSideVector,
    const double k,
    const double penalty,
    const double gap,
    const double LM,
    const double integration_weight,
    const Matrix& rNs,
    const Matrix& rNm,
    const Vector& rN_LM, // Phi
    const Vector& rNormal, // slave normal vector
    const bool active_contact
)
{
    KRATOS_TRY

    if (active_contact) {
        const double augmented_lm = k * LM + penalty * gap; // contact pressure at the integration point
        // NOTE: In here we change the sign of the RHS contributions since Kratos residual is f_ext - f_int, and the contact contributions are internal forces
        noalias(project(rRightHandSideVector, range(0, 3))) -= integration_weight * k * gap * rN_LM; // LM dofs
        noalias(project(rRightHandSideVector, range(3, 12))) += integration_weight * augmented_lm * prod(trans(rNs), rNormal); // slave displacement dofs
        noalias(project(rRightHandSideVector, range(12, 21))) -= integration_weight * augmented_lm * prod(trans(rNm), rNormal); // master displacement dofs
    } else { // Inactive
        noalias(project(rRightHandSideVector, range(0, 3))) += integration_weight * (k * k / penalty) * LM * rN_LM;
    }

    KRATOS_CATCH("AddRightHandSideContribution");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::AddLeftHandSideContribution(
    MatrixType& rLeftHandSideMatrix,
    const double k,
    const double penalty,
    const double integration_weight,
    const double AugmentedLM,
    const Matrix& rNs,
    const Matrix& rNm,
    const Vector& rN_LM, // Phi
    const Vector& rNormal, // slave normal vector
    const Vector& rDeltaX, // Xs - Xm
    const bool active_contact,
    const Matrix& rDnda_slave
)
{
    KRATOS_TRY

    if (active_contact) {
        const Matrix n_n = outer_prod(rNormal, rNormal); // normal matrix
        const double integration_k = integration_weight * k;
        const double integration_penalty = integration_weight * penalty;

        //
        noalias(project(rLeftHandSideMatrix, range(0, 3), range(3, 12))) -= integration_k * outer_prod(rN_LM, Vector(prod(trans(rNormal), rNs))); // LM - slave
        noalias(project(rLeftHandSideMatrix, range(0, 3), range(3, 12))) -= integration_k * outer_prod(rN_LM, Vector(prod(trans(rDeltaX), rDnda_slave))); // LM - slave
        //
        noalias(project(rLeftHandSideMatrix, range(0, 3), range(12, 21))) += integration_k * outer_prod(rN_LM, Vector(prod(trans(rNormal), rNm))); // LM - master
        //
        noalias(project(rLeftHandSideMatrix, range(3, 12), range(0, 3))) -= integration_k * outer_prod(Vector(prod(trans(rNs), rNormal)), rN_LM); // slave - LM
        //
        noalias(project(rLeftHandSideMatrix, range(3, 12), range(3, 12))) += integration_penalty * prod(trans(rNs), Matrix(prod(n_n, rNs))); // slave - slave
        noalias(project(rLeftHandSideMatrix, range(3, 12), range(3, 12))) -= integration_weight * AugmentedLM * prod(trans(rNs), rDnda_slave); // slave - slave
        noalias(project(rLeftHandSideMatrix, range(3, 12), range(3, 12))) += integration_penalty * outer_prod(Vector(prod(trans(rDeltaX), rDnda_slave)), Vector(prod(trans(rNs), rNormal))); // slave - slave
        //
        noalias(project(rLeftHandSideMatrix, range(3, 12), range(12, 21))) -= integration_penalty * prod(trans(rNs), Matrix(prod(n_n, rNm))); // slave - master
        //
        noalias(project(rLeftHandSideMatrix, range(12, 21), range(0, 3))) += integration_k * outer_prod(Vector(prod(trans(rNm), rNormal)), rN_LM); // master - LM
        //
        noalias(project(rLeftHandSideMatrix, range(12, 21), range(3, 12))) -= integration_penalty * prod(trans(rNm), Matrix(prod(n_n, rNs))); // master - slave
        noalias(project(rLeftHandSideMatrix, range(12, 21), range(3, 12))) -= integration_penalty * outer_prod(Vector(prod(trans(rDeltaX), rDnda_slave)), Vector(prod(trans(rNm), rNormal))); // master - slave
        noalias(project(rLeftHandSideMatrix, range(12, 21), range(3, 12))) += integration_weight * AugmentedLM * prod(trans(rNm), rDnda_slave); // master - slave
        //
        noalias(project(rLeftHandSideMatrix, range(12, 21), range(12, 21))) += integration_penalty * prod(trans(rNm), Matrix(prod(n_n, rNm))); // master - master

    } else { // Inactive
        noalias(project(rLeftHandSideMatrix, range(0, 3), range(0, 3))) -= integration_weight * (k * k / penalty) * outer_prod(rN_LM, rN_LM); // LM dofs
    }

    KRATOS_CATCH("AddLeftHandSideContribution");
}

/***********************************************************************************/
/***********************************************************************************/

bool ALM3dMortarFrictionlessCondition::CheckIsolatedElement(
    const double DeltaTime,
    const bool HalfJump
    )
{
    if (this->Is(ISOLATED))
        return true;

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& CurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const GeometryType& r_current_master = this->GetPairedGeometry();
    const GeometryType& r_current_slave = this->GetParentGeometry();
    const SizeType system_size = this->GetSystemSize(); // 3 displacements per node (slave and master) + 1 lagrange multiplier per slave node

    if (rResult.size() != system_size)
        rResult.resize(system_size, false);

    IndexType index = 0;

    /* ORDER - [ LAMBDA, SLAVE, MASTER ] */

    // Slave Nodes  Lambda Equation IDs
    for (IndexType i_slave = 0; i_slave < r_current_slave.PointsNumber(); ++i_slave) {
        const Node &r_slave_node = r_current_slave[i_slave];
        rResult[index++] = r_slave_node.GetDof(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE).EquationId();
    }

    // Slave Nodes Displacement Equation IDs
    for (IndexType i_slave = 0; i_slave < r_current_slave.PointsNumber(); ++i_slave) {
        const Node &r_slave_node = r_current_slave[i_slave];
        rResult[index++] = r_slave_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = r_slave_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = r_slave_node.GetDof(DISPLACEMENT_Z).EquationId();
    }

    // Master Nodes Displacement Equation IDs
    for (IndexType i_master = 0; i_master < r_current_master.PointsNumber(); ++i_master) { // NOTE: Assuming same number of nodes for master and slave
        const Node &r_master_node = r_current_master[i_master];
        rResult[index++] = r_master_node.GetDof(DISPLACEMENT_X).EquationId();
        rResult[index++] = r_master_node.GetDof(DISPLACEMENT_Y).EquationId();
        rResult[index++] = r_master_node.GetDof(DISPLACEMENT_Z).EquationId();
    }

    KRATOS_CATCH("EquationIdVector");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::GetDofList(
    DofsVectorType& rConditionalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY

    const GeometryType& r_current_master = this->GetPairedGeometry();
    const GeometryType& r_current_slave = this->GetParentGeometry();
    const SizeType system_size = this->GetSystemSize(); // 3 displacements per node (slave and master) + 1 lagrange multiplier per slave node

    if (rConditionalDofList.size() != system_size)
        rConditionalDofList.resize(system_size);

    IndexType index = 0;

    // Slave Nodes Lambda Equation IDs
    for (IndexType i_slave = 0; i_slave < r_current_slave.PointsNumber(); ++i_slave)  {
        const Node &r_slave_node = r_current_slave[i_slave];
        rConditionalDofList[index++] = r_slave_node.pGetDof(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
    }

    // Slave Nodes Displacement Equation IDs
    for (IndexType i_slave = 0; i_slave < r_current_slave.PointsNumber(); ++i_slave) {
        const Node& r_slave_node = r_current_slave[ i_slave ];
        rConditionalDofList[index++] = r_slave_node.pGetDof(DISPLACEMENT_X);
        rConditionalDofList[index++] = r_slave_node.pGetDof(DISPLACEMENT_Y);
        rConditionalDofList[index++] = r_slave_node.pGetDof(DISPLACEMENT_Z);
    }

    // Master Nodes Displacement Equation IDs
    for (IndexType i_master = 0; i_master < r_current_master.PointsNumber(); ++i_master) { 
        const Node &r_master_node = r_current_master[i_master];
        rConditionalDofList[index++] = r_master_node.pGetDof(DISPLACEMENT_X);
        rConditionalDofList[index++] = r_master_node.pGetDof(DISPLACEMENT_Y);
        rConditionalDofList[index++] = r_master_node.pGetDof(DISPLACEMENT_Z);
    }

    KRATOS_CATCH("GetDofList");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>> &rVariable,
    std::vector<array_1d<double, 3>> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

void ALM3dMortarFrictionlessCondition::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    KRATOS_CATCH("CalculateOnIntegrationPoints");
}

/***********************************************************************************/
/***********************************************************************************/

int ALM3dMortarFrictionlessCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check_result = BaseType::Check(rCurrentProcessInfo);

    return check_result;

    KRATOS_CATCH("Check")
}

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
