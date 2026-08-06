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

void ALM3dMortarFrictionlessCondition::CalculateConditionSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLHS,
    const bool ComputeRHS
    )
{
    KRATOS_TRY

    const auto& r_slave_geometry = GetParentGeometry(); // where we perrofm the integration
    const auto& r_master_geometry = GetPairedGeometry();
    const SizeType system_size = GetSystemSize(); // 3 displacements per node (slave and master) + 1 lagrange multiplier per slave node

    if (ComputeLHS) {
        if (rLeftHandSideMatrix.size1() != system_size || rLeftHandSideMatrix.size2() != system_size) {
            rLeftHandSideMatrix.resize(system_size, system_size, false);
        }
        rLeftHandSideMatrix.clear();
    }
    
    if (ComputeRHS) {
        if (rRightHandSideVector.size() != system_size) {
            rRightHandSideVector.resize(system_size, false);
        }
        rRightHandSideVector.clear();
    }
    
    const auto slave_normal = r_slave_geometry.UnitNormal(0);
    const auto master_normal = r_master_geometry.UnitNormal(0);
    
    const auto& r_properties = this->GetProperties();
    
    const IndexType integration_order = r_properties.Has(INTEGRATION_ORDER_CONTACT) ? r_properties.GetValue(INTEGRATION_ORDER_CONTACT) : 2;
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e0;
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    
    // The utility that performs the local clipping of the projected master to the slave plane
    IntegrationUtility integration_utility = IntegrationUtility(integration_order, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);
    
    ConditionArrayListType conditions_points_slave;
    // This fills the conditions_points_slave, which are LOCAL coordinates on the slave geometry
    const bool is_inside = CheckIsolatedElement(rCurrentProcessInfo[DELTA_TIME]) ? false : 
    integration_utility.GetExactIntegration(r_slave_geometry, slave_normal, r_master_geometry, master_normal, conditions_points_slave);
    
    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);
    

    if (is_inside && ((integration_area / r_slave_geometry.Area()) > MinIntegrationAreaRatioTolerance)) {

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

                VectorType shape_functions_slave;
                VectorType shape_functions_master;
                VectorType dual_shape_functions; // Discretization of the LM


                array_1d<double, 3> zero = ZeroVector(3);
                r_slave_geometry.ShapeFunctionsValues(shape_functions_slave, zero);
                KRATOS_WATCH(shape_functions_slave)
                // ...................

            } // if valid segmented geometry

        } // loop over segmented surfaces
    } else {
        this->Set(ISOLATED, true); // We set the corresponding flag
        rLeftHandSideMatrix.clear();
        rRightHandSideVector.clear();
    }

    KRATOS_CATCH("CalculateConditionSystem");
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
    // KRATOS_ERROR << "You are calling to the base class method GetDofList, check your condition definition" << std::endl;
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
