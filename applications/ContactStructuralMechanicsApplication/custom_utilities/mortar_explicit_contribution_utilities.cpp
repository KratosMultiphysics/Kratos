// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_utilities/mortar_explicit_contribution_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{
template< const SizeType TDim, const SizeType TNumNodes, const FrictionalCase TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
typename MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::MortarConditionMatrices MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContributionOfMortarCondition(
    PairedCondition* pCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const IndexType IntegrationOrder,
    const bool AxisymmetricCase,
    const bool ComputeNodalArea,
    const bool ComputeDualLM,
    const Variable<double>& rAreaVariable
    )
{
    KRATOS_TRY

    // The slave geometry
    GeometryType& r_slave_geometry = pCondition->GetParentGeometry();
    const array_1d<double, 3>& r_normal_slave = pCondition->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables kinematic_variables;

    // Create the dual LM operator
    BoundedMatrix<double, TNumNodes, TNumNodes> Ae;

    // Create the mortar operators
    MortarConditionMatrices this_mortar_condition_matrices;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e0;
    const auto& r_properties = pCondition->GetProperties();
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    IntegrationUtility integration_utility = IntegrationUtility (IntegrationOrder, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    // The master geometry
    GeometryType& r_master_geometry = pCondition->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& r_normal_master = pCondition->GetPairedNormal();

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    // Dual LM
    bool dual_LM = false;

    const double geometry_area = r_slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-5)) {
        IntegrationMethod this_integration_method = pCondition->GetIntegrationMethod();

        // Initialize general variables for the current master element
        kinematic_variables.Initialize();

        // Initialize the mortar operators
        this_mortar_condition_matrices.Initialize();

        if (ComputeDualLM) {
            const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
            dual_LM = ExplicitCalculateAe(r_slave_geometry, kinematic_variables, conditions_points_slave, Ae, this_integration_method, axisymmetric_coefficient);
        }

        PointerVector< PointType > points_array(TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            PointType global_point;
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            }

            DecompositionType decomp_geom( points_array );

            bool bad_shape;
            if constexpr (TDim == 2) {
                bad_shape = MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * CheckThresholdCoefficient);
            } else { 
                bad_shape = MortarUtilities::HeronCheck(decomp_geom);
            }

            if (!bad_shape) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                PointType local_point_parent, gp_global;
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates
                    const PointType local_point_decomp = PointType(integration_points_slave[point_number].Coordinates());
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    ExplicitCalculateKinematics(pCondition, kinematic_variables, Ae, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
                    const double integration_weight = integration_points_slave[point_number].Weight() * axisymmetric_coefficient;

                    this_mortar_condition_matrices.CalculateMortarOperators(kinematic_variables, integration_weight);
                }
            }
        }

        // Setting the weighted gap
        // Mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = this_mortar_condition_matrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = this_mortar_condition_matrices.MOperator;

        // Computing contribution of the rAreaVariable
        if (ComputeNodalArea && dual_LM) {
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                double& r_nodal_area = r_slave_geometry[i_node].GetValue(rAreaVariable);
                AtomicAdd(r_nodal_area, DOperator(i_node, i_node));
            }
        }

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry);
        const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry);

        const BoundedMatrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(DOperator, x1) - prod(MOperator, x2);

        array_1d<double, TDim> aux_normal, aux_array;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                aux_normal[i_dim] = r_normal[i_dim];
            }
            noalias(aux_array) = row(D_x1_M_x2, i_node);

            double& r_weighted_gap = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);

            AtomicAdd(r_weighted_gap, inner_prod(aux_array, - aux_normal));
        }

        // We reset the flag
        pCondition->Set(ISOLATED, false);
    } else {
        // We set the flag
        pCondition->Set(ISOLATED, true);
    }

    return this_mortar_condition_matrices;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes, const FrictionalCase TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
typename MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::MortarConditionMatrices MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::AddExplicitContributionOfMortarFrictionalCondition(
    PairedCondition* pCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const MortarOperator<TNumNodes, TNumNodesMaster>& rPreviousMortarOperators,
    const IndexType IntegrationOrder,
    const bool AxisymmetricCase,
    const bool ComputeNodalArea,
    const bool ComputeDualLM,
    const Variable<double>& rAreaVariable,
    const bool ConsiderObjetiveFormulation
    )
{
    KRATOS_TRY;

    // The slave geometry
    GeometryType& r_slave_geometry = pCondition->GetParentGeometry();
    const array_1d<double, 3>& r_normal_slave = pCondition->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables kinematic_variables;

    // Create the dual LM operator
    BoundedMatrix<double, TNumNodes, TNumNodes> Ae;

    // Create the mortar operators
    MortarConditionMatrices this_mortar_condition_matrices;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e0;
    const auto& r_properties = pCondition->GetProperties();
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    IntegrationUtility integration_utility = IntegrationUtility (IntegrationOrder, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    // The master geometry
    const GeometryType& r_master_geometry = pCondition->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& r_normal_master = pCondition->GetPairedNormal();

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    // Dual LM
    bool dual_LM = false;

    const double geometry_area = r_slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-5)) {
        IntegrationMethod this_integration_method = pCondition->GetIntegrationMethod();

        // Initialize general variables for the current master element
        kinematic_variables.Initialize();

        // Initialize the mortar operators
        this_mortar_condition_matrices.Initialize();

        if (ComputeDualLM) {
            const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
            dual_LM = ExplicitCalculateAe(r_slave_geometry, kinematic_variables, conditions_points_slave, Ae, this_integration_method, axisymmetric_coefficient);
        }

        // Iterate over segmented conditions
        PointerVector<PointType> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of pCondition points
        PointType global_point;
        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
            }

            DecompositionType decomp_geom( points_array );

            bool bad_shape;
            if constexpr (TDim == 2) {
                bad_shape = MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * CheckThresholdCoefficient);
            } else { 
                bad_shape = MortarUtilities::HeronCheck(decomp_geom);
            }

            if (!bad_shape) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                PointType local_point_parent, gp_global;
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates
                    const PointType local_point_decomp = PointType(integration_points_slave[point_number].Coordinates());
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    ExplicitCalculateKinematics(pCondition, kinematic_variables, Ae, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
                    const double integration_weight = integration_points_slave[point_number].Weight() * axisymmetric_coefficient;

                    this_mortar_condition_matrices.CalculateMortarOperators(kinematic_variables, integration_weight);
                }
            }
        }

        // Setting the weighted gap
        // Mortar condition matrices - DOperator and MOperator
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = this_mortar_condition_matrices.DOperator;
        const BoundedMatrix<double, TNumNodes, TNumNodesMaster>& MOperator = this_mortar_condition_matrices.MOperator;

        // Computing contribution of the rAreaVariable
        if (ComputeNodalArea && dual_LM) {
            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                double& r_nodal_area = r_slave_geometry[i_node].GetValue(rAreaVariable);
                AtomicAdd(r_nodal_area, DOperator(i_node, i_node));
            }
        }

        // Current coordinates
        const BoundedMatrix<double, TNumNodes, TDim> x1 = MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry);
        const BoundedMatrix<double, TNumNodesMaster, TDim> x2 = MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry);

        const BoundedMatrix<double, TNumNodes, TDim> D_x1_M_x2 = prod(DOperator, x1) - prod(MOperator, x2);

        // Setting the weighted slip
        // The increment of time
        const double delta_time = rCurrentProcessInfo.Has(DELTA_TIME) ? rCurrentProcessInfo[DELTA_TIME] : 1.0;

        // The estimation of the slip time derivative
        BoundedMatrix<double, TNumNodes, TDim> slip_time_derivative;
        const bool objective_formulation = ConsiderObjetiveFormulation ? true : pCondition->IsDefined(MODIFIED) ? pCondition->IsNot(MODIFIED) : true;

        if (objective_formulation) {
            // Delta mortar condition matrices - DOperator and MOperator
            const BoundedMatrix<double, TNumNodes, TNumNodes> DeltaDOperator = DOperator - rPreviousMortarOperators.DOperator;
            const BoundedMatrix<double, TNumNodes, TNumNodesMaster> DeltaMOperator = MOperator - rPreviousMortarOperators.MOperator;

            // Delta objetive gap and slip
            noalias(slip_time_derivative)  = (prod(DeltaDOperator, x1) - prod(DeltaMOperator, x2))/delta_time;
        } else {
            const BoundedMatrix<double, TNumNodes, TDim> delta_x1 = x1 - MortarUtilities::GetCoordinates<TDim,TNumNodes>(r_slave_geometry, false, 1);
            const BoundedMatrix<double, TNumNodesMaster, TDim> delta_x2 = x2 - MortarUtilities::GetCoordinates<TDim,TNumNodesMaster>(r_master_geometry, false, 1);

            // Delta non-objetive gap and slip
            noalias(slip_time_derivative)  = - (prod(DOperator, delta_x1) - prod(MOperator, delta_x2))/delta_time;
        }

        array_1d<double, TDim> normal, aux_array;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const array_1d<double, 3>& r_normal = r_slave_geometry[i_node].FastGetSolutionStepValue(NORMAL);
            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                normal[i_dim] = r_normal[i_dim];
            }
            noalias(aux_array) = row(D_x1_M_x2, i_node);

            double& r_weighted_gap = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_GAP);

            AtomicAdd(r_weighted_gap, inner_prod(aux_array, - normal));

            // We compute the tangent component
            const array_1d<double, TDim>& r_slip_time_derivative_node = row(slip_time_derivative, i_node);
            const array_1d<double, TDim> slip_node = delta_time * (r_slip_time_derivative_node - inner_prod(normal, r_slip_time_derivative_node) * normal);

            // The weighted slip
            array_1d<double, 3>& r_weighted_slip = r_slave_geometry[i_node].FastGetSolutionStepValue(WEIGHTED_SLIP);

            for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
                AtomicAdd(r_weighted_slip[i_dim], slip_node[i_dim]);
            }
        }

        // We reset the flag
        pCondition->Set(ISOLATED, false);
    } else {
        // We set the flag
        pCondition->Set(ISOLATED, true);
    }

    return this_mortar_condition_matrices;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes, const FrictionalCase TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
void MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::ComputeNodalArea(
    PairedCondition* pCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const Variable<double>& rAreaVariable,
    const IndexType IntegrationOrder,
    const bool AxisymmetricCase
    )
{
    MortarOperator<TNumNodes, TNumNodesMaster> mortar_operator;
    mortar_operator.Initialize();
    ComputePreviousMortarOperators(pCondition, rCurrentProcessInfo, mortar_operator, IntegrationOrder, AxisymmetricCase, true, true, rAreaVariable);
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes, const FrictionalCase TFrictional, const bool TNormalVariation, const SizeType TNumNodesMaster>
bool MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::ComputePreviousMortarOperators(
    PairedCondition* pCondition,
    const ProcessInfo& rCurrentProcessInfo,
    MortarOperator<TNumNodes, TNumNodesMaster>& rPreviousMortarOperators,
    const IndexType IntegrationOrder,
    const bool AxisymmetricCase,
    const bool ComputeNodalArea,
    const bool ComputeDualLM,
    const Variable<double>& rAreaVariable
    )
{
    // We "save" the mortar operator for the next step
    // The slave geometry
    GeometryType& r_slave_geometry = pCondition->GetParentGeometry();
    const array_1d<double, 3>& r_normal_slave = pCondition->GetValue(NORMAL);

    // Create and initialize condition variables
    GeneralVariables kinematic_variables;

    // Create the dual LM operator
    BoundedMatrix<double, TNumNodes, TNumNodes> Ae;

    // We call the exact integration utility
    const double distance_threshold = rCurrentProcessInfo.Has(DISTANCE_THRESHOLD) ? rCurrentProcessInfo[DISTANCE_THRESHOLD] : 1.0e24;
    const double zero_tolerance_factor = rCurrentProcessInfo.Has(ZERO_TOLERANCE_FACTOR) ? rCurrentProcessInfo[ZERO_TOLERANCE_FACTOR] : 1.0e0;
    const auto& r_properties = pCondition->GetProperties();
    const bool consider_tessellation = r_properties.Has(CONSIDER_TESSELLATION) ? r_properties[CONSIDER_TESSELLATION] : false;
    IntegrationUtility integration_utility = IntegrationUtility (IntegrationOrder, distance_threshold, 0, zero_tolerance_factor, consider_tessellation);

    // The master geometry
    const GeometryType& r_master_geometry = pCondition->GetPairedGeometry();

    // The normal of the master condition
    const array_1d<double, 3>& r_normal_master = pCondition->GetPairedNormal();

    // Reading integration points
    ConditionArrayListType conditions_points_slave;
    const bool is_inside = integration_utility.GetExactIntegration(r_slave_geometry, r_normal_slave, r_master_geometry, r_normal_master, conditions_points_slave);

    double integration_area;
    integration_utility.GetTotalArea(r_slave_geometry, conditions_points_slave, integration_area);

    // Dual LM
    bool dual_LM = false;

    const double geometry_area = r_slave_geometry.Area();
    if (is_inside && ((integration_area/geometry_area) > 1.0e-5)) {
        IntegrationMethod this_integration_method = pCondition->GetIntegrationMethod();

        // Initialize general variables for the current master element
        kinematic_variables.Initialize();

        // Initialize the mortar operators
        rPreviousMortarOperators.Initialize();

        if (ComputeDualLM) {
            const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
            dual_LM = ExplicitCalculateAe(r_slave_geometry, kinematic_variables, conditions_points_slave, Ae, this_integration_method, axisymmetric_coefficient);
        }

        // Iterate over segmented conditions
        PointerVector<PointType> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
        array_1d<BelongType, TDim> belong_array;
        PointType global_point;
        for (IndexType i_geom = 0; i_geom < conditions_points_slave.size(); ++i_geom) {
            for (IndexType i_node = 0; i_node < TDim; ++i_node) {
                r_slave_geometry.GlobalCoordinates(global_point, conditions_points_slave[i_geom][i_node]);
                points_array(i_node) = Kratos::make_shared<PointType>(PointType(global_point));
                belong_array[i_node] = conditions_points_slave[i_geom][i_node].GetBelong();
            }

            DecompositionType decomp_geom( points_array );

            bool bad_shape;
            if constexpr (TDim == 2) {
                bad_shape = MortarUtilities::LengthCheck(decomp_geom, r_slave_geometry.Length() * CheckThresholdCoefficient);
            } else { 
                bad_shape = MortarUtilities::HeronCheck(decomp_geom);
            }

            if (!bad_shape) {
                const GeometryType::IntegrationPointsArrayType& integration_points_slave = decomp_geom.IntegrationPoints( this_integration_method );

                // Integrating the mortar operators
                for ( IndexType point_number = 0; point_number < integration_points_slave.size(); ++point_number ) {
                    // We compute the local coordinates

                    const PointType local_point_decomp = PointType(integration_points_slave[point_number].Coordinates());
                    PointType local_point_parent;
                    PointType gp_global;
                    decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                    r_slave_geometry.PointLocalCoordinates(local_point_parent, gp_global);

                    // Calculate the kinematic variables
                    ExplicitCalculateKinematics(pCondition, kinematic_variables, Ae, r_normal_master, local_point_decomp, local_point_parent, decomp_geom, dual_LM);

                    const double axisymmetric_coefficient = AxisymmetricCase ? AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(pCondition, kinematic_variables.NSlave) : 1.0;
                    const double integration_weight = integration_points_slave[point_number].Weight() * axisymmetric_coefficient;

                    rPreviousMortarOperators.CalculateMortarOperators(kinematic_variables, integration_weight);
                }
            }
        }
    }

    // Computing contribution of the rAreaVariable
    if (ComputeNodalArea && dual_LM) {
        const BoundedMatrix<double, TNumNodes, TNumNodes>& DOperator = rPreviousMortarOperators.DOperator;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            double& r_nodal_area = r_slave_geometry[i_node].GetValue(rAreaVariable);
            AtomicAdd(r_nodal_area, DOperator(i_node, i_node));
        }
    }

    return dual_LM;
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
bool MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::ExplicitCalculateAe(
    const GeometryType& rSlaveGeometry,
    GeneralVariables& rVariables,
    const ConditionArrayListType& rConditionsPointsSlave,
    BoundedMatrix<double, TNumNodes, TNumNodes>& rAe,
    const IntegrationMethod& rIntegrationMethod,
    const double AxiSymCoeff
    )
{
    // Initialize general variables for the current master element
    rVariables.Initialize();

    // We initialize the Ae components
    AeData Ae_data;
    Ae_data.Initialize();

    // Iterate over segmented conditions
    std::vector<PointType::Pointer> points_array (TDim); // The points are stored as local coordinates, we calculate the global coordinates of this points
    PointType global_point;
    for (IndexType i_geom = 0; i_geom < rConditionsPointsSlave.size(); ++i_geom) {
        for (IndexType i_node = 0; i_node < TDim; ++i_node) {
            rSlaveGeometry.GlobalCoordinates(global_point, rConditionsPointsSlave[i_geom][i_node]);
            points_array[i_node] = Kratos::make_shared<PointType>( global_point );
        }

        DecompositionType decomp_geom( PointerVector<PointType>{points_array} );

        bool bad_shape;
        if constexpr (TDim == 2) {
            bad_shape = MortarUtilities::LengthCheck(decomp_geom, rSlaveGeometry.Length() * CheckThresholdCoefficient);
        } else { 
            bad_shape = MortarUtilities::HeronCheck(decomp_geom);
        }

        if (!bad_shape) {
            const GeometryType::IntegrationPointsArrayType& r_integration_points_slave = decomp_geom.IntegrationPoints( rIntegrationMethod );

            // Integrating the mortar operators
            PointType local_point_parent, gp_global;
            for ( IndexType point_number = 0; point_number < r_integration_points_slave.size(); ++point_number ) {
                const PointType local_point_decomp{r_integration_points_slave[point_number].Coordinates()};
                decomp_geom.GlobalCoordinates(gp_global, local_point_decomp);
                rSlaveGeometry.PointLocalCoordinates(local_point_parent, gp_global);

                rSlaveGeometry.ShapeFunctionsValues( rVariables.NSlave, local_point_parent.Coordinates() );
                rVariables.DetjSlave = decomp_geom.DeterminantOfJacobian( local_point_decomp );

                // Integrate
                const double integration_weight = AxiSymCoeff * r_integration_points_slave[point_number].Weight();

                Ae_data.CalculateAeComponents(rVariables, integration_weight);
            }
        }
    }

    return Ae_data.CalculateAe(rAe);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::ExplicitCalculateKinematics(
    const PairedCondition* pCondition,
    GeneralVariables& rVariables,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rAe,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    const GeometryPointType& rGeometryDecomp,
    const bool DualLM
    )
{
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    const auto& r_geometry = pCondition->GetParentGeometry();
    r_geometry.ShapeFunctionsValues( rVariables.NSlave, rLocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM) ? prod(rAe, rVariables.NSlave) : rVariables.NSlave;

    /* SHAPE FUNCTION DERIVATIVES */
    r_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, rLocalPointParent );

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = rGeometryDecomp.Jacobian( rVariables.jSlave, rLocalPointDecomp.Coordinates());
    rVariables.DetjSlave = rGeometryDecomp.DeterminantOfJacobian( rLocalPointDecomp );

    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "ERROR:: CONDITION ID: " << pCondition->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;

    /// MASTER CONDITION ///
    MasterShapeFunctionValue(pCondition, rVariables, rNormalMaster, rLocalPointParent);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::CalculateKinematics(
    const PairedCondition* pCondition,
    GeneralVariables& rVariables,
    const DerivativeDataType& rDerivativeData,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPointDecomp,
    const PointType& rLocalPointParent,
    const GeometryPointType& rGeometryDecomp,
    const bool DualLM
    )
{
    /// SLAVE CONDITION ///
    /* SHAPE FUNCTIONS */
    const auto& r_geometry = pCondition->GetParentGeometry();
    r_geometry.ShapeFunctionsValues( rVariables.NSlave, rLocalPointParent.Coordinates() );
    rVariables.PhiLagrangeMultipliers = (DualLM) ? prod(rDerivativeData.Ae, rVariables.NSlave) : rVariables.NSlave;

    /* SHAPE FUNCTION DERIVATIVES */
    r_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeSlave, rLocalPointParent );

    /* CALCULATE JACOBIAN AND JACOBIAN DETERMINANT */
    rVariables.jSlave = rGeometryDecomp.Jacobian( rVariables.jSlave, rLocalPointDecomp.Coordinates());
    rVariables.DetjSlave = rGeometryDecomp.DeterminantOfJacobian( rLocalPointDecomp );

    KRATOS_ERROR_IF(rVariables.DetjSlave < 0.0) << "ERROR:: CONDITION ID: " << pCondition->Id() << " INVERTED. DETJ: " << rVariables.DetjSlave << std::endl;

    /// MASTER CONDITION ///
    MasterShapeFunctionValue(pCondition, rVariables, rNormalMaster, rLocalPointParent);
}

/***********************************************************************************/
/***********************************************************************************/

template< SizeType TDim, SizeType TNumNodes, FrictionalCase TFrictional, bool TNormalVariation, SizeType TNumNodesMaster>
void MortarExplicitContributionUtilities<TDim,TNumNodes,TFrictional, TNormalVariation,TNumNodesMaster>::MasterShapeFunctionValue(
    const PairedCondition* pCondition,
    GeneralVariables& rVariables,
    const array_1d<double, 3>& rNormalMaster,
    const PointType& rLocalPoint
    )
{
    const GeometryType& r_slave_geometry = pCondition->GetParentGeometry();
    const GeometryType& r_master_geometry = pCondition->GetPairedGeometry();

    PointType projected_gp_global;
    const array_1d<double,3> gp_normal = MortarUtilities::GaussPointUnitNormal(rVariables.NSlave, r_slave_geometry);

    GeometryType::CoordinatesArrayType slave_gp_global;
    r_slave_geometry.GlobalCoordinates( slave_gp_global, rLocalPoint );
    const PointType slave_gp_global_point = PointType(slave_gp_global);
    GeometricalProjectionUtilities::FastProjectDirection( r_master_geometry, slave_gp_global_point, projected_gp_global, rNormalMaster, -gp_normal ); // The opposite direction

    GeometryType::CoordinatesArrayType projected_gp_local;

    r_master_geometry.PointLocalCoordinates( projected_gp_local, projected_gp_global.Coordinates( ) ) ;

    // SHAPE FUNCTIONS
    r_master_geometry.ShapeFunctionsValues(         rVariables.NMaster,    projected_gp_local );
    r_master_geometry.ShapeFunctionsLocalGradients( rVariables.DNDeMaster, projected_gp_local );

    // JACOBIAN
    rVariables.jMaster = r_master_geometry.Jacobian( rVariables.jMaster, projected_gp_local);
}

/***********************************************************************************/
/***********************************************************************************/

double AuxiliaryOperationsUtilities::GetAxisymmetricCoefficient(
    const PairedCondition* pCondition,
    const Vector& rNSlave
    )
{
    const double radius = AuxiliaryOperationsUtilities::CalculateRadius(pCondition, rNSlave);
    const double thickness = pCondition->GetProperties()[THICKNESS];
    return (2.0 * Globals::Pi * radius/thickness);
}

/***********************************************************************************/
/***********************************************************************************/

double AuxiliaryOperationsUtilities::CalculateRadius(
    const PairedCondition* pCondition,
    const Vector& rNSlave
    )
{
    KRATOS_TRY;

    double current_radius = 0.0;

    const auto& r_geometry = pCondition->GetParentGeometry();
    for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
        // Displacement from the reference to the current configuration
        const array_1d<double, 3 >& r_current_position = r_geometry[i_node].Coordinates();
        current_radius   += r_current_position[0] * rNSlave[i_node];
    }

    return current_radius;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

// Frictionless cases
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS, false, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS, false, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS, false, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS, false, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS, false, 3>;
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS, true, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS, true, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS, true, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS, true, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS, true, 3>;

// Frictionless components cases
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, false, 3>;
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_COMPONENTS, true, 3>;

// Frictional cases
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONAL, false, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL, false, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL, false, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL, false, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL, false, 3>;
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONAL, true, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL, true, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL, true, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL, true, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL, true, 3>;

// Frictionless penalty cases
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, false, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, false, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, false, 3>;
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONLESS_PENALTY, true, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONLESS_PENALTY, true, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONLESS_PENALTY, true, 3>;

// Frictional penalty cases
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONAL_PENALTY, false, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL_PENALTY, false, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL_PENALTY, false, 3>;
template class MortarExplicitContributionUtilities<2, 2, FrictionalCase::FRICTIONAL_PENALTY, true, 2>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
template class MortarExplicitContributionUtilities<3, 3, FrictionalCase::FRICTIONAL_PENALTY, true, 4>;
template class MortarExplicitContributionUtilities<3, 4, FrictionalCase::FRICTIONAL_PENALTY, true, 3>;

} // namespace Kratos
