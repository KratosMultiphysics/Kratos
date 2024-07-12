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
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/contact_utilities.h"
#include "custom_processes/advanced_contact_search_process.h"
#include "contact_structural_mechanics_application_variables.h"

namespace Kratos
{
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::AdvancedContactSearchProcess(
    ModelPart & rMainModelPart,
    Parameters ThisParameters,
    Properties::Pointer pPairedProperties
    ) : BaseType(rMainModelPart, ThisParameters, pPairedProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CheckPairing(
    ModelPart& rComputingModelPart,
    IndexType& rConditionId
    )
{
    // Getting the corresponding submodelparts
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = this->IsNotMultipleSearchs() ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub" + BaseType::mThisParameters["id_name"].GetString());
    ModelPart& r_master_model_part = r_sub_contact_model_part.GetSubModelPart("MasterSubModelPart" + BaseType::mThisParameters["id_name"].GetString());
    ModelPart& r_slave_model_part = r_sub_contact_model_part.GetSubModelPart("SlaveSubModelPart" + BaseType::mThisParameters["id_name"].GetString());

    // We compute the maximal nodal h and some auxiliary values  // TODO: Think about this criteria
    const double distance_threshold = std::max(ContactUtilities::CalculateMeanNodalH(r_slave_model_part), ContactUtilities::CalculateMeanNodalH(r_master_model_part));
//     const double distance_threshold = 2.0/3.0 * ContactUtilities::CalculateMeanNodalH(mrMainModelPart);
//     const double distance_threshold = ContactUtilities::CalculateMeanNodalH(mrMainModelPart);
//     const double distance_threshold = ContactUtilities::CalculateMaxNodalH(mrMainModelPart)();
//     const double distance_threshold = mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR] * ContactUtilities::CalculateMaxNodalH(mrMainModelPart)();

    // Updating the distance distance threshold
    BaseType::mrMainModelPart.GetProcessInfo().SetValue(DISTANCE_THRESHOLD, distance_threshold);

    // Calling base class
    BaseType::CheckPairing(rComputingModelPart, rConditionId);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::ComputeActiveInactiveNodes()
{
    // We get the process info
    const ProcessInfo& r_process_info = BaseType::mrMainModelPart.GetProcessInfo();

    // The penalty value
    const double common_epsilon = r_process_info[INITIAL_PENALTY];

    // Some auxiliary values
    const double active_check_factor = BaseType::mrMainModelPart.GetProcessInfo()[ACTIVE_CHECK_FACTOR];
    const double distance_threshold = BaseType::mrMainModelPart.GetProcessInfo()[DISTANCE_THRESHOLD];
    double reference_auxiliary_length = distance_threshold * active_check_factor;
    double auxiliary_length = reference_auxiliary_length;

    // Compute linear regression
    const bool predict_correct_lagrange_multiplier = BaseType::mThisParameters["predict_correct_lagrange_multiplier"].GetBool();
    double a = 0.0;
    double b = 0.0;
    if (predict_correct_lagrange_multiplier) {
        ComputeLinearRegressionGapPressure(a, b);
    }

    // Iterate over the nodes
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = this->IsNotMultipleSearchs() ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub" + BaseType::mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();

    // If we do an static check
    const bool static_check_movement = BaseType::mThisParameters["static_check_movement"].GetBool();
    const bool consider_gap_threshold = BaseType::mThisParameters["consider_gap_threshold"].GetBool();

    // We compute now the normal gap and set the nodes under certain threshold as active
    struct auxiliary {bool auxiliary_check = false; bool has_weighted_gap = false; double auxiliary_length = 0.0;};
    block_for_each(r_nodes_array, auxiliary(), [this, &reference_auxiliary_length, &consider_gap_threshold, &static_check_movement, &predict_correct_lagrange_multiplier, &a, &b, &common_epsilon](Node& rNode, auxiliary& aux) {
        if (rNode.Is(SLAVE) == this->IsNotInvertedSearch()) {
            aux.auxiliary_check = false;
            aux.auxiliary_length = reference_auxiliary_length;
            aux.has_weighted_gap = rNode.SolutionStepsDataHas(WEIGHTED_GAP);
            const double weighted_gap = aux.has_weighted_gap ? rNode.FastGetSolutionStepValue(WEIGHTED_GAP) : 0.0;
            if (rNode.Is(ACTIVE) && aux.has_weighted_gap) {
                const double nodal_area = rNode.Has(NODAL_AREA) ? rNode.GetValue(NODAL_AREA) : 1.0;
                aux.auxiliary_check = (weighted_gap/nodal_area < aux.auxiliary_length) ? true : false;
            }
            // We do the static check
            if (static_check_movement && aux.has_weighted_gap) {
                const double movement_check = weighted_gap - rNode.FastGetSolutionStepValue(WEIGHTED_GAP, 1);
                if (movement_check > -(std::abs(weighted_gap) * 1.0e-3)) {
                    // If we consider the threshold
                    if (consider_gap_threshold) {
                        if (aux.auxiliary_length < GapThreshold) {
                            aux.auxiliary_length = GapThreshold;
                        }
                    }
                }
            }
            if ((rNode.GetValue(NORMAL_GAP) < aux.auxiliary_length) || aux.auxiliary_check) {
                if (predict_correct_lagrange_multiplier) {
                    SetActiveNodeWithRegression(rNode, a, b);
                } else { // We just mark it
                    this->SetActiveNode(rNode, common_epsilon);
                }
            } else {
            #ifdef KRATOS_DEBUG
                KRATOS_WARNING_IF("BaseContactSearch", rNode.Is(ACTIVE)) << "WARNING: A node that used to be active is not active anymore. Check that. Node ID: " << rNode.Id() << std::endl;
            #endif
                this->SetInactiveNode(rNode);
            }
        }
    });

    // In case of debug mode
    if (BaseType::mThisParameters["debug_mode"].GetBool()) {
        std::filebuf debug_buffer;
        debug_buffer.open("gap_active_nodes_debug_" + r_sub_contact_model_part.Name() + "_step=" + std::to_string( r_sub_contact_model_part.GetProcessInfo()[STEP]) + ".out",std::ios::out);
        std::ostream os(&debug_buffer);
        for (const auto& r_node : r_nodes_array) {
            os << "Node " << r_node.Id() << "\tNORMAL GAP:" << "\t" << r_node.GetValue(NORMAL_GAP) << "\tTHRESHOLD:\t" << auxiliary_length;
            os << "\tNORMAL:" << "\t" << r_node.FastGetSolutionStepValue(NORMAL_X) << "\t" << r_node.FastGetSolutionStepValue(NORMAL_Y) << "\t" << r_node.FastGetSolutionStepValue(NORMAL_Z);
//             os << "\tDISPLACEMENT:" << "\t" << r_node.FastGetSolutionStepValue(DISPLACEMENT_X) << "\t" << r_node.FastGetSolutionStepValue(DISPLACEMENT_Y) << "\t" << r_node.FastGetSolutionStepValue(DISPLACEMENT_Z);
            os << "\tAUXILIAR_COORDINATES:" << "\t" << r_node.GetValue(AUXILIAR_COORDINATES_X) << "\t" << r_node.GetValue(AUXILIAR_COORDINATES_Y) << "\t" << r_node.GetValue(AUXILIAR_COORDINATES_Z);
            os << "\tWEIGHTED_GAP: " << r_node.FastGetSolutionStepValue(WEIGHTED_GAP);
            os << "\tNODAL_AREA: " << r_node.GetValue(NODAL_AREA);
            os << "\tACTIVE: " << r_node.Is(ACTIVE) <<"\n";
        }
        debug_buffer.close();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::ComputeLinearRegressionGapPressure(
    double& a,
    double& b
    )
{
    // Initialize the n
    SizeType n = 0;

    // Initialize the values
    struct AuxValues {
        double xi = 0.0;
        double yi = 0.0;
    };
    double sum_x, sum_xsq, sum_y, sum_xy;

    sum_x = 0.0;
    sum_xsq = 0.0;
    sum_y = 0.0;
    sum_xy = 0.0;

    // Iterate over the nodes
    ModelPart& r_contact_model_part = BaseType::mrMainModelPart.GetSubModelPart("Contact");
    ModelPart& r_sub_contact_model_part = this->IsNotMultipleSearchs() ? r_contact_model_part : r_contact_model_part.GetSubModelPart("ContactSub" + BaseType::mThisParameters["id_name"].GetString());
    NodesArrayType& r_nodes_array = r_sub_contact_model_part.Nodes();

    // We compute now the normal gap and set the nodes under certain threshold as active
    using FiveReduction = CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<double>, SumReduction<SizeType>>;
    std::tie(sum_x,sum_xsq,sum_y,sum_xy,n) = block_for_each<FiveReduction>(r_nodes_array, AuxValues(), [&](Node& rNode, AuxValues& aux_values) {
        // We get the condition geometry
        if (rNode.Is(ACTIVE)) {
            aux_values.xi = rNode.FastGetSolutionStepValue(WEIGHTED_GAP);
            if (BaseType::mTypeSolution == BaseType::TypeSolution::NormalContactStress) {
                aux_values.yi = rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
            } else {
                const array_1d<double, 3>& lm = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                const array_1d<double, 3>& normal = rNode.FastGetSolutionStepValue(NORMAL);
                aux_values.yi = inner_prod(lm, normal);
            }
            return std::make_tuple(aux_values.xi, std::pow(aux_values.xi, 2), aux_values.yi, aux_values.xi * aux_values.yi, 1);
        }
        return std::make_tuple(0.0,0.0,0.0,0.0,0);
    });

    const double size = static_cast<double>(n);
    const double denom = size * sum_xsq - std::pow(sum_x, 2);
    a = (sum_y * sum_xsq - sum_x * sum_xy) / denom;
    b = (size * sum_xy - sum_x * sum_y) / denom;

    // In case of debug mode
    if (BaseType::mThisParameters["debug_mode"].GetBool()) {
        std::filebuf debug_buffer;
        debug_buffer.open("linear_regression_gap_" + r_sub_contact_model_part.Name() + "_step=" + std::to_string( r_sub_contact_model_part.GetProcessInfo()[STEP]) + ".out",std::ios::out);
        std::ostream os(&debug_buffer);
        double contact_pressure = 0.0;
        for (const auto& r_node : r_nodes_array) {
            if (r_node.Is(ACTIVE)) {
                os << "Node " << r_node.Id() << "\tWEIGHTED GAP:" << "\t" << r_node.GetValue(WEIGHTED_GAP);
                os << "\tNORMAL GAP:" << "\t" << r_node.GetValue(NORMAL_GAP);
                os << "\tNODAL_AREA: " << r_node.GetValue(NODAL_AREA);
                if (BaseType::mTypeSolution == BaseType::TypeSolution::NormalContactStress) {
                    contact_pressure = r_node.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                } else {
                    const array_1d<double, 3>& lm = r_node.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
                    const array_1d<double, 3>& normal = r_node.FastGetSolutionStepValue(NORMAL);
                    contact_pressure = inner_prod(lm, normal);
                }
                os << "\tNORMAL PRESSURE:" << "\t" << contact_pressure<< "\n";
            }
        }
        os << "\n\nREGRESSION VALUES:\ta = " << a << "\tb = " << b;
        debug_buffer.close();
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::SetActiveNodeWithRegression(
    Node& rNode,
    const double a,
    const double b
    )
{
    if (rNode.Is(ACTIVE) ) {
        switch(BaseType::mTypeSolution) {
            case BaseType::TypeSolution::VectorLagrangeMultiplier :
                if (BaseType::mrMainModelPart.Is(SLIP))
                    CorrectALMFrictionalMortarLM(rNode, a, b);
                else if (BaseType::mrMainModelPart.Is(CONTACT))
                    CorrectALMFrictionlessComponentsMortarLM(rNode, a, b);
                else
                    CorrectComponentsMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::ScalarLagrangeMultiplier :
                CorrectScalarMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::NormalContactStress :
                CorrectALMFrictionlessMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::FrictionlessPenaltyMethod :
                break;
            case BaseType::TypeSolution::FrictionalPenaltyMethod :
                break;
            case BaseType::TypeSolution::OtherFrictionless :
                break;
            case BaseType::TypeSolution::OtherFrictional :
                break;
        }
    } else {
        rNode.Set(ACTIVE, true);
        switch(BaseType::mTypeSolution) {
            case BaseType::TypeSolution::VectorLagrangeMultiplier :
                if (BaseType::mrMainModelPart.Is(SLIP))
                    PredictALMFrictionalMortarLM(rNode, a, b);
                else if (BaseType::mrMainModelPart.Is(CONTACT))
                    PredictALMFrictionlessComponentsMortarLM(rNode, a, b);
                else
                    PredictComponentsMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::ScalarLagrangeMultiplier :
                PredictScalarMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::NormalContactStress :
                PredictALMFrictionlessMortarLM(rNode, a, b);
                break;
            case BaseType::TypeSolution::FrictionlessPenaltyMethod :
                break;
            case BaseType::TypeSolution::FrictionalPenaltyMethod :
                break;
            case BaseType::TypeSolution::OtherFrictionless :
                break;
            case BaseType::TypeSolution::OtherFrictional :
                break;
        }
    }

    // If finally is active we set as visited
    if (rNode.Is(ACTIVE)) rNode.Set(MARKER);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CorrectScalarMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    // TODO: Add correction
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CorrectComponentsMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    // TODO: Add correction
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CorrectALMFrictionlessMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
//     const double old_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP, 1);
    const double current_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

    double& r_current_contact_stress = rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);

    // Apply linear regression
    const double aux_press = a + current_weighted_gap * b;
    r_current_contact_stress = (aux_press < 0.0) ? aux_press : 0.0;

//     const bool old_penetration = (old_weighted_gap < 0.0) ? true : false;
//     const bool current_penetration = (current_weighted_gap < 0.0) ? true : false;
//
//     // If both are positive or negative we just interpolate the current values
//     if (old_penetration == current_penetration) {
//         if (old_penetration) { // Penetrating
//             if (std::abs(old_weighted_gap) > ZeroTolerance)
//                 r_current_contact_stress *= current_weighted_gap/old_weighted_gap;
//         } else { // Not penetrating
//             if (old_weighted_gap > ZeroTolerance)
//                 r_current_contact_stress /= current_weighted_gap/old_weighted_gap;
//         }
//     } else if (old_penetration && !current_penetration) { // We had penenetration and we don't have it anymore
//         const double gap_variation = (current_weighted_gap + old_weighted_gap);
//         if (std::abs(gap_variation) > ZeroTolerance) {
//             if (gap_variation > 0.0)
//                 r_current_contact_stress *= - old_weighted_gap/gap_variation;
//             else
//                 r_current_contact_stress *= - current_weighted_gap/gap_variation;
//         }
//     } else { // We have penetration and we haven't before
//         const double gap_variation = (old_weighted_gap - current_weighted_gap);
//         if (std::abs(gap_variation) > ZeroTolerance) {
//             if (gap_variation > 0.0)
//                 r_current_contact_stress /= - current_weighted_gap/gap_variation;
//             else
//                 r_current_contact_stress /= - old_weighted_gap/gap_variation;
//         }
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CorrectALMFrictionlessComponentsMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
//     const double old_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP, 1);
    const double current_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP);

    array_1d<double, 3>& r_current_lm = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
    const array_1d<double, 3>& r_current_normal = rNode.FastGetSolutionStepValue(NORMAL);

    // Apply linear regression
    const double aux_press = a + current_weighted_gap * b;
    const array_1d<double, 3> zero_vector = ZeroVector(3);
    noalias(r_current_lm) = (aux_press < 0.0) ? aux_press * r_current_normal : zero_vector;

//     const bool old_penetration = (old_weighted_gap < 0.0) ? true : false;
//     const bool current_penetration = (current_weighted_gap < 0.0) ? true : false;
//
//     // If both are positive or negative we just interpolate the current values
//     if (old_penetration == current_penetration) {
//         if (old_penetration) { // Penetrating
//             if (std::abs(old_weighted_gap) > ZeroTolerance)
//                 r_current_lm *= current_weighted_gap/old_weighted_gap;
//         } else { // Not penetrating
//             if (old_weighted_gap > ZeroTolerance)
//                 r_current_lm /= current_weighted_gap/old_weighted_gap;
//         }
//     } else if (old_penetration && !current_penetration) { // We had penenetration and we don't have it anymore
//         const double gap_variation = (current_weighted_gap + old_weighted_gap);
//         if (std::abs(gap_variation) > ZeroTolerance) {
//             if (gap_variation > 0.0)
//                 r_current_lm *= - old_weighted_gap/gap_variation;
//             else
//                 r_current_lm *= - current_weighted_gap/gap_variation;
//         }
//     } else { // We have penetration and we haven't before
//         const double gap_variation = (old_weighted_gap - current_weighted_gap);
//         if (std::abs(gap_variation) > ZeroTolerance) {
//             if (gap_variation > 0.0)
//                 r_current_lm /= - current_weighted_gap/gap_variation;
//             else
//                 r_current_lm /= - old_weighted_gap/gap_variation;
//         }
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::CorrectALMFrictionalMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    if (norm_2(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
        if (rNode.GetValue(FRICTION_COEFFICIENT) < ZeroTolerance || this->IsPureSlip()) {
            rNode.Set(SLIP, true);
        } else {
            rNode.Set(SLIP, false);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::PredictScalarMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    // TODO: Add correction
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::PredictComponentsMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    // TODO: Add correction
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::PredictALMFrictionlessMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
//     // The penalty to be use (TODO: think about use the nodal penalty)
//     ProcessInfo& current_process_info = mrMainModelPart.GetProcessInfo(); // TODO: Avoid call the process info each time
//     const double initial_penalty = current_process_info[INITIAL_PENALTY];
//     const double distance_threshold = current_process_info[DISTANCE_THRESHOLD];
//
//     // If we have penetration
    const double current_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP);
//     const double nodal_area = rNode.GetValue(NODAL_AREA);
//     const bool current_penetration = (current_weighted_gap < 0.0) ? true : false;
//
    double& r_current_contact_stress = rNode.FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);

    // Apply linear regression
    const double aux_press = a + current_weighted_gap * b;
    r_current_contact_stress = (aux_press < 0.0) ? aux_press : 0.0;
//
//     // We have penetration so we just basically approximate the solution with the traditional
//     if (current_penetration) {
//         r_current_contact_stress = initial_penalty * current_weighted_gap;
//     } else { // We don't have penetration, we do a simpler approach
//         const double relative_gap = (current_weighted_gap - distance_threshold * nodal_area);
//         r_current_contact_stress = (relative_gap < 0.0) ? initial_penalty * relative_gap : 0.0;
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::PredictALMFrictionlessComponentsMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
//     // The penalty to be use (TODO: think about use the nodal penalty)
//     ProcessInfo& current_process_info = mrMainModelPart.GetProcessInfo(); // TODO: Avoid call the process info each time
//     const double initial_penalty = current_process_info[INITIAL_PENALTY];
//     const double distance_threshold = current_process_info[DISTANCE_THRESHOLD];
//
//     // If we have penetration
    const double current_weighted_gap = rNode.FastGetSolutionStepValue(WEIGHTED_GAP);
//     const double nodal_area = rNode.GetValue(NODAL_AREA);
//     const bool current_penetration = (current_weighted_gap < 0.0) ? true : false;
//
    array_1d<double, 3>& r_current_lm = rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER);
    const array_1d<double, 3>& r_current_normal = rNode.FastGetSolutionStepValue(NORMAL);

    // Apply linear regression
    const double aux_press = a + current_weighted_gap * b;
    const array_1d<double, 3> zero_vector = ZeroVector(3);
    noalias(r_current_lm) = (aux_press < 0.0) ? aux_press * r_current_normal : zero_vector;
//
//     // We have penetration so we just basically approximate the solution with the traditional
//     if (current_penetration) {
//         noalias(r_current_lm) = r_current_normal * initial_penalty * current_weighted_gap;
//     } else { // We don't have penetration, we do a simpler approach
//         const double relative_gap = (current_weighted_gap - distance_threshold * nodal_area);
//         noalias(r_current_lm) = (relative_gap < 0.0) ? r_current_normal * initial_penalty * relative_gap : zero_vector;
//     }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
void AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>::PredictALMFrictionalMortarLM(
    Node& rNode,
    const double a,
    const double b
    )
{
    if (norm_2(rNode.FastGetSolutionStepValue(VECTOR_LAGRANGE_MULTIPLIER)) < ZeroTolerance) {
        if (rNode.GetValue(FRICTION_COEFFICIENT) < ZeroTolerance || this->IsPureSlip()) {
            rNode.Set(SLIP, true);
        } else {
            rNode.Set(SLIP, false);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class AdvancedContactSearchProcess<2, 2>;
template class AdvancedContactSearchProcess<3, 3>;
template class AdvancedContactSearchProcess<3, 4>;
template class AdvancedContactSearchProcess<3, 3, 4>;
template class AdvancedContactSearchProcess<3, 4, 3>;

}  // namespace Kratos.
