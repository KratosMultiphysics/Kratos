//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//


// System includes
#include <tuple>
#include <sstream>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "custom_utilities/geometrical/opt_app_model_part_utils.h"
#include "optimization_application_variables.h"

// Include base h
#include "max_overhang_response_utils.h"

namespace Kratos
{

double MaxOverhangAngleResponseUtils::CalculateValue(
    const std::vector<ModelPart*>& rModelParts,
    const Parameters ResponseSettings)
{
    double value = 0.0;
    for (auto p_model_part : rModelParts) {
        const double local_condition_value = block_for_each<SumReduction<double>>(p_model_part->Conditions(), [&](const auto& rCondition) {
            return CalculateConditionValue(rCondition,ResponseSettings);
        });
        value += p_model_part->GetCommunicator().GetDataCommunicator().SumAll(local_condition_value);
    }

    double total_area = 0.0;
    for (auto p_model_part : rModelParts) {
        const double local_condition_area = block_for_each<SumReduction<double>>(p_model_part->Conditions(), [&](const auto& rCondition) {
            return rCondition.GetGeometry().Area();
        });
        total_area += p_model_part->GetCommunicator().GetDataCommunicator().SumAll(local_condition_area);
    }

    return value/total_area;
}

double MaxOverhangAngleResponseUtils::CalculateConditionValue(
    const Condition& rCondition,
    const Parameters ResponseSettings)
{
    array_3d print_direction = ResponseSettings["print_direction"].GetVector();
    KRATOS_ERROR_IF_NOT(norm_2(print_direction) > std::numeric_limits<double>::epsilon())
    << "MaxOverhangAngleResponseUtils:print_direction: "<<print_direction<<" has size close to zero\n";
    print_direction /= norm_2(print_direction);

    const double max_angle = ResponseSettings["max_angle"].GetDouble();
    KRATOS_ERROR_IF(!(max_angle >= 0.0 && max_angle <= 90.0))
    << "MaxOverhangAngleResponseUtils: max_angle: "<<max_angle<<" should be in range [0,90]\n";
    const double sin_max_angle = std::sin(max_angle * Globals::Pi / 180);


    const double heaviside_beta = ResponseSettings["heaviside_beta"].GetDouble();
    KRATOS_ERROR_IF(std::signbit(heaviside_beta))
    << "MaxOverhangAngleResponseUtils: heaviside_beta: "<<heaviside_beta<<" should be positive\n";

    const double penalty_factor = ResponseSettings["penalty_factor"].GetDouble();
    KRATOS_ERROR_IF(std::signbit(penalty_factor))
    << "MaxOverhangAngleResponseUtils: penalty_factor: "<<penalty_factor<<" should be positive\n";

    const array_3d local_coords = ZeroVector(3);
    const array_3d& normal = rCondition.GetGeometry().UnitNormal(local_coords);
    const double area = rCondition.GetGeometry().Area();
    const double ratio = -inner_prod(print_direction, normal)/sin_max_angle;
    double pow_val = -2.0 * heaviside_beta * (ratio-1);
    pow_val = std::clamp(pow_val, -700.0, 700.0);

    return area * (1.0/(1+std::exp(pow_val))) * std::pow(ratio,penalty_factor);
}

void MaxOverhangAngleResponseUtils::CalculateSensitivity(
    const std::vector<ModelPart*>& rEvaluatedModelParts,
    const SensitivityVariableModelPartsListMap& rSensitivityVariableModelPartInfo,
    const Parameters ResponseSettings)
{
    KRATOS_TRY

    // calculate sensitivities for each and every model part w.r.t. their sensitivity variables list
    for (const auto& it : rSensitivityVariableModelPartInfo) {
        std::visit([&](auto&& r_variable) {
            const auto& r_sensitivity_model_parts = OptAppModelPartUtils::GetModelPartsWithCommonReferenceEntities(
                it.second, rEvaluatedModelParts, false, true, false, false, 0);

            // reset nodal common interface values
            for (auto p_sensitivity_model_part : r_sensitivity_model_parts) {
                if (*r_variable == SHAPE_SENSITIVITY) {
                    VariableUtils().SetNonHistoricalVariablesToZero(p_sensitivity_model_part->Nodes(), SHAPE_SENSITIVITY);
                }
            }

            // now compute sensitivities on the variables
            for (auto p_sensitivity_model_part : r_sensitivity_model_parts) {
                if (*r_variable == SHAPE_SENSITIVITY) {
                    CalculateFiniteDifferenceShapeSensitivity(*p_sensitivity_model_part, ResponseSettings, SHAPE_SENSITIVITY);
                } else {
                    KRATOS_ERROR
                        << "Unsupported sensitivity w.r.t. " << r_variable->Name()
                        << " requested. Followings are supported sensitivity variables:"
                        << "\n\t" << SHAPE_SENSITIVITY.Name();
                }
            }
        }, it.first);
    }

    KRATOS_CATCH("");
}

void MaxOverhangAngleResponseUtils::CalculateFiniteDifferenceShapeSensitivity(
    ModelPart& rModelPart,
    const Parameters ResponseSettings,
    const Variable<array_3d>& rOutputSensitivityVariable)
{
    KRATOS_TRY

    using tls_condition = Condition::Pointer;

    const int max_id = block_for_each<MaxReduction<int>>(rModelPart.GetRootModelPart().Nodes(), [](const auto& rNode) {
        return rNode.Id();
    });

    std::vector<std::string> model_part_names;

    block_for_each(rModelPart.Conditions(), tls_condition(), [&](auto& rCondition, tls_condition& rTLS) {
        CalculateConditionFiniteDifferenceShapeSensitivity(
            rCondition,
            rTLS,
            rModelPart,
            model_part_names,
            ResponseSettings,
            max_id + ParallelUtilities::GetNumThreads() * 1000, // no element or condition is not suppose to have 1000 nodes per element or condition. This is done to avoid re-calculating the max id for conditions
            rOutputSensitivityVariable);
    });

    // now clear the temp model part
    VariableUtils().SetFlag(TO_ERASE, false, rModelPart.Nodes());
    for (const auto& r_model_part_name : model_part_names) {
        auto& tls_model_part = rModelPart.GetSubModelPart(r_model_part_name);
        for (auto& r_node : tls_model_part.Nodes()) {
            r_node.Set(TO_ERASE, true);
        }
    }

    // remove temporary nodes
    rModelPart.RemoveNodesFromAllLevels(TO_ERASE);

    // remove temporary model parts
    for (const auto& r_model_part_name : model_part_names) {
        rModelPart.RemoveSubModelPart(r_model_part_name);
    }

    // Assemble nodal result
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rOutputSensitivityVariable);

    KRATOS_CATCH("");
}

void MaxOverhangAngleResponseUtils::CalculateConditionFiniteDifferenceShapeSensitivity(
    Condition& rCondition,
    Condition::Pointer& pThreadLocalCondition,
    ModelPart& rModelPart,
    std::vector<std::string>& rModelPartNames,
    const Parameters ResponseSettings,
    const IndexType MaxNodeId,
    const Variable<array_3d>& rOutputSensitivityVariable)
{
    KRATOS_TRY

    const auto& r_process_info = rModelPart.GetProcessInfo();
    auto& r_geometry = rCondition.GetGeometry();
    const auto domain_size = r_geometry.WorkingSpaceDimension();
    const double pert_size =  ResponseSettings["gradient_settings"]["step_size"].GetDouble();

    // calculate the reference value
    const double ref_val = CalculateConditionValue(rCondition,ResponseSettings);

    // initialize dummy element for parallelized perturbation based sensitivity calculation
    // in each thread seperately
    if (!pThreadLocalCondition) {

        std::stringstream name;
        const int thread_id = OpenMPUtils::ThisThread();
        name << rModelPart.Name() << "_temp_" << thread_id;
        name << "_condition";


        GeometryType::PointsArrayType nodes;
        {
            KRATOS_CRITICAL_SECTION

            rModelPartNames.push_back(name.str());
            auto& tls_model_part = rModelPart.CreateSubModelPart(name.str());

            for (IndexType i = 0; i < r_geometry.size(); ++i) {
                nodes.push_back(tls_model_part.CreateNewNode((thread_id + 1) * MaxNodeId + i + 1, r_geometry[i]));
            }
        }

        pThreadLocalCondition = rCondition.Create(1, nodes, rCondition.pGetProperties());
        pThreadLocalCondition->Initialize(r_process_info);
        pThreadLocalCondition->InitializeSolutionStep(r_process_info);
    }

    *pThreadLocalCondition = static_cast<const Condition&>(rCondition);

    // TODO: For some reason, assignment operator of the element
    //       assigns current position to initial position of the lhs element
    //       hence these are corrected manually. But need to be checked
    //       in the element/node/geometry assignment operators.
    pThreadLocalCondition->GetGeometry().SetData(r_geometry.GetData());
    pThreadLocalCondition->SetFlags(rCondition.GetFlags());
    for (IndexType i = 0; i < r_geometry.size(); ++i) {
        auto& r_node = pThreadLocalCondition->GetGeometry()[i];
        const auto& r_orig_node = r_geometry[i];
        r_node = r_orig_node;
    }

    // now calculate perturbed
    for (IndexType i = 0; i < r_geometry.size(); ++i) {
        auto& r_orig_node_sensitivity = r_geometry[i].GetValue(rOutputSensitivityVariable);

        auto& r_node = pThreadLocalCondition->GetGeometry()[i];
        auto& r_coordinates = r_node.Coordinates();
        auto& r_initial_coordintes = r_node.GetInitialPosition();

        r_initial_coordintes[0] += pert_size;
        r_coordinates[0] += pert_size;
        double pert_val = CalculateConditionValue(*pThreadLocalCondition,ResponseSettings);
        r_initial_coordintes[0] -= pert_size;
        r_coordinates[0] -= pert_size;
        AtomicAdd<double>(r_orig_node_sensitivity[0], (pert_val - ref_val) / pert_size);

        r_initial_coordintes[1] += pert_size;
        r_coordinates[1] += pert_size;
        pert_val = CalculateConditionValue(*pThreadLocalCondition,ResponseSettings);
        r_initial_coordintes[1] -= pert_size;
        r_coordinates[1] -= pert_size;
        AtomicAdd<double>(r_orig_node_sensitivity[1], (pert_val - ref_val) / pert_size);


        if (domain_size == 3) {
            r_initial_coordintes[2] += pert_size;
            r_coordinates[2] += pert_size;
            pert_val = CalculateConditionValue(*pThreadLocalCondition,ResponseSettings);
            r_initial_coordintes[2] -= pert_size;
            r_coordinates[2] -= pert_size;
            AtomicAdd<double>(r_orig_node_sensitivity[2], (pert_val - ref_val) / pert_size);
        }
    }

    KRATOS_CATCH("");
}

}