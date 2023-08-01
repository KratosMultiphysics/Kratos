//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
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
#include "utilities/string_utilities.h"
#include "utilities/model_part_operation_utilities.h"

// Application includes
#include "optimization_application_variables.h"

// Include base h
#include "linear_strain_energy_response_utils.h"

namespace Kratos
{

template<class TEntityType>
double LinearStrainEnergyResponseUtils::CalculateEntityStrainEnergy(
    TEntityType& rEntity,
    Matrix& rLHS,
    Vector& rRHS,
    Vector& rX,
    const ProcessInfo& rProcessInfo)
{
    if (rEntity.IsActive()) {
        rEntity.GetValuesVector(rX);
        rEntity.CalculateLocalSystem(rLHS, rRHS, rProcessInfo);

        // Compute strain energy
        return 0.5 * inner_prod(rX, prod(rLHS, rX));
    } else {
        return 0.0;
    }
}

double LinearStrainEnergyResponseUtils::CalculateValue(ModelPart& rEvaluatedModelPart)
{
    KRATOS_TRY

    using tls_type = std::tuple<Matrix, Vector, Vector>;

    const double local_element_value = block_for_each<SumReduction<double>>(rEvaluatedModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        Matrix& r_lhs = std::get<0>(rTLS);
        Vector& r_rhs = std::get<1>(rTLS);
        Vector& r_u = std::get<2>(rTLS);
        return CalculateEntityStrainEnergy(rElement, r_lhs, r_rhs, r_u, rEvaluatedModelPart.GetProcessInfo());
    });

    const double local_condition_value = block_for_each<SumReduction<double>>(rEvaluatedModelPart.Conditions(), tls_type(), [&](auto& rCondition, tls_type& rTLS) {
        Matrix& r_lhs = std::get<0>(rTLS);
        Vector& r_rhs = std::get<1>(rTLS);
        Vector& r_u = std::get<2>(rTLS);
        return CalculateEntityStrainEnergy(rCondition, r_lhs, r_rhs, r_u, rEvaluatedModelPart.GetProcessInfo());
    });

    return rEvaluatedModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_element_value + local_condition_value);

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtils::CalculateGradient(
    const std::vector<GradientFieldVariableTypes>& rListOfGradientVariables,
    const std::vector<ModelPart*>& rListOfGradientRequiredModelParts,
    const std::vector<ModelPart*>& rListOfGradientComputedModelParts,
    const double PerturbationSize)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rListOfGradientVariables.size() !=
                    rListOfGradientRequiredModelParts.size())
        << "Number of gradient variables and required model parts mismatch.";
    KRATOS_ERROR_IF(rListOfGradientVariables.size() !=
                    rListOfGradientComputedModelParts.size())
        << "Number of gradient variables and computed model parts mismatch.";

    for (IndexType i = 0; i < rListOfGradientVariables.size(); ++i) {
        std::visit([&](auto p_variable) {
            if (*p_variable == YOUNG_MODULUS) {
                block_for_each(rListOfGradientRequiredModelParts[i]->Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(YOUNG_MODULUS_SENSITIVITY, 0.0); });
                CalculateStrainEnergyLinearlyDependentPropertyGradient(*rListOfGradientComputedModelParts[i], YOUNG_MODULUS, YOUNG_MODULUS_SENSITIVITY);
            } else if (*p_variable == THICKNESS) {
                block_for_each(rListOfGradientRequiredModelParts[i]->Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(THICKNESS_SENSITIVITY, 0.0); });
                CalculateStrainEnergyLinearlyDependentPropertyGradient(*rListOfGradientComputedModelParts[i], THICKNESS, THICKNESS_SENSITIVITY);
            } else if (*p_variable == POISSON_RATIO) {
                block_for_each(rListOfGradientRequiredModelParts[i]->Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(POISSON_RATIO_SENSITIVITY, 0.0); });
                CalculateStrainEnergySemiAnalyticPropertyGradient(*rListOfGradientComputedModelParts[i], PerturbationSize, POISSON_RATIO, POISSON_RATIO_SENSITIVITY);
            } else if (*p_variable == SHAPE) {
                VariableUtils().SetNonHistoricalVariableToZero(SHAPE_SENSITIVITY, rListOfGradientRequiredModelParts[i]->Nodes());
                CalculateStrainEnergySemiAnalyticShapeGradient(*rListOfGradientComputedModelParts[i], PerturbationSize, SHAPE_SENSITIVITY);
            } else {
                KRATOS_ERROR
                    << "Unsupported sensitivity w.r.t. " << p_variable->Name()
                    << " requested. Followings are supported sensitivity variables:"
                    << "\n\t" << YOUNG_MODULUS.Name()
                    << "\n\t" << THICKNESS.Name()
                    << "\n\t" << POISSON_RATIO.Name()
                    << "\n\t" << SHAPE.Name();
            }
        }, rListOfGradientVariables[i]);
    }

    KRATOS_CATCH("");
}

template<class TEntityType>
void LinearStrainEnergyResponseUtils::CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
    TEntityType& rEntity,
    Vector& rX,
    Vector& rRefRHS,
    Vector& rPerturbedRHS,
    typename TEntityType::Pointer& pThreadLocalEntity,
    ModelPart& rModelPart,
    std::vector<std::string>& rModelPartNames,
    const double Delta,
    const IndexType MaxNodeId,
    const Variable<array_1d<double, 3>>& rOutputGradientVariable)
{
    KRATOS_TRY

    if (rEntity.IsActive()) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        auto& r_geometry = rEntity.GetGeometry();
        const auto domain_size = r_geometry.WorkingSpaceDimension();

        rEntity.GetValuesVector(rX);

        // calculate the reference value
        rEntity.CalculateRightHandSide(rRefRHS, r_process_info);

        // initialize dummy element for parallelized perturbation based sensitivity calculation
        // in each thread seperately
        if (!pThreadLocalEntity) {

            std::stringstream name;
            const int thread_id = OpenMPUtils::ThisThread();
            name << rModelPart.Name() << "_temp_" << thread_id;
            if constexpr(std::is_same_v<TEntityType, ModelPart::ConditionType>) {
                name << "_condition";
            } else if constexpr(std::is_same_v<TEntityType, ModelPart::ElementType>) {
                name << "_element";
            }
            GeometryType::PointsArrayType nodes;
            {
                KRATOS_CRITICAL_SECTION

                rModelPartNames.push_back(name.str());
                auto& tls_model_part = rModelPart.CreateSubModelPart(name.str());

                for (IndexType i = 0; i < r_geometry.size(); ++i) {
                    nodes.push_back(tls_model_part.CreateNewNode((thread_id + 1) * MaxNodeId + i + 1, r_geometry[i]));
                }
            }

            pThreadLocalEntity = rEntity.Create(1, nodes, rEntity.pGetProperties());
            pThreadLocalEntity->Initialize(r_process_info);
            pThreadLocalEntity->InitializeSolutionStep(r_process_info);
        }

        *pThreadLocalEntity = static_cast<const TEntityType&>(rEntity);

        // TODO: For some reason, assignment operator of the element
        //       assigns current position to initial position of the lhs element
        //       hence these are corrected manually. But need to be checked
        //       in the element/node/geometry assignment operators.
        pThreadLocalEntity->GetGeometry().SetData(r_geometry.GetData());
        pThreadLocalEntity->SetFlags(rEntity.GetFlags());
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            auto& r_node = pThreadLocalEntity->GetGeometry()[i];
            const auto& r_orig_node = r_geometry[i];
            r_node = r_orig_node;
        }

        // now calculate perturbed
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            auto& r_orig_node_sensitivity = r_geometry[i].GetValue(rOutputGradientVariable);

            auto& r_node = pThreadLocalEntity->GetGeometry()[i];
            auto& r_coordinates = r_node.Coordinates();
            auto& r_initial_coordintes = r_node.GetInitialPosition();

            r_initial_coordintes[0] += Delta;
            r_coordinates[0] += Delta;
            pThreadLocalEntity->CalculateRightHandSide(rPerturbedRHS, r_process_info);
            r_initial_coordintes[0] -= Delta;
            r_coordinates[0] -= Delta;
            AtomicAdd<double>(r_orig_node_sensitivity[0], 0.5 * inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);

            r_initial_coordintes[1] += Delta;
            r_coordinates[1] += Delta;
            pThreadLocalEntity->CalculateRightHandSide(rPerturbedRHS, r_process_info);
            r_initial_coordintes[1] -= Delta;
            r_coordinates[1] -= Delta;
            AtomicAdd<double>(r_orig_node_sensitivity[1], 0.5 * inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);

            if (domain_size == 3) {
                r_initial_coordintes[2] += Delta;
                r_coordinates[2] += Delta;
                pThreadLocalEntity->CalculateRightHandSide(rPerturbedRHS, r_process_info);
                r_initial_coordintes[2] -= Delta;
                r_coordinates[2] -= Delta;
                AtomicAdd<double>(r_orig_node_sensitivity[2], 0.5 * inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);
            }
        }
    }

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtils::CalculateStrainEnergySemiAnalyticShapeGradient(
    ModelPart& rModelPart,
    const double Delta,
    const Variable<array_1d<double, 3>>& rOutputGradientVariable)
{
    KRATOS_TRY

    using tls_element_type = std::tuple<Vector, Vector, Vector, Element::Pointer>;

    using tls_condition_type = std::tuple<Vector, Vector, Vector, Condition::Pointer>;

    const int max_id = block_for_each<MaxReduction<int>>(rModelPart.GetRootModelPart().Nodes(), [](const auto& rNode) {
        return rNode.Id();
    });

    std::vector<std::string> model_part_names;
    VariableUtils().SetNonHistoricalVariableToZero(rOutputGradientVariable, rModelPart.Nodes());

    block_for_each(rModelPart.Elements(), tls_element_type(), [&](auto& rElement, tls_element_type& rTLS) {
        CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
            rElement,
            std::get<0>(rTLS),
            std::get<1>(rTLS),
            std::get<2>(rTLS),
            std::get<3>(rTLS),
            rModelPart,
            model_part_names,
            Delta,
            max_id,
            rOutputGradientVariable);
    });

    block_for_each(rModelPart.Conditions(), tls_condition_type(), [&](auto& rCondition, tls_condition_type& rTLS) {
        CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
            rCondition,
            std::get<0>(rTLS),
            std::get<1>(rTLS),
            std::get<2>(rTLS),
            std::get<3>(rTLS),
            rModelPart,
            model_part_names,
            Delta,
            max_id + ParallelUtilities::GetNumThreads() * 1000, // no element or condition is not suppose to have 1000 nodes per element or condition. This is done to avoid re-calculating the max id for conditions
            rOutputGradientVariable);
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
    rModelPart.GetCommunicator().AssembleNonHistoricalData(rOutputGradientVariable);

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtils::CalculateStrainEnergyLinearlyDependentPropertyGradient(
    ModelPart& rModelPart,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rOutputGradientVariable)
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Vector>;

    const auto& r_process_info = rModelPart.GetProcessInfo();

    block_for_each(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        if (rElement.IsActive()) {
            Vector& r_u = std::get<0>(rTLS);
            Vector& r_sensitivity = std::get<1>(rTLS);

            rElement.GetValuesVector(r_u);

            auto& r_properties = rElement.GetProperties();
            const double current_value = r_properties[rPrimalVariable];

            r_properties[rPrimalVariable] = 1.0;
            rElement.CalculateRightHandSide(r_sensitivity, r_process_info);
            r_properties[rPrimalVariable] = current_value;

            // now calculate the sensitivity
            rElement.GetProperties().SetValue(rOutputGradientVariable, 0.5 * inner_prod(r_u, r_sensitivity));
        } else {
            rElement.GetProperties().SetValue(rOutputGradientVariable, 0.0);
        }
    });

    KRATOS_CATCH("");
}

void LinearStrainEnergyResponseUtils::CalculateStrainEnergySemiAnalyticPropertyGradient(
    ModelPart& rModelPart,
    const double Delta,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rOutputGradientVariable)
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Vector, Vector>;

    const auto& r_process_info = rModelPart.GetProcessInfo();

    block_for_each(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        if (rElement.IsActive()) {
            Vector& r_u = std::get<0>(rTLS);
            Vector& r_ref_rhs = std::get<1>(rTLS);
            Vector& r_perturbed_rhs = std::get<2>(rTLS);

            rElement.GetValuesVector(r_u);

            // calculate the reference value
            rElement.CalculateRightHandSide(r_ref_rhs, r_process_info);

            // now calculate perturbed
            auto& r_properties = rElement.GetProperties();
            r_properties[rPrimalVariable] += Delta;
            rElement.CalculateRightHandSide(r_perturbed_rhs, r_process_info);
            r_properties[rPrimalVariable] -= Delta;

            // now calculate the sensitivity
            rElement.GetProperties().SetValue(rOutputGradientVariable, 0.5 * inner_prod(r_u, r_perturbed_rhs - r_ref_rhs) / Delta);
        } else {
            rElement.GetProperties().SetValue(rOutputGradientVariable, 0.0);
        }
    });

    KRATOS_CATCH("");
}

}