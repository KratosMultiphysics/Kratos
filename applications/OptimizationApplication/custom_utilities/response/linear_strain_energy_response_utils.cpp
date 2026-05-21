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
#include "expression/variable_expression_io.h"
#include "custom_utilities/properties_variable_expression_io.h"
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
    const PhysicalFieldVariableTypes& rPhysicalVariable,
    ModelPart& rGradientRequiredModelPart,
    ModelPart& rGradientComputedModelPart,
    std::vector<ContainerExpressionType>& rListOfContainerExpressions,
    const double PerturbationSize)
{
    KRATOS_TRY

    std::visit([&](auto pVariable) {
        if (*pVariable == YOUNG_MODULUS) {
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(YOUNG_MODULUS_SENSITIVITY, 0.0); });
            CalculateStrainEnergyLinearlyDependentPropertyGradient(rGradientComputedModelPart, YOUNG_MODULUS, YOUNG_MODULUS_SENSITIVITY);
        } else if (*pVariable == THICKNESS) {
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(THICKNESS_SENSITIVITY, 0.0); });
            CalculateStrainEnergySemiAnalyticPropertyGradient(rGradientComputedModelPart, PerturbationSize, THICKNESS, THICKNESS_SENSITIVITY);
        } else if (*pVariable == POISSON_RATIO) {
            block_for_each(rGradientRequiredModelPart.Elements(), [](auto& rElement) { rElement.GetProperties().SetValue(POISSON_RATIO_SENSITIVITY, 0.0); });
            CalculateStrainEnergySemiAnalyticPropertyGradient(rGradientComputedModelPart, PerturbationSize, POISSON_RATIO, POISSON_RATIO_SENSITIVITY);
        } else if (*pVariable == SHAPE) {
            VariableUtils().SetNonHistoricalVariableToZero(SHAPE_SENSITIVITY, rGradientRequiredModelPart.Nodes());
            CalculateStrainEnergySemiAnalyticShapeGradient(rGradientComputedModelPart, PerturbationSize, SHAPE_SENSITIVITY);
        } else {
            KRATOS_ERROR
                << "Unsupported sensitivity w.r.t. " << pVariable->Name()
                << " requested. Followings are supported sensitivity variables:"
                << "\n\t" << YOUNG_MODULUS.Name()
                << "\n\t" << THICKNESS.Name()
                << "\n\t" << POISSON_RATIO.Name()
                << "\n\t" << SHAPE.Name();
        }

        // now fill the container expressions
        for (auto& p_container_expression : rListOfContainerExpressions) {
            std::visit([pVariable](auto& pContainerExpression){
                using container_type = std::decay_t<decltype(*pContainerExpression)>;

                if (*pVariable == SHAPE) {
                    if constexpr(std::is_same_v<container_type, ContainerExpression<ModelPart::NodesContainerType>>) {
                        VariableExpressionIO::Read(*pContainerExpression, &SHAPE_SENSITIVITY, false);
                    } else {
                        KRATOS_ERROR << "Requesting sensitivity w.r.t. "
                                        "SHAPE for a Container expression "
                                        "which is not a NodalExpression. [ "
                                        "Requested container expression = "
                                        << *pContainerExpression << " ].\n";
                    }
                } else {
                    if constexpr(std::is_same_v<container_type, ContainerExpression<ModelPart::ElementsContainerType>>) {
                        const auto& sensitivity_variable = KratosComponents<Variable<double>>::Get(pVariable->Name() + "_SENSITIVITY");
                        PropertiesVariableExpressionIO::Read(*pContainerExpression, &sensitivity_variable);
                    } else {
                        KRATOS_ERROR << "Requesting sensitivity w.r.t. "
                                     << pVariable->Name()
                                     << " for a Container expression "
                                        "which is not an ElementExpression. [ "
                                        "Requested container expression = "
                                     << *pContainerExpression << " ].\n";
                    }
                }


            }, p_container_expression);
        }
    }, rPhysicalVariable);

    KRATOS_CATCH("");
}

template<class TEntityType>
void LinearStrainEnergyResponseUtils::CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
    TEntityType& rEntity,
    Vector& rX,
    Vector& rRefRHS,
    Vector& rPerturbedRHS,
    Node::Pointer& pThreadLocalNode,
    ModelPart& rModelPart,
    const double Delta,
    const Variable<array_1d<double, 3>>& rOutputGradientVariable)
{
    KRATOS_TRY

    if (rEntity.IsActive()) {
        const auto& r_process_info = rModelPart.GetProcessInfo();
        auto& r_geometry = rEntity.GetGeometry();
        const auto domain_size = r_geometry.WorkingSpaceDimension();

        rEntity.GetValuesVector(rX);

        rX /= 2.0;

        // calculate the reference value
        rEntity.CalculateRightHandSide(rRefRHS, r_process_info);

        // initialize dummy node for parallelized perturbation based sensitivity calculation
        // in each thread separately
        if (!pThreadLocalNode) {
            pThreadLocalNode = Kratos::make_intrusive<Node>(1, 0, 0, 0);
        }

        // now calculate perturbed
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            auto& r_orig_node_sensitivity = r_geometry[i].GetValue(rOutputGradientVariable);

            // get the geometry node pointer
            auto& p_node = rEntity.GetGeometry()(i);

            // now copy the node data to thread local node using the operator= in Node
            (*pThreadLocalNode) = (*p_node);

            // now swap entity node with the thread local node
            std::swap(p_node, pThreadLocalNode);

            // now do the finite difference computation.
            auto& r_coordinates = p_node->Coordinates();
            auto& r_initial_coordinates = p_node->GetInitialPosition();

            r_initial_coordinates[0] += Delta;
            r_coordinates[0] += Delta;
            rEntity.CalculateRightHandSide(rPerturbedRHS, r_process_info);
            r_initial_coordinates[0] -= Delta;
            r_coordinates[0] -= Delta;
            AtomicAdd<double>(r_orig_node_sensitivity[0], inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);

            r_initial_coordinates[1] += Delta;
            r_coordinates[1] += Delta;
            rEntity.CalculateRightHandSide(rPerturbedRHS, r_process_info);
            r_initial_coordinates[1] -= Delta;
            r_coordinates[1] -= Delta;
            AtomicAdd<double>(r_orig_node_sensitivity[1], inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);

            if (domain_size == 3) {
                r_initial_coordinates[2] += Delta;
                r_coordinates[2] += Delta;
                rEntity.CalculateRightHandSide(rPerturbedRHS, r_process_info);
                r_initial_coordinates[2] -= Delta;
                r_coordinates[2] -= Delta;
                AtomicAdd<double>(r_orig_node_sensitivity[2], inner_prod(rX, rPerturbedRHS - rRefRHS) / Delta);
            }

            // revert back the node change.
            std::swap(p_node, pThreadLocalNode);
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

    using tls_type = std::tuple<Vector, Vector, Vector, Node::Pointer>;

    VariableUtils().SetNonHistoricalVariableToZero(rOutputGradientVariable, rModelPart.Nodes());

    block_for_each(rModelPart.Elements(), tls_type(), [&](auto& rElement, tls_type& rTLS) {
        CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
            rElement,
            std::get<0>(rTLS),
            std::get<1>(rTLS),
            std::get<2>(rTLS),
            std::get<3>(rTLS),
            rModelPart,
            Delta,
            rOutputGradientVariable);
    });

    block_for_each(rModelPart.Conditions(), tls_type(), [&](auto& rCondition, tls_type& rTLS) {
        CalculateStrainEnergyEntitySemiAnalyticShapeGradient(
            rCondition,
            std::get<0>(rTLS),
            std::get<1>(rTLS),
            std::get<2>(rTLS),
            std::get<3>(rTLS),
            rModelPart,
            Delta,
            rOutputGradientVariable);
    });

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