// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "custom_utilities/node_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geomechanics_time_integration_scheme.hpp"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class BackwardEulerScheme : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    using BaseType          = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType     = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    BackwardEulerScheme(const std::vector<FirstOrderScalarVariable>&  rFirstOrderScalarVariables,
                        const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(rFirstOrderScalarVariables, rSecondOrderVectorVariables)
    {
    }

    void Predict(ModelPart& rModelPart, DofsArrayType&, TSystemMatrixType&, TSystemVectorType&, TSystemVectorType&) override
    {
        this->PredictVariables(rModelPart);
        this->UpdateVariablesDerivatives(rModelPart);
    }

    void PredictVariables(const ModelPart& rModelPart)
    {
        block_for_each(rModelPart.Nodes(), [this](Node& rNode) { PredictVariablesForNode(rNode); });
    }

    void PredictVariablesForNode(Node& rNode)
    {
        for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
            if (!rNode.SolutionStepsDataHas(r_first_order_scalar_variable.instance)) continue;
            PredictFirstOrderScalarVariableForNode(rNode, r_first_order_scalar_variable);
        }

        // Only for a dynamic scheme, we would need to add the prediction of the second
        // order vector variables.
    }

    void PredictFirstOrderScalarVariableForNode(Node& rNode, const FirstOrderScalarVariable& rFirstOrderScalarVariable)
    {
        const double previous_variable =
            rNode.FastGetSolutionStepValue(rFirstOrderScalarVariable.instance, 1);
        const double previous_first_time_derivative =
            rNode.FastGetSolutionStepValue(rFirstOrderScalarVariable.first_time_derivative, 1);
        const double current_first_time_derivative =
            rNode.FastGetSolutionStepValue(rFirstOrderScalarVariable.first_time_derivative, 0);

        if (rNode.IsFixed(rFirstOrderScalarVariable.first_time_derivative)) {
            rNode.FastGetSolutionStepValue(rFirstOrderScalarVariable.instance) =
                previous_variable + this->GetDeltaTime() * current_first_time_derivative;
        } else if (!rNode.IsFixed(rFirstOrderScalarVariable.instance)) {
            rNode.FastGetSolutionStepValue(rFirstOrderScalarVariable.instance) =
                previous_variable + this->GetDeltaTime() * previous_first_time_derivative;
        }
    }

protected:
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::SetTimeFactors(rModelPart);

        for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
            rModelPart.GetProcessInfo()[r_first_order_scalar_variable.delta_time_coefficient] =
                1.0 / this->GetDeltaTime();
        }

        if (!this->GetSecondOrderVectorVariables().empty()) {
            rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] = 1.0 / this->GetDeltaTime();
        }

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            UpdateVectorTimeDerivatives(rNode);

            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                if (rNode.IsFixed(r_first_order_scalar_variable.first_time_derivative)) continue;

                rNode.FastGetSolutionStepValue(r_first_order_scalar_variable.first_time_derivative) =
                    CalculateDerivative(r_first_order_scalar_variable.instance, rNode);
            }
        });

        KRATOS_CATCH("")
    }

    void UpdateVectorTimeDerivatives(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            const auto updated_first_time_derivative =
                CalculateDerivative(r_second_order_vector_variable.instance, rNode);

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.first_time_derivative, updated_first_time_derivative);

            // Make sure that setting the second_time_derivative is done
            // after setting the first_time_derivative.
            const auto updated_second_time_derivative =
                CalculateDerivative(r_second_order_vector_variable.first_time_derivative, rNode);

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.second_time_derivative, updated_second_time_derivative);
        }
    }

    template <class T>
    T CalculateDerivative(const Variable<T>& rInstanceVariable, Node& rNode) const
    {
        return (rNode.FastGetSolutionStepValue(rInstanceVariable, 0) -
                rNode.FastGetSolutionStepValue(rInstanceVariable, 1)) /
               this->GetDeltaTime();
    }
};

} // namespace Kratos
