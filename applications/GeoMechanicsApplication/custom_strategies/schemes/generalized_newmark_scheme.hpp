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
#include "custom_utilities/variables_utilities.hpp"
#include "geomechanics_time_integration_scheme.hpp"
#include "includes/model_part.h"
#include <optional>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkScheme : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    using BaseType          = Scheme<TSparseSpace, TDenseSpace>;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables,
                             const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables,
                             std::optional<double> beta,
                             std::optional<double> gamma,
                             std::optional<double> theta)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(rFirstOrderScalarVariables,
                                                                       rSecondOrderVectorVariables),
          mBeta(beta),
          mGamma(gamma),
          mTheta(theta)
    {
        KRATOS_ERROR_IF(!rSecondOrderVectorVariables.empty() && !mBeta.has_value())
            << "Beta must be set when second order vector variables are used\n";
        KRATOS_ERROR_IF(!rSecondOrderVectorVariables.empty() && !mGamma.has_value())
            << "Gamma must be set when second order vector variables are used\n";
        KRATOS_ERROR_IF(!rFirstOrderScalarVariables.empty() && !mTheta.has_value())
            << "Theta must be set when first order scalar variables are used\n";

        KRATOS_ERROR_IF(mBeta.has_value() && mBeta.value() <= 0)
            << "Beta must be larger than zero, but got " << mBeta.value() << "\n";
        KRATOS_ERROR_IF(mGamma.has_value() && mGamma.value() <= 0)
            << "Gamma must be larger than zero, but got " << mGamma.value() << "\n";
        KRATOS_ERROR_IF(mTheta.has_value() && mTheta.value() <= 0)
            << "Theta must be larger than zero, but got " << mTheta.value() << "\n";
    }

    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderScalarVariables, {}, std::nullopt, std::nullopt, theta)
    {
    }

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb) override
    {
        KRATOS_TRY

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING]) {
            // Clear nodal variables
            block_for_each(rModelPart.Nodes(),
                           [](Node& rNode) { rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0; });

            this->FinalizeSolutionStepActiveEntities(rModelPart, rA, rDx, rb);
        } else {
            this->FinalizeSolutionStepActiveEntities(rModelPart, rA, rDx, rb);
        }

        KRATOS_CATCH("")
    }

protected:
    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::SetTimeFactors(rModelPart);

        for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
            rModelPart.GetProcessInfo()[r_first_order_scalar_variable.delta_time_coefficient] =
                1.0 / (GetTheta() * this->GetDeltaTime());
        }

        if (!this->GetSecondOrderVectorVariables().empty()) {
            rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] =
                GetGamma() / (GetBeta() * this->GetDeltaTime());
        }

        KRATOS_CATCH("")
    }

    [[nodiscard]] double GetBeta() const { return mBeta.value(); }

    [[nodiscard]] double GetGamma() const { return mGamma.value(); }

    [[nodiscard]] double GetTheta() const { return mTheta.value(); }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // For the Newmark schemes the second derivatives should be updated before calculating the first derivatives
            this->UpdateVectorSecondTimeDerivative(rNode);
            this->UpdateVectorFirstTimeDerivative(rNode);

            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                this->UpdateScalarTimeDerivative(rNode, r_first_order_scalar_variable.instance,
                                                 r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

    void UpdateScalarTimeDerivative(Node& rNode, const Variable<double>& variable, const Variable<double>& dt_variable) const
    {
        if (rNode.IsFixed(dt_variable)) return;

        const auto delta_variable =
            rNode.FastGetSolutionStepValue(variable, 0) - rNode.FastGetSolutionStepValue(variable, 1);
        const auto previous_dt_variable = rNode.FastGetSolutionStepValue(dt_variable, 1);

        rNode.FastGetSolutionStepValue(dt_variable, 0) =
            (delta_variable - (1.0 - GetTheta()) * this->GetDeltaTime() * previous_dt_variable) /
            (GetTheta() * this->GetDeltaTime());
    }

    void UpdateVectorFirstTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            const array_1d<double, 3> updated_first_derivative =
                rNode.FastGetSolutionStepValue(r_second_order_vector_variable.first_time_derivative, 1) +
                (1.0 - GetGamma()) * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1) +
                GetGamma() * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 0);

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.first_time_derivative, updated_first_derivative);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            const array_1d<double, 3> updated_second_time_derivative =
                ((rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 0) -
                  rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 1)) -
                 this->GetDeltaTime() * rNode.FastGetSolutionStepValue(
                                            r_second_order_vector_variable.first_time_derivative, 1) -
                 (0.5 - GetBeta()) * this->GetDeltaTime() * this->GetDeltaTime() *
                     rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1)) /
                (GetBeta() * this->GetDeltaTime() * this->GetDeltaTime());

            NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(
                rNode, r_second_order_vector_variable.second_time_derivative, updated_second_time_derivative);
        }
    }

private:
    std::optional<double> mBeta;
    std::optional<double> mGamma;
    std::optional<double> mTheta;
};

} // namespace Kratos
