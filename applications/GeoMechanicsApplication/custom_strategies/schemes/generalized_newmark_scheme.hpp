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

#include "geomechanics_time_integration_scheme.hpp"
#include <optional>

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkScheme : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables,
                             const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables,
                             std::optional<double> theta,
                             std::optional<double> beta,
                             std::optional<double> gamma)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(rFirstOrderScalarVariables,
                                                                       rSecondOrderVectorVariables),
          mTheta(theta),
          mBeta(beta),
          mGamma(gamma)
    {
        KRATOS_ERROR_IF(!rFirstOrderScalarVariables.empty() && !mTheta.has_value())
            << "Theta must be set when first order scalar variables are used\n";
        KRATOS_ERROR_IF(!rSecondOrderVectorVariables.empty() && !mBeta.has_value())
            << "Beta must be set when second order vector variables are used\n";
        KRATOS_ERROR_IF(!rSecondOrderVectorVariables.empty() && !mGamma.has_value())
            << "Gamma must be set when second order vector variables are used\n";

        KRATOS_ERROR_IF(mTheta.has_value() && mTheta.value() <= 0)
            << "Theta must be larger than zero, but got " << mTheta.value() << "\n";
        KRATOS_ERROR_IF(mBeta.has_value() && mBeta.value() <= 0)
            << "Beta must be larger than zero, but got " << mBeta.value() << "\n";
        KRATOS_ERROR_IF(mGamma.has_value() && mGamma.value() <= 0)
            << "Gamma must be larger than zero, but got " << mGamma.value() << "\n";
    }

    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderScalarVariables, {}, theta, std::nullopt, std::nullopt)
    {
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // For the Newmark schemes the second derivatives should be updated before calculating the first derivatives
            UpdateVectorSecondTimeDerivative(rNode);
            UpdateVectorFirstTimeDerivative(rNode);

            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                UpdateScalarTimeDerivative(rNode, r_first_order_scalar_variable.instance,
                                           r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

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

    double GetBeta() const { return mBeta.value(); }

    double GetGamma() const { return mGamma.value(); }

    double GetTheta() const { return mTheta.value(); }

private:
    void UpdateScalarTimeDerivative(Node& rNode, const Variable<double>& variable, const Variable<double>& dt_variable) const
    {
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

            noalias(rNode.FastGetSolutionStepValue(r_second_order_vector_variable.first_time_derivative, 0)) =
                rNode.FastGetSolutionStepValue(r_second_order_vector_variable.first_time_derivative, 1) +
                (1.0 - GetGamma()) * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1) +
                GetGamma() * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 0);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;

            noalias(rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 0)) =
                ((rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 0) -
                  rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 1)) -
                 this->GetDeltaTime() * rNode.FastGetSolutionStepValue(
                                            r_second_order_vector_variable.first_time_derivative, 1) -
                 (0.5 - GetBeta()) * this->GetDeltaTime() * this->GetDeltaTime() *
                     rNode.FastGetSolutionStepValue(r_second_order_vector_variable.second_time_derivative, 1)) /
                (GetBeta() * this->GetDeltaTime() * this->GetDeltaTime());
        }
    }

    std::optional<double> mTheta;
    std::optional<double> mBeta;
    std::optional<double> mGamma;
};

} // namespace Kratos
