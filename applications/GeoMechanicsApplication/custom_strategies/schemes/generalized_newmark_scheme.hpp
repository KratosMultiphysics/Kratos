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

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class GeneralizedNewmarkScheme
    : public GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>
{
public:
    GeneralizedNewmarkScheme(double theta,
                             const std::vector<FirstOrderVariable>& rFirstOrderVariables,
                             const std::vector<SecondOrderVariable> rSecondOrderVariables,
                             double beta = 0.25,
                             double gamma = 0.5)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderVariables, rSecondOrderVariables),
          mTheta(theta),
          mBeta(beta),
          mGamma(gamma)
    {
        KRATOS_ERROR_IF(this->mTheta <= 0)
            << "Theta must be larger than zero, but got " << this->mTheta << "\n";
        KRATOS_ERROR_IF(mBeta <= 0)
            << "Beta must be larger than zero, but got " << mBeta << "\n";
        KRATOS_ERROR_IF(mGamma <= 0)
            << "Gamma must be larger than zero, but got " << mGamma << "\n";
    }

protected:
    void UpdateVectorFirstTimeDerivative(Node& rNode) const override
    {
        for (const auto& variable_derivative : this->GetSecondOrderVariables())
        {
            if (!rNode.SolutionStepsDataHas(variable_derivative.instance))
                continue;

            noalias(rNode.FastGetSolutionStepValue(
                variable_derivative.first_time_derivative, 0)) =
                rNode.FastGetSolutionStepValue(variable_derivative.first_time_derivative, 1) +
                (1.0 - mGamma) * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(
                        variable_derivative.second_time_derivative, 1) +
                mGamma * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(
                        variable_derivative.second_time_derivative, 0);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const override
    {
        for (const auto& variable_derivative : this->GetSecondOrderVariables())
        {
            if (!rNode.SolutionStepsDataHas(variable_derivative.instance))
                continue;

            noalias(rNode.FastGetSolutionStepValue(
                variable_derivative.second_time_derivative, 0)) =
                ((rNode.FastGetSolutionStepValue(variable_derivative.instance, 0) -
                  rNode.FastGetSolutionStepValue(variable_derivative.instance, 1)) -
                 this->GetDeltaTime() * rNode.FastGetSolutionStepValue(
                                            variable_derivative.first_time_derivative, 1) -
                 (0.5 - mBeta) * this->GetDeltaTime() * this->GetDeltaTime() *
                     rNode.FastGetSolutionStepValue(
                         variable_derivative.second_time_derivative, 1)) /
                (mBeta * this->GetDeltaTime() * this->GetDeltaTime());
        }
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const auto delta_time = rModelPart.GetProcessInfo()[DELTA_TIME];
        this->SetDeltaTime(delta_time);

        for (const auto& r_first_order_variable : this->GetFirstOrderVariables())
        {
            rModelPart.GetProcessInfo()[r_first_order_variable.delta_time_coefficient] =
                1.0 / (mTheta * delta_time);
        }

        if (!this->GetSecondOrderVariables().empty())
        {
            rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] =
                mGamma / (mBeta * this->GetDeltaTime());
        }


        KRATOS_CATCH("")
    }

    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const override
    {
        const auto delta_variable = rNode.FastGetSolutionStepValue(variable, 0) -
                                    rNode.FastGetSolutionStepValue(variable, 1);
        const auto previous_dt_variable = rNode.FastGetSolutionStepValue(dt_variable, 1);

        rNode.FastGetSolutionStepValue(dt_variable, 0) =
            (delta_variable - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_variable) /
            (mTheta * this->GetDeltaTime());
    }

    double GetBeta() const
    {
        return mBeta;
    }

    double GetGamma() const
    {
        return mGamma;
    }

private:
    double mTheta = 1.0;
    double mBeta = 0.25;
    double mGamma = 0.5;
};

} // namespace Kratos
