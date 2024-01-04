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
    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables,
                             const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables,
                             double theta,
                             double beta,
                             double gamma)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderScalarVariables, rSecondOrderVectorVariables),
          mTheta(theta),
          mBeta(beta),
          mGamma(gamma)
    {
        KRATOS_ERROR_IF(mTheta <= 0)
            << "Theta must be larger than zero, but got " << mTheta << "\n";
        KRATOS_ERROR_IF(mBeta <= 0)
            << "Beta must be larger than zero, but got " << mBeta << "\n";
        KRATOS_ERROR_IF(mGamma <= 0)
            << "Gamma must be larger than zero, but got " << mGamma << "\n";
    }

    GeneralizedNewmarkScheme(const std::vector<FirstOrderScalarVariable>& rFirstOrderScalarVariables,
                             double theta)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              rFirstOrderScalarVariables, {}),
          mTheta(theta)
    {
        KRATOS_ERROR_IF(mTheta <= 0)
            << "Theta must be larger than zero, but got " << mTheta << "\n";
    }

    GeneralizedNewmarkScheme(const std::vector<SecondOrderVectorVariable>& rSecondOrderVectorVariables,
                             double beta,
                             double gamma)
        : GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>(
              {}, rSecondOrderVectorVariables),
          mBeta(beta),
          mGamma(gamma)
    {
        KRATOS_ERROR_IF(mBeta <= 0)
            << "Beta must be larger than zero, but got " << mBeta << "\n";
        KRATOS_ERROR_IF(mGamma <= 0)
            << "Gamma must be larger than zero, but got " << mGamma << "\n";
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode)
        {
            // For the Newmark schemes the second derivatives should be updated before calculating the first derivatives
            UpdateVectorSecondTimeDerivative(rNode);
            UpdateVectorFirstTimeDerivative(rNode);

            for (const auto& r_first_order_scalar_variable :
                 this->GetFirstOrderScalarVariables())
            {
                UpdateScalarTimeDerivative(
                    rNode, r_first_order_scalar_variable.instance,
                    r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        GeoMechanicsTimeIntegrationScheme<TSparseSpace, TDenseSpace>::SetTimeFactors(rModelPart);

        for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables())
        {
            rModelPart.GetProcessInfo()[r_first_order_scalar_variable.delta_time_coefficient] =
                1.0 / (mTheta * this->GetDeltaTime());
        }

        if (!this->GetSecondOrderVectorVariables().empty())
        {
            rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] =
                mGamma / (mBeta * this->GetDeltaTime());
        }

        KRATOS_CATCH("")
    }

    double GetBeta() const { return mBeta; }
    double GetGamma() const { return mGamma; }

private:
    void UpdateScalarTimeDerivative(Node& rNode,
                                    const Variable<double>& variable,
                                    const Variable<double>& dt_variable) const
    {
        const auto delta_variable = rNode.FastGetSolutionStepValue(variable, 0) -
                                    rNode.FastGetSolutionStepValue(variable, 1);
        const auto previous_dt_variable = rNode.FastGetSolutionStepValue(dt_variable, 1);

        rNode.FastGetSolutionStepValue(dt_variable, 0) =
            (delta_variable - (1.0 - mTheta) * this->GetDeltaTime() * previous_dt_variable) /
            (mTheta * this->GetDeltaTime());
    }

    void UpdateVectorFirstTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable :
             this->GetSecondOrderVectorVariables())
        {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance))
                continue;

            noalias(rNode.FastGetSolutionStepValue(
                r_second_order_vector_variable.first_time_derivative, 0)) =
                rNode.FastGetSolutionStepValue(
                    r_second_order_vector_variable.first_time_derivative, 1) +
                (1.0 - mGamma) * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(
                        r_second_order_vector_variable.second_time_derivative, 1) +
                mGamma * this->GetDeltaTime() *
                    rNode.FastGetSolutionStepValue(
                        r_second_order_vector_variable.second_time_derivative, 0);
        }
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        for (const auto& r_second_order_vector_variable :
             this->GetSecondOrderVectorVariables())
        {
            if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance))
                continue;

            noalias(rNode.FastGetSolutionStepValue(
                r_second_order_vector_variable.second_time_derivative, 0)) =
                ((rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 0) -
                  rNode.FastGetSolutionStepValue(r_second_order_vector_variable.instance, 1)) -
                 this->GetDeltaTime() *
                     rNode.FastGetSolutionStepValue(
                         r_second_order_vector_variable.first_time_derivative, 1) -
                 (0.5 - mBeta) * this->GetDeltaTime() * this->GetDeltaTime() *
                     rNode.FastGetSolutionStepValue(
                         r_second_order_vector_variable.second_time_derivative, 1)) /
                (mBeta * this->GetDeltaTime() * this->GetDeltaTime());
        }
    }

    double mTheta = 1.0;
    double mBeta = 0.25;
    double mGamma = 0.5;
};

} // namespace Kratos
