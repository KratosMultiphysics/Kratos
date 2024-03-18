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

    void FinalizeSolutionStep(ModelPart& rModelPart, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING]) {
            const unsigned int dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
            const auto stress_tensor_size = dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D;

            // Clear nodal variables
            block_for_each(rModelPart.Nodes(), [&stress_tensor_size](Node& rNode) {
                rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& r_nodal_stress = rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if (r_nodal_stress.size1() != stress_tensor_size)
                    r_nodal_stress.resize(stress_tensor_size, stress_tensor_size, false);
                noalias(r_nodal_stress) = ZeroMatrix(stress_tensor_size, stress_tensor_size);
                rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA)      = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH)     = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE)    = 0.0;
            });

            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);

            // Compute smoothed nodal variables
            block_for_each(rModelPart.Nodes(), [](Node& rNode) {
                if (const double& nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA); nodal_area > 1.0e-20) {
                    const double inv_nodal_area = 1.0 / nodal_area;
                    rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) *= inv_nodal_area;
                    rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) *= inv_nodal_area;
                }

                if (const double& nodal_joint_area = rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA);
                    nodal_joint_area > 1.0e-20) {
                    const double inv_nodal_joint_area = 1.0 / nodal_joint_area;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) *= inv_nodal_joint_area;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) *= inv_nodal_joint_area;
                }
            });
        } else {
            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
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

private:
    std::optional<double> mBeta;
    std::optional<double> mGamma;
    std::optional<double> mTheta;
};

} // namespace Kratos
