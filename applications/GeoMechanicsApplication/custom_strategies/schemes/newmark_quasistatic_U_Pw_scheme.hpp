// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include "generalized_newmark_scheme.hpp"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

struct VariableWithTimeDerivatives
{
    Variable<array_1d<double, 3>> instance;
    Variable<array_1d<double, 3>> first_time_derivative;
    Variable<array_1d<double, 3>> second_time_derivative;

    explicit VariableWithTimeDerivatives(const Variable<array_1d<double, 3>>& instance)
        : instance(instance),
          first_time_derivative(instance.GetTimeDerivative()),
          second_time_derivative(first_time_derivative.GetTimeDerivative())
    {
    }
};

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticUPwScheme
    : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticUPwScheme);

    using BaseType = Scheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              theta, WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT),
          mBeta(beta),
          mGamma(gamma)
    {
        KRATOS_ERROR_IF(mBeta <= 0)
            << "Beta must be larger than zero, but got " << mBeta << "\n";
        KRATOS_ERROR_IF(mGamma <= 0)
            << "Gamma must be larger than zero, but got " << mGamma << "\n";
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              TSystemMatrixType& A,
                              TSystemVectorType& Dx,
                              TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING])
        {
            const unsigned int dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
            const auto stress_tensor_size =
                dim == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D;

            // Clear nodal variables
            block_for_each(
                rModelPart.Nodes(),
                [&stress_tensor_size](Node& rNode)
                {
                    rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                    Matrix& r_nodal_stress =
                        rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    if (r_nodal_stress.size1() != stress_tensor_size)
                        r_nodal_stress.resize(stress_tensor_size, stress_tensor_size, false);
                    noalias(r_nodal_stress) =
                        ZeroMatrix(stress_tensor_size, stress_tensor_size);
                    rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
                });

            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);

            // Compute smoothed nodal variables
            block_for_each(
                rModelPart.Nodes(),
                [](Node& rNode)
                {
                    if (const double& nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                        nodal_area > 1.0e-20)
                    {
                        const double inv_nodal_area = 1.0 / nodal_area;
                        Matrix& r_nodal_stress =
                            rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                        for (unsigned int i = 0; i < r_nodal_stress.size1(); ++i)
                        {
                            for (unsigned int j = 0; j < r_nodal_stress.size2(); ++j)
                            {
                                r_nodal_stress(i, j) *= inv_nodal_area;
                            }
                        }
                        rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) *= inv_nodal_area;
                    }

                    if (const double& nodal_joint_area =
                            rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA);
                        nodal_joint_area > 1.0e-20)
                    {
                        const double inv_nodal_joint_area = 1.0 / nodal_joint_area;
                        rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) *= inv_nodal_joint_area;
                        rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) *= inv_nodal_joint_area;
                    }
                });
        }
        else
        {
            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
        }

        KRATOS_CATCH("")
    }

protected:
    double mBeta = 0.25;
    double mGamma = 0.5;

    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>::CheckAllocatedVariables(rModelPart);

        for (const auto& r_node : rModelPart.Nodes())
        {
            for (const auto& variable_derivative : mVariableDerivatives)
            {
                if (!rModelPart.HasNodalSolutionStepVariable(variable_derivative.instance))
                    continue;

                this->CheckSolutionStepsData(r_node, variable_derivative.instance);
                this->CheckSolutionStepsData(r_node, variable_derivative.first_time_derivative);
                this->CheckSolutionStepsData(r_node, variable_derivative.second_time_derivative);

                std::vector<std::string> components{"X", "Y", "Z"};
                for (const auto& component : components)
                {
                    const auto& variable_component = GetComponentFromVectorVariable(
                        variable_derivative.instance, component);
                    this->CheckDof(r_node, variable_component);
                }
            }
        }
    }

    inline void SetTimeFactors(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>::SetTimeFactors(rModelPart);
        rModelPart.GetProcessInfo()[VELOCITY_COEFFICIENT] =
            mGamma / (mBeta * this->GetDeltaTime());

        KRATOS_CATCH("")
    }

    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        // Update Acceleration, Velocity and DtPressure
        block_for_each(rModelPart.Nodes(),
                       [this](Node& rNode)
                       {
                           UpdateVectorSecondTimeDerivative(rNode);
                           UpdateVectorFirstTimeDerivative(rNode);
                           this->UpdateScalarTimeDerivative(
                               rNode, WATER_PRESSURE, DT_WATER_PRESSURE);
                       });

        KRATOS_CATCH("")
    }

    const Variable<double>& GetComponentFromVectorVariable(
        const Variable<array_1d<double, 3>>& rSource, const std::string& rComponent) const
    {
        return KratosComponents<Variable<double>>::Get(rSource.Name() + "_" + rComponent);
    }

    const std::vector<VariableWithTimeDerivatives>& GetVariableDerivatives() const
    {
        return mVariableDerivatives;
    }

private:
    void UpdateVectorFirstTimeDerivative(Node& rNode) const
    {
        for (const auto& variable_derivative : mVariableDerivatives)
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

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        for (const auto& variable_derivative : mVariableDerivatives)
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

    std::vector<VariableWithTimeDerivatives> mVariableDerivatives{
        VariableWithTimeDerivatives(DISPLACEMENT), VariableWithTimeDerivatives{ROTATION}};

}; // Class NewmarkQuasistaticUPwScheme

} // namespace Kratos