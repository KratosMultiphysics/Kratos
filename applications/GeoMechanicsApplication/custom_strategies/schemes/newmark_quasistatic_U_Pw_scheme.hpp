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

namespace Kratos {

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticUPwScheme
    : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace> {
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

protected:
    void CheckAllocatedVariables(const ModelPart& rModelPart) const override
    {
        GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>::CheckAllocatedVariables(rModelPart);

        for (const auto& r_node : rModelPart.Nodes()) {
            this->CheckSolutionStepsData(r_node, DISPLACEMENT);
            this->CheckSolutionStepsData(r_node, VELOCITY);
            this->CheckSolutionStepsData(r_node, ACCELERATION);

            this->CheckDof(r_node, DISPLACEMENT_X);
            this->CheckDof(r_node, DISPLACEMENT_Y);
            this->CheckDof(r_node, DISPLACEMENT_Z);
        }
    }

public:
    void FinalizeSolutionStep(ModelPart& rModelPart,
                              TSystemMatrixType& A,
                              TSystemVectorType& Dx,
                              TSystemVectorType& b) override
    {
        KRATOS_TRY

        if (rModelPart.GetProcessInfo()[NODAL_SMOOTHING]) {
            unsigned int dim = rModelPart.GetProcessInfo()[DOMAIN_SIZE];

            SizeType stress_tensor_size = STRESS_TENSOR_SIZE_2D;
            if (dim == N_DIM_3D)
                stress_tensor_size = STRESS_TENSOR_SIZE_3D;

            // Clear nodal variables
            block_for_each(rModelPart.Nodes(), [&stress_tensor_size](Node& rNode) {
                rNode.FastGetSolutionStepValue(NODAL_AREA) = 0.0;
                Matrix& r_nodal_stress =
                    rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                if (r_nodal_stress.size1() != stress_tensor_size)
                    r_nodal_stress.resize(stress_tensor_size, stress_tensor_size, false);
                noalias(r_nodal_stress) = ZeroMatrix(stress_tensor_size, stress_tensor_size);
                rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) = 0.0;
                rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) = 0.0;
            });

            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);

            // Compute smoothed nodal variables
            block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
                const double& nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                if (nodal_area > 1.0e-20) {
                    const double inv_nodal_area = 1.0 / nodal_area;
                    Matrix& r_nodal_stress =
                        rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR);
                    for (unsigned int i = 0; i < r_nodal_stress.size1(); ++i) {
                        for (unsigned int j = 0; j < r_nodal_stress.size2(); ++j) {
                            r_nodal_stress(i, j) *= inv_nodal_area;
                        }
                    }
                    rNode.FastGetSolutionStepValue(NODAL_DAMAGE_VARIABLE) *= inv_nodal_area;
                }

                const double& nodal_joint_area =
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_AREA);
                if (nodal_joint_area > 1.0e-20) {
                    const double inv_nodal_joint_area = 1.0 / nodal_joint_area;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_WIDTH) *= inv_nodal_joint_area;
                    rNode.FastGetSolutionStepValue(NODAL_JOINT_DAMAGE) *= inv_nodal_joint_area;
                }
            });
        }
        else {
            this->FinalizeSolutionStepActiveEntities(rModelPart, A, Dx, b);
        }

        KRATOS_CATCH("")
    }

protected:
    double mBeta = 0.25;
    double mGamma = 0.5;

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
        block_for_each(rModelPart.Nodes(), [&](Node& rNode) {
            UpdateVectorSecondTimeDerivative(rNode);
            UpdateVectorFirstTimeDerivative(rNode);
            this->UpdateScalarTimeDerivative(rNode, WATER_PRESSURE, DT_WATER_PRESSURE);
        });

        KRATOS_CATCH("")
    }

private:
    void UpdateVectorFirstTimeDerivative(Node& rNode) const
    {
        noalias(rNode.FastGetSolutionStepValue(VELOCITY, 0)) =
            rNode.FastGetSolutionStepValue(VELOCITY, 1) +
            (1.0 - mGamma) * this->GetDeltaTime() *
                rNode.FastGetSolutionStepValue(ACCELERATION, 1) +
            mGamma * this->GetDeltaTime() * rNode.FastGetSolutionStepValue(ACCELERATION, 0);
    }

    void UpdateVectorSecondTimeDerivative(Node& rNode) const
    {
        noalias(rNode.FastGetSolutionStepValue(ACCELERATION, 0)) =
            ((rNode.FastGetSolutionStepValue(DISPLACEMENT, 0) -
              rNode.FastGetSolutionStepValue(DISPLACEMENT, 1)) -
             this->GetDeltaTime() * rNode.FastGetSolutionStepValue(VELOCITY, 1) -
             (0.5 - mBeta) * this->GetDeltaTime() * this->GetDeltaTime() *
                 rNode.FastGetSolutionStepValue(ACCELERATION, 1)) /
            (mBeta * this->GetDeltaTime() * this->GetDeltaTime());
    }
}; // Class NewmarkQuasistaticUPwScheme

} // namespace Kratos