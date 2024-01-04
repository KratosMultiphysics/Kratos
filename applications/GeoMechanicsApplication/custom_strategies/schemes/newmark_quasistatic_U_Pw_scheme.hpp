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
              theta,
              {FirstOrderVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {VariableWithTimeDerivatives(DISPLACEMENT),
               VariableWithTimeDerivatives{ROTATION}}, beta, gamma)
    {
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
            block_for_each(rModelPart.Nodes(), [&stress_tensor_size](Node& rNode)
            {
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
            block_for_each(rModelPart.Nodes(), [](Node& rNode)
            {
                if (const double& nodal_area = rNode.FastGetSolutionStepValue(NODAL_AREA);
                    nodal_area > 1.0e-20)
                {
                    const double inv_nodal_area = 1.0 / nodal_area;
                    rNode.FastGetSolutionStepValue(NODAL_CAUCHY_STRESS_TENSOR) *= inv_nodal_area;
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
}; // Class NewmarkQuasistaticUPwScheme

} // namespace Kratos