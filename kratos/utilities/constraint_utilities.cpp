//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/constraint_utilities.h"

namespace Kratos
{
namespace ConstraintUtilities
{
void ResetSlaveDofs(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of constraints
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Setting to zero the slave dofs
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                it_const->ResetSlaveDofs(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ApplyConstraints(ModelPart& rModelPart)
{
    KRATOS_TRY

    // The number of constraints
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

    // The current process info
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

    // Adding MPC contribution
    #pragma omp parallel
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = rModelPart.MasterSlaveConstraints().begin() + i_const;

            // Detect if the constraint is active or not. If the user did not make any choice the constraint
            // It is active by default
            bool constraint_is_active = true;
            if (it_const->IsDefined(ACTIVE))
                constraint_is_active = it_const->Is(ACTIVE);

            if (constraint_is_active) {
                it_const->Apply(r_current_process_info);
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void PreComputeExplicitConstraintConstribution(
    ModelPart& rModelPart,
    const std::string DofVariableName,
    const std::string ResidualDofVariableName
    )
{
    KRATOS_TRY

    typedef ModelPart::VariableComponentType VariableComponentType;

    // Getting the dof to check
    const VariableComponentType& r_check_dof_x = KratosComponents<VariableComponentType>::Get(DofVariableName + "_X");
    const VariableComponentType& r_check_dof_y = KratosComponents<VariableComponentType>::Get(DofVariableName + "_Y");
    const VariableComponentType& r_check_dof_z = KratosComponents<VariableComponentType>::Get(DofVariableName + "_Z");

    // Getting the residual dofs
    const VariableComponentType& r_residual_dof_x = KratosComponents<VariableComponentType>::Get(ResidualDofVariableName + "_X");
    const VariableComponentType& r_residual_dof_y = KratosComponents<VariableComponentType>::Get(ResidualDofVariableName + "_Y");
    const VariableComponentType& r_residual_dof_z = KratosComponents<VariableComponentType>::Get(ResidualDofVariableName + "_Z");

    // Getting auxiliar variables
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    const auto it_cont_begin = rModelPart.MasterSlaveConstraints().begin();

    // Auxiliar values
    Matrix transformation_matrix;
    Vector constant_vector;

    #pragma omp parallel firstprivate(transformation_matrix, constant_vector)
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = it_cont_begin + i_const;

            for (auto& r_dof_j: it_const->GetSlaveDofsVector()) {
                it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);

                if (r_dof_j->GetVariable() == r_check_dof_x) {
                    const double r_current_force_residual_x = rModelPart.pGetNode(r_dof_j->Id())->FastGetSolutionStepValue(r_residual_dof_x);

                    for (auto& r_master_dof: it_const->GetMasterDofsVector()) {
                        if (r_master_dof->GetVariable() == r_check_dof_x) {
                            double& current_master_force_residual_x = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_x);
                            #pragma omp atomic
                            current_master_force_residual_x += r_current_force_residual_x * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_y) {
                            double& current_master_force_residual_y = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_y);
                            #pragma omp atomic
                            current_master_force_residual_y += r_current_force_residual_x * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_z) {
                            double& current_master_force_residual_z = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_z);
                            #pragma omp atomic
                            current_master_force_residual_z += r_current_force_residual_x * std::abs(transformation_matrix(0,0));
                        } else continue;
                    }
                } else if (r_dof_j->GetVariable() == r_check_dof_y) {
                    const double r_current_force_residual_y = rModelPart.pGetNode(r_dof_j->Id())->FastGetSolutionStepValue(r_residual_dof_y);

                    for (auto& r_master_dof: it_const->GetMasterDofsVector()) {
                        if (r_master_dof->GetVariable() == r_check_dof_x) {
                            double& current_master_force_residual_x = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_x);
                            #pragma omp atomic
                            current_master_force_residual_x += r_current_force_residual_y * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_y) {
                            double& current_master_force_residual_y = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_y);
                            #pragma omp atomic
                            current_master_force_residual_y +=  r_current_force_residual_y * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_z) {
                            double& current_master_force_residual_z = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_z);
                            #pragma omp atomic
                            current_master_force_residual_z += r_current_force_residual_y * std::abs(transformation_matrix(0,0));
                        } else continue;
                    }
                } else if (r_dof_j->GetVariable() == r_check_dof_z) {
                    const double r_current_force_residual_z = rModelPart.pGetNode(r_dof_j->Id())->FastGetSolutionStepValue(r_residual_dof_z);

                    for (auto& r_master_dof: it_const->GetMasterDofsVector()) {
                        if (r_master_dof->GetVariable() == r_check_dof_x) {
                            double& current_master_force_residual_x = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_x);
                            #pragma omp atomic
                            current_master_force_residual_x += r_current_force_residual_z * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_y) {
                            double& current_master_force_residual_y = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_y);
                            #pragma omp atomic
                            current_master_force_residual_y += r_current_force_residual_z * std::abs(transformation_matrix(0,0));
                        } else if (r_master_dof->GetVariable() == r_check_dof_z) {
                            double& current_master_force_residual_z = rModelPart.pGetNode(r_master_dof->Id())->FastGetSolutionStepValue(r_residual_dof_z);
                            #pragma omp atomic
                            current_master_force_residual_z += r_current_force_residual_z * std::abs(transformation_matrix(0,0));
                        } else continue;
                    }
                } else continue;
            }
        }
    }

    KRATOS_CATCH("")
}

} // namespace ConstraintUtilities
} // namespace Kratos
