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
    const std::vector<std::string>& rDofVariableNames,
    const std::vector<std::string>& rResidualDofVariableNames
    )
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rDofVariableNames.size() == rResidualDofVariableNames.size()) << "PreComputeExplicitConstraintConstribution not properly defined variables" << std::endl;

    // Defining variable maps
    std::unordered_map<std::size_t, Variable<double>> double_variable_map;
    typedef ModelPart::VariableComponentType VariableComponentType;
    std::unordered_map<std::size_t, VariableComponentType> components_variable_map;

    std::size_t counter = 0;
    for (auto& r_dof_variable_name : rDofVariableNames) {
        const std::string& r_reaction_variable_name = rResidualDofVariableNames[counter];

        if (KratosComponents<Variable<double>>::Has(r_dof_variable_name)) {
            double_variable_map.insert(std::pair<std::size_t, Variable<double>>(KratosComponents<Variable<double>>::Get(r_dof_variable_name).Key(), KratosComponents<Variable<double>>::Get(r_reaction_variable_name)));
        } else if (KratosComponents<VariableComponentType>::Has(r_dof_variable_name)) {
            // Getting the dof to check
            const VariableComponentType& r_check_dof_x = KratosComponents<VariableComponentType>::Get(r_dof_variable_name + "_X");
            const VariableComponentType& r_check_dof_y = KratosComponents<VariableComponentType>::Get(r_dof_variable_name + "_Y");
            const VariableComponentType& r_check_dof_z = KratosComponents<VariableComponentType>::Get(r_dof_variable_name + "_Z");

            // Getting the residual dofs
            const VariableComponentType& r_residual_dof_x = KratosComponents<VariableComponentType>::Get(r_reaction_variable_name + "_X");
            const VariableComponentType& r_residual_dof_y = KratosComponents<VariableComponentType>::Get(r_reaction_variable_name + "_Y");
            const VariableComponentType& r_residual_dof_z = KratosComponents<VariableComponentType>::Get(r_reaction_variable_name + "_Z");

            components_variable_map.insert(std::pair<std::size_t, VariableComponentType>(r_check_dof_x.Key(), r_residual_dof_x));
            components_variable_map.insert(std::pair<std::size_t, VariableComponentType>(r_check_dof_y.Key(), r_residual_dof_y));
            components_variable_map.insert(std::pair<std::size_t, VariableComponentType>(r_check_dof_z.Key(), r_residual_dof_z));
        } else {
            KRATOS_ERROR << "Variable is not a component or a double" << std::endl;
        }

        ++counter;
    }

    // Getting auxiliar variables
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());
    const auto it_cont_begin = rModelPart.MasterSlaveConstraints().begin();

    // Auxiliar values
    Matrix transformation_matrix;
    Vector constant_vector, slave_solution_vector, master_solution_vector;

    #pragma omp parallel firstprivate(transformation_matrix, constant_vector, slave_solution_vector, master_solution_vector)
    {
        #pragma omp for schedule(guided, 512)
        for (int i_const = 0; i_const < number_of_constraints; ++i_const) {
            auto it_const = it_cont_begin + i_const;

            // Getting the transformation matrix and constant vector
            it_const->CalculateLocalSystem(transformation_matrix, constant_vector, r_current_process_info);

            // Resizing vectors
            if (slave_solution_vector.size() != transformation_matrix.size1())
                slave_solution_vector.resize(transformation_matrix.size1());
            if (master_solution_vector.size() != transformation_matrix.size2())
                master_solution_vector.resize(transformation_matrix.size2());

            std::size_t counter = 0;
            for (auto& r_dof_master: it_const->GetMasterDofsVector()) {
                const std::size_t master_variable_key = r_dof_master->GetVariable().Key();

                if (double_variable_map.find(master_variable_key) != double_variable_map.end()) {
                    const auto& r_aux_var = double_variable_map.find(master_variable_key)->second;
                    master_solution_vector[counter] = rModelPart.pGetNode(r_dof_master->Id())->FastGetSolutionStepValue(r_aux_var);
                } else if (components_variable_map.find(master_variable_key) != components_variable_map.end()) {
                    const auto& r_aux_var = components_variable_map.find(master_variable_key)->second;
                    master_solution_vector[counter] = rModelPart.pGetNode(r_dof_master->Id())->FastGetSolutionStepValue(r_aux_var);
                } else {
                    KRATOS_ERROR << "Dof variable is not defined" << std::endl;
                }

                ++counter;
            }

            // Computing transfered solution
            noalias(slave_solution_vector) = prod(transformation_matrix, master_solution_vector);

            counter = 0;
            for (auto& r_dof_slave: it_const->GetSlaveDofsVector()) {
                const std::size_t slave_variable_key = r_dof_slave->GetVariable().Key();

                if (double_variable_map.find(slave_variable_key) != double_variable_map.end()) {
                    const auto& r_aux_var = double_variable_map.find(slave_variable_key)->second;
                    double& aux_value = rModelPart.pGetNode(r_dof_slave->Id())->FastGetSolutionStepValue(r_aux_var);
                    #pragma omp atomic
                    aux_value += slave_solution_vector[counter];
                } else if (components_variable_map.find(slave_variable_key) != components_variable_map.end()) {
                    const auto& r_aux_var = components_variable_map.find(slave_variable_key)->second;
                    double& aux_value = rModelPart.pGetNode(r_dof_slave->Id())->FastGetSolutionStepValue(r_aux_var);
                    #pragma omp atomic
                    aux_value += slave_solution_vector[counter];
                } else {
                    KRATOS_ERROR << "Dof variable is not defined" << std::endl;
                }

                ++counter;
            }
        }
    }

    KRATOS_CATCH("")
}

} // namespace ConstraintUtilities
} // namespace Kratos
