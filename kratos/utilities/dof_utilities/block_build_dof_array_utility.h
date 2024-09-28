//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class BlockBuildDofArrayUtility
 * @ingroup KratosCore
 * @brief Utility class for handling the block build of DOFs
 * @details This static class collects all the methods required to
 * handle the DOFs when doing a block build.
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) BlockBuildDofArrayUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockBuildDofArrayUtility
    KRATOS_CLASS_POINTER_DEFINITION(BlockBuildDofArrayUtility);

    /// DOFs array type definition from ModelPart (note that this is PointerVectorSet<DofType>)
    using DofsArrayType = ModelPart::DofsArrayType;

    /// DOFs vector type definition from ModelPart (note that this is std::vector<DofType::Pointer>)s
    using DofsVectorType = ModelPart::DofsVectorType;

    /// Auxilary DOF set type definition
    using AuxiliaryDofsSetType = std::unordered_set<Node::DofType::Pointer, DofPointerHasher>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BlockBuildDofArrayUtility() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetUpDofArray(
        const ModelPart& rModelPart,
        DofsArrayType& rDofArray,
        const unsigned int EchoLevel = 0,
        const bool CheckReactionDofs = false)
    {
        KRATOS_TRY;

        // Check that the provided DOF set is empty
        if (rDofArray.size() != 0) {
            KRATOS_WARNING_IF("BlockBuildDofArrayUtility", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Provided DOF set is not empty. About to clear it." << std::endl;
            rDofArray.clear();
        }

        // Allocate auxiliary arrays
        DofsVectorType dof_list;
        DofsVectorType second_dof_list; // NOTE: The second dof list is only used on constraints to include master/slave relations

        AuxiliaryDofsSetType dof_global_set;
        dof_global_set.reserve(rModelPart.NumberOfElements() * 20);

        KRATOS_INFO_IF("BlockBuildDofArrayUtility", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the DOFs" << std::endl;
        KRATOS_INFO_IF("BlockBuildDofArrayUtility", (EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Number of threads: " << ParallelUtilities::GetNumThreads() << std::endl;

        // Fill the global DOF set array
        #pragma omp parallel firstprivate(dof_list, second_dof_list)
        {
            const auto& r_current_process_info = rModelPart.GetProcessInfo();

            // We create the temporal set in current thread and we reserve some space on it
            AuxiliaryDofsSetType dofs_tmp_set;
            dofs_tmp_set.reserve(20000);

            // Add the DOFs from the model part elements
            const auto& r_elements_array = rModelPart.Elements();
            const int number_of_elements = static_cast<int>(r_elements_array.size());
            KRATOS_INFO_IF("BlockBuildDofArrayUtility", EchoLevel > 2) << "Initializing elements loop" << std::endl;
            #pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < number_of_elements; ++i) {
                // Get current element iterator
                const auto it_elem = r_elements_array.cbegin() + i;

                // Gets list of DOF involved on every element
                it_elem->GetDofList(dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Add the DOFs from the model part conditions
            const auto& r_conditions_array = rModelPart.Conditions();
            const int number_of_conditions = static_cast<int>(r_conditions_array.size());
            KRATOS_INFO_IF("BlockBuildDofArrayUtility", EchoLevel > 2) << "Initializing conditions loop" << std::endl;
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_conditions; ++i) {
                // Get current condition iterator
                const auto it_cond = r_conditions_array.cbegin() + i;

                // Gets list of DOF involved on every condition
                it_cond->GetDofList(dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
            }

            // Add the DOFs from the model part constraints
            const auto& r_constraints_array = rModelPart.MasterSlaveConstraints();
            const int number_of_constraints = static_cast<int>(r_constraints_array.size());
            KRATOS_INFO_IF("BlockBuildDofArrayUtility", EchoLevel > 2) << "Initializing constraints loop" << std::endl;
            #pragma omp for  schedule(guided, 512) nowait
            for (int i = 0; i < number_of_constraints; ++i) {
                // Get current constraint iterator
                const auto it_const = r_constraints_array.cbegin() + i;

                // Gets list of DOF involved on every constraint
                it_const->GetDofList(dof_list, second_dof_list, r_current_process_info);
                dofs_tmp_set.insert(dof_list.begin(), dof_list.end());
                dofs_tmp_set.insert(second_dof_list.begin(), second_dof_list.end());
            }

            // Merge all the sets in one thread
            #pragma omp critical
            {
                dof_global_set.insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
            }
        }

        // Fill and sort the provided DOF array from the auxiliary global DOFs set
        KRATOS_INFO_IF("BlockBuildDofArrayUtility", EchoLevel > 2) << "Initializing ordered array filling" << std::endl;
        rDofArray.reserve(dof_global_set.size());
        for (auto it = dof_global_set.begin(); it!= dof_global_set.end(); it++) {
            rDofArray.push_back(*it);
        }
        rDofArray.Sort();

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(rDofArray.size() == 0) << "No degrees of freedom in model part: " << rModelPart.FullName() << std::endl;
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolver", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the DOFs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is to be done only in debug mode
        if (CheckReactionDofs) {
            for (auto it_dof : rDofArray) {
                KRATOS_ERROR_IF_NOT(it_dof->HasReaction())
                    << "Reaction variable not set for the following: " << std::endl
                    << "- Node : "<< it_dof->Id() << std::endl;
                    << "- DOF: " << (*it_dof) << std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
}; // Class BlockBuildDofArrayUtility

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
