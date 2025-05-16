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
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class EliminationBuildDofArrayUtility
 * @ingroup KratosCore
 * @brief Utility class for handling the elimination build of DOFs
 * @details This static class collects all the methods required to
 * handle the DOFs when doing an elimination build.
 * @author Ruben Zorrilla
 */
class EliminationBuildDofArrayUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EliminationBuildDofArrayUtility
    KRATOS_CLASS_POINTER_DEFINITION(EliminationBuildDofArrayUtility);

    /// Definition of the size type
    using SizeType = std::size_t;

    /// Definition of the index type
    using IndexType = std::size_t;

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
    EliminationBuildDofArrayUtility() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetUpDofArray(
        const ModelPart& rModelPart,
        DofsArrayType& rDofArray,
#ifdef USE_LOCKS_IN_ASSEMBLY
        std::vector<omp_lock_t>& rLockArray,
#endif
        const unsigned int EchoLevel = 0,
        const bool CheckReactionDofs = false)
    {
        KRATOS_TRY;

        // Check that the provided DOF set is empty
        if (rDofArray.size() != 0) {
            KRATOS_WARNING_IF("EliminationBuildDofArrayUtility", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Provided DOF set is not empty. About to clear it." << std::endl;
            rDofArray.clear();
        }

        // Allocate auxiliary arrays
        DofsVectorType dof_list;
        SizeType n_threads = ParallelUtilities::GetNumThreads();
        std::vector<AuxiliaryDofsSetType> dofs_aux_list(n_threads);
        for (IndexType i = 0; i < n_threads; ++i) {
            dofs_aux_list[i].reserve(rModelPart.NumberOfElements());
        }

        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the DOFs" << std::endl;
        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", (EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Number of threads: " << ParallelUtilities::GetNumThreads() << std::endl;

        // Add the DOFs from the model part elements
        const auto& r_elements_array = rModelPart.Elements();
        const SizeType n_elements = rModelPart.NumberOfElements();
        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", EchoLevel > 2) << "Initializing elements loop" << std::endl;
        IndexPartition<SizeType>(n_elements).for_each(dof_list, [&](IndexType Index, DofsVectorType& tls_dof_list){
            // Get current element iterator
            const auto it_elem = r_elements_array.begin() + Index;

            // Gets list of Dof involved on every element
            it_elem->GetDofList(tls_dof_list, rModelPart.GetProcessInfo());
            dofs_aux_list[OpenMPUtils::ThisThread()].insert(tls_dof_list.begin(), tls_dof_list.end());
        });

        // Add the DOFs from the model part conditions
        const auto& r_conditions_array = rModelPart.Conditions();
        const SizeType n_conditions = rModelPart.NumberOfConditions();
        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", EchoLevel > 2) << "Initializing conditions loop" << std::endl;
        IndexPartition<SizeType>(n_conditions).for_each(dof_list, [&](IndexType Index, DofsVectorType& tls_dof_list){
            // Get current condition iterator
            const auto it_cond = r_conditions_array.begin() + Index;

            // Gets list of Dof involved on every condition
            it_cond->GetDofList(tls_dof_list, rModelPart.GetProcessInfo());
            dofs_aux_list[OpenMPUtils::ThisThread()].insert(tls_dof_list.begin(), tls_dof_list.end());
        });

        // Here we do a reduction in a tree so to have everything on thread 0
        SizeType old_max = n_threads;
        SizeType new_max = ceil(0.5*static_cast<double>(old_max));
        while (new_max >= 1 && new_max != old_max) {
            IndexPartition<SizeType>(new_max).for_each([&](SizeType Index){
                if (Index + new_max < old_max) {
                    dofs_aux_list[Index].insert(dofs_aux_list[Index + new_max].begin(), dofs_aux_list[Index + new_max].end());
                    dofs_aux_list[Index + new_max].clear();
                }
            });

            old_max = new_max;
            new_max = ceil(0.5*static_cast<double>(old_max));
        }

        // Fill and sort the provided DOF array from the auxiliary global DOFs set
        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", EchoLevel > 2) << "Initializing ordered array filling" << std::endl;
        rDofArray.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); ++it) {
            rDofArray.push_back(*it);
        }
        rDofArray.Sort();

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(rDofArray.size() == 0) << "No degrees of freedom in model part: " << rModelPart.FullName() << std::endl;
        KRATOS_INFO_IF("EliminationBuildDofArrayUtility", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the DOFs" << std::endl;

#ifdef KRATOS_DEBUG
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        // This is to be done only in debug mode
        if (CheckReactionDofs) {
            for (auto& r_dof : rDofArray) {
                KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
                    << "Reaction variable not set for the following: " << std::endl
                    << "- Node: "<< r_dof.Id() << std::endl
                    << "- DOF: " << r_dof << std::endl;
            }
        }
#endif

        KRATOS_CATCH("");
    }

#ifdef USE_LOCKS_IN_ASSEMBLY
    static void SetUpDofArray(
        const ModelPart& rModelPart,
        DofsArrayType& rDofArray,
        std::vector<omp_lock_t>& rLockArray,
        const unsigned int EchoLevel = 0,
        const bool CheckReactionDofs = false)
    {

        SetUpDofArray(rModelPart, rDofArray, EchoLevel, CheckReactionDofs);

        // If required by the architecture, reset the locks array
        if (rLockArray.size() != 0) {
            for (int i = 0; i < static_cast<int>(rLockArray.size()); i++) {
                omp_destroy_lock(&rLockArray[i]);
            }
        }

        rLockArray.resize(rDofArray.size());

        for (int i = 0; i < static_cast<int>(mLockArray.size()); i++) {
            omp_init_lock(&rLockArray[i]);
        }
    }
#endif

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
}; // Class EliminationBuildDofArrayUtility

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
