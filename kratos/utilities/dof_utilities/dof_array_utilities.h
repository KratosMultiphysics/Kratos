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
 * @class DofArrayUtilities
 * @ingroup KratosCore
 * @brief Utility class for handling the build of DOFs
 * @details This static class collects some functions to handle the DOFs.
 * @author Ruben Zorrilla
 */
class DofArrayUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DofArrayUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DofArrayUtilities);

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

    /// Auxiliary slave to master(s) map for slave DOFs
    using SlaveToMasterDofsMap = std::unordered_map<typename Node::DofType::Pointer, DofsVectorType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DofArrayUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    static void SetUpDofArray(
        const ModelPart& rModelPart,
        DofsArrayType& rDofArray,
        const unsigned int EchoLevel = 0)
    {
        KRATOS_TRY;

        // Check that the provided DOF set is empty
        if (rDofArray.size() != 0) {
            KRATOS_WARNING_IF("DofArrayUtilities", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Provided DOF set is not empty. About to clear it." << std::endl;
            rDofArray.clear();
        }

        // Allocate auxiliary arrays
        DofsVectorType dof_list;
        SizeType n_threads = ParallelUtilities::GetNumThreads();
        std::vector<AuxiliaryDofsSetType> dofs_aux_list(n_threads);
        for (IndexType i = 0; i < n_threads; ++i) {
            dofs_aux_list[i].reserve(rModelPart.NumberOfElements());
        }

        KRATOS_INFO_IF("DofArrayUtilities", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the DOFs" << std::endl;
        KRATOS_INFO_IF("DofArrayUtilities", (EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Number of threads: " << ParallelUtilities::GetNumThreads() << std::endl;

        // Add the DOFs from the model part elements
        const auto& r_elements_array = rModelPart.Elements();
        const SizeType n_elements = rModelPart.NumberOfElements();
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing elements loop" << std::endl;
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
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing conditions loop" << std::endl;
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
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing ordered array filling" << std::endl;
        rDofArray.reserve(dofs_aux_list[0].size());
        for (auto it = dofs_aux_list[0].begin(); it != dofs_aux_list[0].end(); ++it) {
            rDofArray.push_back(*it);
        }
        rDofArray.Sort();

        // Throw an exception if there are no DOFs involved in the analysis
        KRATOS_ERROR_IF(rDofArray.size() == 0) << "No degrees of freedom in model part: " << rModelPart.FullName() << std::endl;
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the DOFs" << std::endl;

        KRATOS_CATCH("");
    }

    static void SetUpEffectiveDofArray(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofArray,
        DofsArrayType& rEffectiveDofArray,
        SlaveToMasterDofsMap& rSlaveToMasterDofsMap,
        const unsigned int EchoLevel = 0)
    {
        KRATOS_TRY;

        // Clear the provided containers
        KRATOS_WARNING_IF("Builder", !rEffectiveDofArray.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
        rEffectiveDofArray.clear();
        KRATOS_WARNING_IF("Builder", !rSlaveToMasterDofsMap.empty()) << "Provided slave to master DOFs map is not empty. About to clear it." << std::endl;
        rSlaveToMasterDofsMap.clear();

        //FIXME: In here we should check if there are ACTIVE constraints among all processes
        // Check if there are constraints to build the effective DOFs map and the corresponding arrays
        const std::size_t n_constraints = rModelPart.NumberOfMasterSlaveConstraints();
        if (n_constraints) {
            // Auxiliary set to store the unordered effective DOFs (masters from constraints and standard ones)
            std::unordered_set<typename Node::DofType::Pointer> effective_dofs_set;

            // Get the master / slave DOFs from the constraints
            KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Getting master and slave DOFs from constraints." << std::endl;
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            for (IndexType i_const = 0; i_const < n_constraints; ++i_const) {
                // Get current constraint master and slave DOFs
                auto it_const = it_const_begin + i_const;
                const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
                const auto& r_master_dofs = it_const->GetMasterDofsVector();

                // Add the slave DOFs to the slave map
                for (auto& rp_slave : r_slave_dofs) {
                    rSlaveToMasterDofsMap.insert(std::make_pair(rp_slave, r_master_dofs));
                }

                // Add the master DOFs to the effective DOFs set
                // Note that we initialize the system ids to zero as these will be overwritten later
                for (auto& rp_master : r_master_dofs) {
                    effective_dofs_set.insert(rp_master);
                }
            }

            // Loop the elements and conditions DOFs container to get the DOFs that are not slave
            KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Getting non-slave DOFs from elements and conditions." << std::endl;
            for (IndexType i_dof = 0; i_dof < rDofArray.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofArray.ptr_begin() + i_dof);

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered in the resolution of the system
                // Note that this includes masters DOFs or and standard DOFs (those not involved in any constraint)
                if (rSlaveToMasterDofsMap.find(p_dof) == rSlaveToMasterDofsMap.end()) {
                    // Add current DOF to the effective DOFs set (note that the std::unordered_set guarantees uniqueness)
                    effective_dofs_set.insert(p_dof);
                }
            }

            // Sort the effective DOFs before setting the equation ids
            // Note that we dereference the DOF pointers in order to use the greater operator from dof.h
            KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Sorting the effective DOFs auxiliary container." << std::endl;
            std::vector<typename Node::DofType::Pointer> ordered_eff_dofs_vector(effective_dofs_set.begin(), effective_dofs_set.end());
            std::sort(
                ordered_eff_dofs_vector.begin(),
                ordered_eff_dofs_vector.end(),
                [](const typename Node::DofType::Pointer& pA, const typename Node::DofType::Pointer& pB){return *pA > *pB;});

            // Fill the effective DOFs PVS with the sorted effective DOFs container
            KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Creating the effective DOFs array from the auxiliary sorted container." << std::endl;
            rEffectiveDofArray = DofsArrayType(ordered_eff_dofs_vector);
        } else {
            // If there are no constraints the effective DOF set is the standard one
            KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "There are no constraints. The effective DOF array matches the standard one." << std::endl;
            rEffectiveDofArray = rDofArray;
        }

        KRATOS_CATCH("");
    }

    static void SetUpDofIds(const DofsArrayType& rDofArray)
    {
        // Set up the DOFs equation global ids
        IndexPartition<IndexType>(rDofArray.size()).for_each([&](IndexType Index) {
            auto it_dof = rDofArray.begin() + Index;
            it_dof->SetEquationId(Index);
        });
    }

    static void SetUpEffectiveDofIds(
        const DofsArrayType& rDofArray,
        DofsArrayType& rEffectiveDofArray)
    {
        // Check if the effective and "standard" containers are the same
        // We do it with the addresses to avoid checking the content (i.e., each DOF one-by-one)
        if (&rEffectiveDofArray == &rDofArray) {
            // Set the DOFs' effective equation global ids to match the standard ones
            IndexPartition<IndexType>(rEffectiveDofArray.size()).for_each([&](IndexType Index) {
                auto it_dof = rEffectiveDofArray.begin() + Index;
                it_dof->SetEffectiveEquationId(it_dof->EquationId());
            });
        } else {
            // Set the effective DOFs equation ids
            // Note that in here we assume the effective DOFs to be already sorted
            IndexPartition<IndexType>(rEffectiveDofArray.size()).for_each([&](IndexType Index) {
                auto it_dof = rEffectiveDofArray.begin() + Index;
                it_dof->SetEffectiveEquationId(Index);
            });
        }
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
}; // Class DofArrayUtilities

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
