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

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "dof_array_utilities.h"

namespace Kratos
{

void DofArrayUtilities::SetUpDofArray(
    const ModelPart& rModelPart,
    DofsArrayType& rDofArray,
    const unsigned int EchoLevel)
{
    KRATOS_TRY;

    // Check that the provided DOF set is empty
    if (rDofArray.size() != 0) {
        KRATOS_WARNING_IF("DofArrayUtilities", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Provided DOF set is not empty. About to clear it." << std::endl;
        rDofArray.clear();
    }

    // Allocate auxiliary arrays
    DofsVectorType dof_list;
    unsigned int assumed_dofs_per_element = 20;
    std::size_t n_threads = ParallelUtilities::GetNumThreads();
    std::vector<AuxiliaryDofsSetType> dofs_aux_list(n_threads);
    for (IndexType i = 0; i < n_threads; ++i) {
        dofs_aux_list[i].reserve(assumed_dofs_per_element*rModelPart.NumberOfElements());
    }

    KRATOS_INFO_IF("DofArrayUtilities", (EchoLevel > 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Setting up the DOFs" << std::endl;
    KRATOS_INFO_IF("DofArrayUtilities", (EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Number of threads: " << ParallelUtilities::GetNumThreads() << std::endl;

    // Add the DOFs from the model part elements
    const auto& r_elements_array = rModelPart.Elements();
    const std::size_t n_elements = rModelPart.NumberOfElements();
    KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing elements loop" << std::endl;
    IndexPartition<std::size_t>(n_elements).for_each(dof_list, [&](IndexType Index, DofsVectorType& tls_dof_list){
        // Get current element iterator
        const auto it_elem = r_elements_array.begin() + Index;

        // Gets list of Dof involved on every element
        it_elem->GetDofList(tls_dof_list, rModelPart.GetProcessInfo());
        dofs_aux_list[OpenMPUtils::ThisThread()].insert(tls_dof_list.begin(), tls_dof_list.end());
    });

    // Add the DOFs from the model part conditions
    const auto& r_conditions_array = rModelPart.Conditions();
    const std::size_t n_conditions = rModelPart.NumberOfConditions();
    KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing conditions loop" << std::endl;
    IndexPartition<std::size_t>(n_conditions).for_each(dof_list, [&](IndexType Index, DofsVectorType& tls_dof_list){
        // Get current condition iterator
        const auto it_cond = r_conditions_array.begin() + Index;

        // Gets list of Dof involved on every condition
        it_cond->GetDofList(tls_dof_list, rModelPart.GetProcessInfo());
        dofs_aux_list[OpenMPUtils::ThisThread()].insert(tls_dof_list.begin(), tls_dof_list.end());
    });

    // Here we do a reduction in a tree so to have everything on thread 0
    std::size_t old_max = n_threads;
    std::size_t new_max = std::ceil(0.5*static_cast<double>(old_max));
    while (new_max >= 1 && new_max != old_max) {
        IndexPartition<std::size_t>(new_max).for_each([&](std::size_t Index){
            if (Index + new_max < old_max) {
                dofs_aux_list[Index].insert(dofs_aux_list[Index + new_max].begin(), dofs_aux_list[Index + new_max].end());
            }
        });

        old_max = new_max;
        new_max = std::ceil(0.5*static_cast<double>(old_max));
    }

    // Fill and sort the provided DOF array from the reduced (i.e., thread 0) auxiliary global DOFs set
    KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "Initializing ordered array filling" << std::endl;
    rDofArray.insert(dofs_aux_list[0].begin(), dofs_aux_list[0].end());

    // Throw an exception if there are no DOFs involved in the analysis
    KRATOS_ERROR_IF(rDofArray.size() == 0) << "No degrees of freedom in model part: " << rModelPart.FullName() << std::endl;
    KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Finished setting up the DOFs" << std::endl;

    KRATOS_CATCH("");
}

void DofArrayUtilities::SetUpEffectiveDofArray(
    const ModelPart &rModelPart,
    const DofsArrayType &rDofArray,
    DofsArrayType &rEffectiveDofArray,
    const unsigned int EchoLevel)
{
    KRATOS_TRY;

    //TODO: MPI implementation
    KRATOS_ERROR_IF(rModelPart.GetCommunicator().IsDistributed()) << "MPI parallelism is not available yet." << std::endl;

    // Clear the provided containers
    KRATOS_WARNING_IF("Builder", !rEffectiveDofArray.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
    rEffectiveDofArray.clear();

    //FIXME: In here we should check if there are ACTIVE constraints among all processes
    // Check if there are constraints to build the effective DOFs map and the corresponding arrays
    const std::size_t n_constraints = rModelPart.NumberOfMasterSlaveConstraints();
    if (n_constraints) {
        // Auxiliary set to store the unordered effective DOFs (masters from constraints and standard ones)
        std::unordered_set<typename Node::DofType::Pointer> effective_dofs_set;

        // Get the master / slave DOFs from the constraints
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2 && rModelPart.GetCommunicator().MyPID() == 0) << "Getting master and slave DOFs from constraints." << std::endl;
        std::unordered_map<Node::DofType::Pointer, ModelPart::DofsVectorType> slave_to_master_dofs_map;
        const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
        for (IndexType i_const = 0; i_const < n_constraints; ++i_const) {
            // Get current constraint master and slave DOFs
            auto it_const = it_const_begin + i_const;
            const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
            const auto& r_master_dofs = it_const->GetMasterDofsVector();

            // Add the slave DOFs to the slave map
            for (auto& rp_slave : r_slave_dofs) {
                slave_to_master_dofs_map.insert(std::make_pair(rp_slave, r_master_dofs));
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
            if (slave_to_master_dofs_map.find(p_dof) == slave_to_master_dofs_map.end()) {
                // Add current DOF to the effective DOFs set (note that the std::unordered_set guarantees uniqueness)
                effective_dofs_set.insert(p_dof);
            }
        }

        // Fill the effective DOFs PVS with the effective DOFs container
        rEffectiveDofArray.insert(effective_dofs_set.begin(), effective_dofs_set.end());
    } else {
        // If there are no constraints the effective DOF set is the standard one
        KRATOS_INFO_IF("DofArrayUtilities", EchoLevel > 2) << "There are no constraints. The effective DOF array matches the standard one." << std::endl;
        rEffectiveDofArray = rDofArray;
    }

    KRATOS_CATCH("");
}

void DofArrayUtilities::SetDofEquationIds(const DofsArrayType &rDofArray)
{
    // Set up the DOFs equation global ids
    IndexPartition<IndexType>(rDofArray.size()).for_each([&](IndexType Index) {
        auto it_dof = rDofArray.begin() + Index;
        it_dof->SetEquationId(Index);
    });
}

void DofArrayUtilities::SetEffectiveDofEquationIds(
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
        // Initialize all DOFs effective equation ids to the maximum allowable value
        // Note that this makes possible to distingish the effective DOFs from the non-effective ones
        IndexPartition<IndexType>(rDofArray.size()).for_each([&](IndexType Index) {
            auto it_dof = rDofArray.begin() + Index;
            it_dof->SetEffectiveEquationId(std::numeric_limits<typename Node::DofType::EquationIdType>::max());
        });

        // Set the effective DOFs equation ids
        // Note that in here we assume the effective DOFs to be already sorted
        IndexPartition<IndexType>(rEffectiveDofArray.size()).for_each([&](IndexType Index) {
            auto it_dof = rEffectiveDofArray.begin() + Index;
            it_dof->SetEffectiveEquationId(Index);
        });
    }
}

}  // namespace Kratos.
