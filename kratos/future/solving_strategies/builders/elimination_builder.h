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
#include "builder.h"
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "containers/sparse_contiguous_row_graph.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/amgcl_csr_conversion_utilities.h"
#include "utilities/amgcl_csr_spmm_utilities.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class EliminationBuilder
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSparseVectorType, class TSparseGraphType>
class EliminationBuilder : public Builder<TThreadLocalStorage, TSparseMatrixType, TSparseVectorType, TSparseGraphType>
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of EliminationBuilder
    KRATOS_CLASS_POINTER_DEFINITION(EliminationBuilder);

    /// Base builder type definition
    using BaseType = Builder<TThreadLocalStorage, TSparseMatrixType, TSparseVectorType, TSparseGraphType>;

    /// Data type definition from sparse matrix
    using DataType = typename TSparseMatrixType::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename TSparseMatrixType::IndexType;

    /// DOF type definition
    using DofType = typename BaseType::DofType;

    /// DOF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Effective DOFs map type definition
    using EffectiveDofsMapType = typename BaseType::EffectiveDofsMapType;

    /// DOF pointer vector type definition
    using DofPointerVectorType = typename BaseType::DofPointerVectorType;

    /// Function type for elements assembly
    using ElementAssemblyFunctionType = typename BaseType::ElementAssemblyFunctionType;

    /// Function type for conditions assembly
    using ConditionAssemblyFunctionType = typename BaseType::ConditionAssemblyFunctionType;

    /// Function type for constraints assembly
    using ConstraintAssemblyFunctionType = typename BaseType::ConstraintAssemblyFunctionType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EliminationBuilder() = delete;

    /// Constructor with model part
    EliminationBuilder(
        const ModelPart &rModelPart,
        Parameters Settings = Parameters(R"({})"))
        : BaseType(rModelPart, Settings)
    {
        Parameters default_parameters( R"({
            "name" : "elimination_builder",
            "scaling_type" : "max_diagonal",
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);
    }

    virtual ~EliminationBuilder() = default;

    ///@}
    ///@name Operations
    ///@{

    void ConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector) override
    {
        // Clear the provided effective DOFs map
        KRATOS_WARNING_IF("EliminationBuilder", !rEffectiveDofSet.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
        KRATOS_WARNING_IF("EliminationBuilder", !rEffectiveDofIdMap.empty()) << "Provided effective DOFs ids map is not empty. About to clear it." << std::endl;
        rEffectiveDofSet.clear();
        rEffectiveDofIdMap.clear();

        //FIXME: Do the IsActiveConstraints in here and set a flag that stays "forever"

        // Check if there are constraints to build the effective DOFs map and the corresponding arrays
        const std::size_t n_constraints = rModelPart.NumberOfMasterSlaveConstraints();
        if (n_constraints) {

            // Auxiliary set to store the unordered effective DOFs (masters from constraints and standard ones)
            std::unordered_set<typename DofType::Pointer> effective_dofs_set;

            // Get the master / slave DOFs from the constraints
            std::unordered_map<typename DofType::Pointer, DofPointerVectorType> constraints_slave_dofs;
            const auto it_const_begin = rModelPart.MasterSlaveConstraints().begin();
            for (IndexType i_const = 0; i_const < n_constraints; ++i_const) {
                // Get current constraint master and slave DOFs
                auto it_const = it_const_begin + i_const;
                const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
                const auto& r_master_dofs = it_const->GetMasterDofsVector();

                // Add the slave DOFs to the slave map
                for (auto& rp_slave : r_slave_dofs) {
                    constraints_slave_dofs.insert(std::make_pair(rp_slave, r_master_dofs));
                }

                // Add the master DOFs to the effective DOFs set
                // Note that we initialize the system ids to zero as these will be overwritten later
                for (auto& rp_master : r_master_dofs) {
                    effective_dofs_set.insert(rp_master);
                }
            }

            // Loop the elements and conditions DOFs container to get the DOFs that are not slave
            for (IndexType i_dof = 0; i_dof < rDofSet.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofSet.ptr_begin() + i_dof);

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered in the resolution of the system
                // Note that this includes masters DOFs or and standard DOFs (those not involved in any constraint)
                if (constraints_slave_dofs.find(p_dof) == constraints_slave_dofs.end()) {
                    // Add current DOF to the effective DOFs set (note that the std::unordered_set guarantees uniqueness)
                    effective_dofs_set.insert(p_dof);
                }
            }

            // Sort the effective DOFs before setting the equation ids
            // Note that we dereference the DOF pointers in order to use the greater operator from dof.h
            std::vector<typename DofType::Pointer> ordered_eff_dofs_vector(effective_dofs_set.begin(), effective_dofs_set.end());
            std::sort(
                ordered_eff_dofs_vector.begin(),
                ordered_eff_dofs_vector.end(),
                [](const typename DofType::Pointer& pA, const typename DofType::Pointer& pB){return *pA > *pB;});

            // Fill the effective DOFs PVS with the sorted effective DOFs container
            rEffectiveDofSet = std::move(DofsArrayType(ordered_eff_dofs_vector));

            // Set the effective DOFs equation ids based on the sorted list
            rEffectiveDofIdMap.reserve(rEffectiveDofSet.size());
            IndexType aux_dof_id = 0;
            for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
                auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);
                rEffectiveDofIdMap.insert(std::make_pair(p_dof, aux_dof_id));
                ++aux_dof_id;
            }

            // Clear the equation ids vectors
            // mSlaveIds.clear();
            // mMasterIds.clear();

            // Set up constraints matrix sparse graph (note that mEquationSystemSize is the DOF set size in the block build)
            KRATOS_ERROR_IF(this->GetEquationSystemSize() == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
            TSparseGraphType constraints_sparse_graph(this->GetEquationSystemSize());

            // Loop the elements and conditions DOFs container to add the slave entries to the graph
            for (IndexType i_dof = 0; i_dof < rDofSet.size(); ++i_dof) {
                // Get current DOF
                auto p_dof = *(rDofSet.ptr_begin() + i_dof);
                const IndexType i_dof_eq_id = p_dof->EquationId();

                // Check if current DOF is slave by checking the slaves DOFs map
                // If not present in the slaves DOFs map it should be considered a "master" DOF
                // Note that here "master" means an actual masters DOF or a DOF that do not involve any constraint
                auto i_dof_slave_find = constraints_slave_dofs.find(p_dof);
                if (i_dof_slave_find != constraints_slave_dofs.end()) { // Slave DOF
                    // // Add current slave DOF to slave equation ids list
                    // mSlaveIds.push_back(i_dof_eq_id);

                    // Add current slave DOF connectivities to the constraints sparse graph
                    // The slave rows eq ids come from the system ones while the column master ones are the above defined
                    for (auto& rp_master : i_dof_slave_find->second) {
                        auto eff_dof_find = rEffectiveDofIdMap.find(rp_master);
                        KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofIdMap.end()) << "Effective DOF cannot be find." << std::endl;
                        constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->second);
                    }
                } else { // Effective DOF
                    auto eff_dof_find = rEffectiveDofIdMap.find(p_dof);
                    KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofIdMap.end()) << "Effective DOF cannot be find." << std::endl;
                    // mMasterIds.push_back(eff_dof_find->second);
                    constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->second);
                }
            }

            // // Loop the effective DOFs container to add the remaining diagonal entries to the graph
            // for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
            //     // Get current effective DOF
            //     auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);
            //     auto eff_dof_find = rEffectiveDofIdMap.find(p_dof);
            //     KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofIdMap.end()) << "Effective DOF cannot be find." << std::endl;
            //     std::cout << "Effective DOF: " << eff_dof_find->second << " - " << eff_dof_find->second << std::endl;
            //     constraints_sparse_graph.AddEntry(eff_dof_find->second, eff_dof_find->second);
            // }

            // Allocate the constraints arrays (note that we are using the move assignment operator in here)
            rConstraintsConstantVector = std::move(TSparseVectorType(this->GetEquationSystemSize()));
            rConstraintsRelationMatrix = std::move(TSparseMatrixType(constraints_sparse_graph));

        } else {
            rEffectiveDofSet = rDofSet; // If there are no constraints the effective DOF set is the standard one
            rEffectiveDofIdMap = EffectiveDofsMapType(); // Create an empty master ids map as the standard DOF equation ids can be used
        }
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
}; // Class EliminationBuilder

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
