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
 * @class BlockBuilder
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class BlockBuilder : public Builder<TThreadLocalStorage, TSparseMatrixType, TSystemVectorType, TSparseGraphType>
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of BlockBuilder
    KRATOS_CLASS_POINTER_DEFINITION(BlockBuilder);

    /// Base builder type definition
    using BaseType = Builder<TThreadLocalStorage, TSparseMatrixType, TSystemVectorType, TSparseGraphType>;

    /// Data type definition from sparse matrix
    using DataType = typename TSparseMatrixType::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename TSparseMatrixType::IndexType;

    /// DOF type definition
    using DofType = typename BaseType::DofType;

    /// DOF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// DOF pointer vector type definition
    using DofPointerVectorType = typename BaseType::DofPointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BlockBuilder() = delete;

    /// Constructor with model part
    BlockBuilder(
        const ModelPart &rModelPart,
        Parameters Settings = Parameters(R"({})"))
        : BaseType(rModelPart, Settings)
    {
        Parameters default_parameters( R"({
            "name" : "block_builder",
            "scaling_type" : "max_diagonal",
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);
    }

    virtual ~BlockBuilder() = default;

    ///@}
    ///@name Operations
    ///@{

    void ConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Clear the provided effective DOFs map
        KRATOS_WARNING_IF("BlockBuilder", !rEffectiveDofSet.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
        rEffectiveDofSet.clear();

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
            rEffectiveDofSet = DofsArrayType(ordered_eff_dofs_vector);

            // Set the effective DOFs equation ids based on the sorted list
            IndexType aux_dof_id = 0;
            for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
                auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);
                p_dof->SetEffectiveEquationId(aux_dof_id);
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
                        auto eff_dof_find = rEffectiveDofSet.find(*rp_master);
                        KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                        constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
                    }
                } else { // Effective DOF
                    auto eff_dof_find = rEffectiveDofSet.find(*p_dof);
                    KRATOS_ERROR_IF(eff_dof_find == rEffectiveDofSet.end()) << "Effective DOF cannot be find." << std::endl;
                    // mMasterIds.push_back(eff_dof_find->second);
                    constraints_sparse_graph.AddEntry(i_dof_eq_id, eff_dof_find->EffectiveEquationId());
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
            auto p_aux_q = Kratos::make_shared<TSystemVectorType>(this->GetEquationSystemSize());
            rLinearSystemContainer.pConstraintsQ.swap(p_aux_q);

            auto p_aux_T = Kratos::make_shared<TSparseMatrixType>(constraints_sparse_graph);
            rLinearSystemContainer.pConstraintsT.swap(p_aux_T);
        } else {
            rEffectiveDofSet = rDofSet; // If there are no constraints the effective DOF set is the standard one
        }
    }

    //FIXME: Do the RHS-only version
    void ApplyLinearSystemConstraints(
        const DofsArrayType& rEffectiveDofArray,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Calculate the effective LHS, RHS and solution vector
        ApplyBlockBuildMasterSlaveConstraints(rLinearSystemContainer);

        // Apply the Dirichlet BCs in a block way by leveraging the CSR matrix implementation
        auto& r_eff_rhs = *(rLinearSystemContainer.pEffectiveRhs);
        auto& r_eff_lhs = *(rLinearSystemContainer.pEffectiveLhs);
        ApplyBlockBuildDirichletConditions(rEffectiveDofArray, r_eff_lhs, r_eff_rhs);
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void ApplyBlockBuildMasterSlaveConstraints(LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer)
    {
        const auto& r_model_part = this->GetModelPart();
        const std::size_t n_constraints = r_model_part.NumberOfMasterSlaveConstraints();
        if (n_constraints) { //FIXME: In here we should check the number of active constraints
            // Get effective arrays
            auto p_eff_dx = rLinearSystemContainer.pEffectiveDx;
            auto p_eff_rhs = rLinearSystemContainer.pEffectiveRhs;

            // Initialize the effective RHS
            KRATOS_ERROR_IF(p_eff_rhs == nullptr) << "Effective RHS vector has not been initialized yet." << std::endl;
            p_eff_rhs->SetValue(0.0);

            // Initialize the effective solution vector
            KRATOS_ERROR_IF(p_eff_rhs == nullptr) << "Effective solution increment vector has not been initialized yet." << std::endl;
            p_eff_dx->SetValue(0.0);

            // Apply constraints to RHS
            auto p_rhs = rLinearSystemContainer.pRhs;
            auto p_constraints_T = rLinearSystemContainer.pConstraintsT;
            p_constraints_T->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            auto p_lhs = rLinearSystemContainer.pLhs;
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *p_constraints_T);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*p_constraints_T);
            rLinearSystemContainer.pEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);

            // // Compute the scale factor value
            // //TODO: think on how to make this user-definable
            // const double scale_factor = rpEffectiveLhs->NormDiagonal();

            // // Apply diagonal values on slave DOFs
            // IndexPartition<IndexType>(mSlaveIds.size()).for_each([&](IndexType Index){
            //     const IndexType slave_eq_id = mSlaveIds[Index];
            //     if (mInactiveSlaveDofs.find(slave_eq_id) == mInactiveSlaveDofs.end()) {
            //         (*rpEffectiveLhs)(slave_eq_id, slave_eq_id) = scale_factor;
            //         (*rpEffectiveRhs)[slave_eq_id] = 0.0;
            //     }
            // });
        } else {
            // If there are no constraints the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rLinearSystemContainer.pEffectiveLhs = rLinearSystemContainer.pLhs;
            rLinearSystemContainer.pEffectiveRhs = rLinearSystemContainer.pRhs;
            rLinearSystemContainer.pEffectiveDx = rLinearSystemContainer.pDx;
        }
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        TSparseMatrixType& rLHS,
        TSystemVectorType& rRHS) const
    {
        // Set the free DOFs vector (0 means fixed / 1 means free)
        // Note that we initialize to 1 so we start assuming all free
        // Also note that the type is uint_8 for the sake of efficiency
        const std::size_t system_size = rLHS.size1();
        std::vector<uint8_t> free_dofs_vector(system_size, 1);

        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_begin = rDofArray.begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            const auto p_dof = dof_begin + Index;
            if (p_dof->IsFixed()) {
                free_dofs_vector[p_dof->EffectiveEquationId()] = 0;
            }
        });

        //TODO: Implement this in the CSR matrix or here?
        // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
        // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        // Get the diagonal scaling factor
        const double diagonal_value = this->GetDiagonalScalingFactor(rLHS);

        // Apply the free DOFs (i.e., fixity) vector to the system arrays
        rLHS.ApplyHomogeneousDirichlet(free_dofs_vector, diagonal_value, rRHS);
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        TSystemVectorType& rRHS) const
    {
        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_begin = rDofArray.begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            auto p_dof = dof_begin + Index;
            if (p_dof->IsFixed()) {
                rRHS[p_dof->EffectiveEquationId()] = 0;
            }
        });
    }

    ///@}
}; // Class BlockBuilder

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
