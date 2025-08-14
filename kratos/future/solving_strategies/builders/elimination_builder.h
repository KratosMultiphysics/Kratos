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
#include "utilities/reduction_utilities.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class EliminationBuilder
 * @ingroup KratosCore
 * @brief Utility class for handling the elimination build
 * @details This helper class extends the base Builder class to consider an elimination type build
 * The elimination type build removes the Dirichlet DOFs from the the linear system of equations to
 * be solved. This is achieved "a la master-slave" by creating an extra constraints relation matrix.
 * This Dirichlet constraints relation matrix is applied to the linear system of equations, resulting
 * in the effective removal of the fixed DOFs. Note that the constraints constant vector associated
 * to the Dirichlet constraints is never build as we always solve for the solution increment and the
 * Dirichlet values are inherently taken into account in the residual database as we store them in
 * the nodal database.
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class EliminationBuilder : public Builder<TThreadLocalStorage, TSparseMatrixType, TSystemVectorType, TSparseGraphType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EliminationBuilder
    KRATOS_CLASS_POINTER_DEFINITION(EliminationBuilder);

    /// Base builder type definition
    using BaseType = Builder<TThreadLocalStorage, TSparseMatrixType, TSystemVectorType, TSparseGraphType>;

    /// Index type definition from sparse matrix
    using IndexType = typename TSparseMatrixType::IndexType;

    /// DOF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

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
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);
    }

    virtual ~EliminationBuilder() = default;

    ///@}
    ///@name Operations
    ///@{

    void AllocateLinearSystemArrays(
        const typename DofsArrayType::Pointer pDofSet,
        const typename DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer) override
    {
        // // Set up the system sparse matrix graph (note that the sparse graph will be destroyed when leaving this scope)
        BuiltinTimer sparse_matrix_graph_time;
        TSparseGraphType sparse_matrix_graph(pDofSet->size());
        this->SetUpSparseMatrixGraph(sparse_matrix_graph);
        KRATOS_INFO_IF("BlockBuilder", this->GetEchoLevel() > 0) << "Set up sparse matrix graph time: " << sparse_matrix_graph_time << std::endl;

        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_dx = Kratos::make_shared<TSystemVectorType>(sparse_matrix_graph);
        rLinearSystemContainer.pDx.swap(p_dx);

        auto p_rhs = Kratos::make_shared<TSystemVectorType>(sparse_matrix_graph);
        rLinearSystemContainer.pRhs.swap(p_rhs);

        auto p_lhs = Kratos::make_shared<TSparseMatrixType>(sparse_matrix_graph);
        rLinearSystemContainer.pLhs.swap(p_lhs);

        // Get the number of free DOFs to allocate the effective system arrays
        const auto dof_begin = pEffectiveDofSet->begin();
        const std::size_t n_free_dofs = IndexPartition<std::size_t>(pEffectiveDofSet->size()).for_each<SumReduction<std::size_t>>([&](IndexType Index) {
            auto p_dof = dof_begin + Index;
            return p_dof->IsFixed() ? 0 : 1;
        });

        // Allocate the effective arrays according to the number of free effective DOFs
        auto p_eff_lhs = Kratos::make_shared<TSparseMatrixType>();
        rLinearSystemContainer.pEffectiveLhs.swap(p_eff_lhs);

        auto p_eff_rhs = Kratos::make_shared<TSystemVectorType>(n_free_dofs);
        rLinearSystemContainer.pEffectiveRhs.swap(p_eff_rhs);

        auto p_eff_dx = Kratos::make_shared<TSystemVectorType>(n_free_dofs);
        rLinearSystemContainer.pEffectiveDx.swap(p_eff_dx);
    }

    void AllocateLinearSystemConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        const DofArrayUtilities::SlaveToMasterDofsMap& rSlaveToMasterDofsMap,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Check if there are master-slave constraints
        if (!rSlaveToMasterDofsMap.empty()) {
            // Fill the master-slave constraints graph
            TSparseGraphType constraints_sparse_graph;
            this->SetUpMasterSlaveConstraintsGraph(rDofSet, rEffectiveDofSet, rSlaveToMasterDofsMap, constraints_sparse_graph);

            // Allocate the constraints arrays (note that we are using the move assignment operator in here)
            auto p_aux_q = Kratos::make_shared<TSystemVectorType>(rDofSet.size());
            rLinearSystemContainer.pConstraintsQ.swap(p_aux_q);

            auto p_aux_T = Kratos::make_shared<TSparseMatrixType>(constraints_sparse_graph);
            rLinearSystemContainer.pConstraintsT.swap(p_aux_T);

            // // Free the constraints sparse matrix graph memory to avoid having two graphs at the same time
            // delete constraints_sparse_graph;
        }

        // Set up Dirichlet matrix sparse graph
        // Note that the row size is the effective DOF set size as the master-slave constraints act over the already effective DOF set
        KRATOS_ERROR_IF(rEffectiveDofSet.empty()) << "Effective DOF set is empty." << std::endl;
        TSparseGraphType dirichlet_sparse_graph(rEffectiveDofSet.size());

        // Loop the effective DOFs to add the free ones to the graph
        // Note that this graph results in a diagonal matrix with zeros in the fixed DOFs
        unsigned int aux_count = 0;
        for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
            // Get current DOF
            auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);

            // Check if current DOF is free and add it to the sparse graph
            if (p_dof->IsFree()) {
                dirichlet_sparse_graph.AddEntry(p_dof->EffectiveEquationId(), aux_count);
                aux_count++;
            }
        }

        // Allocate the Dirichlet constraints relation matrix
        // Note that there is no need to allocate the Dirichlet constraints relation vector as this is never used
        // as we always solve for the solution increment (the Dirichlet values are already in the effective DOF set data)
        auto p_aux_T = Kratos::make_shared<TSparseMatrixType>(dirichlet_sparse_graph);
        rLinearSystemContainer.pDirichletT.swap(p_aux_T);
    }

    void BuildDirichletConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Set ones in the entries of the Dirichlet constraints relation matrix
        rLinearSystemContainer.pDirichletT->SetValue(1.0);
    }

    //FIXME: Do the RHS-only version
    void ApplyLinearSystemConstraints(
        const DofsArrayType& rEffectiveDofArray,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Get effective arrays
        auto p_eff_dx = rLinearSystemContainer.pEffectiveDx;
        auto p_eff_rhs = rLinearSystemContainer.pEffectiveRhs;

        // Initialize the effective RHS
        KRATOS_ERROR_IF(p_eff_rhs == nullptr) << "Effective RHS vector has not been initialized yet." << std::endl;
        p_eff_rhs->SetValue(0.0);

        // Initialize the effective solution vector
        KRATOS_ERROR_IF(p_eff_dx == nullptr) << "Effective solution increment vector has not been initialized yet." << std::endl;
        p_eff_dx->SetValue(0.0);

        // Check if there are master-slave constraints to do the constraints composition
        const auto& r_model_part = this->GetModelPart();
        const std::size_t n_constraints = r_model_part.NumberOfMasterSlaveConstraints();
        if (n_constraints) { //FIXME: In here we should check the number of active constraints
            // Compute the total relation matrix including master-slave and Dirichlet constraints
            auto& r_dirichlet_T = *rLinearSystemContainer.pDirichletT;
            auto& r_constraints_T = *rLinearSystemContainer.pConstraintsT;
            rLinearSystemContainer.pEffectiveT = AmgclCSRSpMMUtilities::SparseMultiply(r_constraints_T, r_dirichlet_T);

            // Apply constraints to RHS
            auto p_rhs = rLinearSystemContainer.pRhs;
            rLinearSystemContainer.pEffectiveT->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            auto p_lhs = rLinearSystemContainer.pLhs;
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *rLinearSystemContainer.pEffectiveT);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rLinearSystemContainer.pEffectiveT);
            rLinearSystemContainer.pEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
        } else {
            // Assign the Dirichlet relation matrix as the effective ones since there are no other constraints
            rLinearSystemContainer.pEffectiveT = rLinearSystemContainer.pDirichletT;

            // Apply Dirichlet constraints to RHS
            auto p_rhs = rLinearSystemContainer.pRhs;
            rLinearSystemContainer.pEffectiveT->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply Dirichlet constraints to LHS
            auto p_lhs = rLinearSystemContainer.pLhs;
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *rLinearSystemContainer.pEffectiveT);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rLinearSystemContainer.pEffectiveT);
            rLinearSystemContainer.pEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
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
