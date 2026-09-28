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
#include "future/containers/linear_system.h"
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
template<class TLinearAlgebra>
class EliminationBuilder : public Builder<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of EliminationBuilder
    KRATOS_CLASS_POINTER_DEFINITION(EliminationBuilder);

    /// Base builder type definition
    using BaseType = Builder<TLinearAlgebra>;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Sparse graph type definition
    using SparseGraphType = typename TLinearAlgebra::SparseGraphType;

    /// Index type definition from sparse matrix
    using IndexType = typename TLinearAlgebra::IndexType;

    /// DOF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Linear system type definition
    using LinearSystemType = typename BaseType::LinearSystemType;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

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

    void SetDofEquationIds(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        // Get the DOFs and effective DOFs containers and check they are not empty
        auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        KRATOS_ERROR_IF(r_dof_set.empty()) << "DOFs set is empty. Set up the DOFs array first." << std::endl;

        // Initialize the free and fixed effective DOFs counters
        // Note that the free DOFs ids start by 0 to position them at beginning of the system
        // while the fixed ones are positioned at the end (in opposite order) by starting the ids from the number of effective DOFs
        std::size_t free_id = 0;
        std::size_t fix_id = r_dof_set.size();

        // Set the DOFs' effective equation global ids
        // The free degrees of freedom are positioned at the beginning of the system, while the fixed ones are at the end (in opposite order)
        // That means that if the EquationId is greater than the fix_id (i.e., the number of free DOFs) then it means that the pointed degree of freedom is restrained
        for (auto it_dof = r_dof_set.begin(); it_dof != r_dof_set.end(); ++it_dof) {
            if (it_dof->IsFixed()) {
                it_dof->SetEquationId(--fix_id);
            } else {
                it_dof->SetEquationId(free_id++);
            }
        }
    }

    void SetDofEffectiveEquationIds(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        KRATOS_TRY

        // Get the DOFs and effective DOFs containers and check they are not empty
        auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
        KRATOS_ERROR_IF(r_dof_set.empty()) << "DOFs set is empty. Set up the DOFs array first." << std::endl;
        KRATOS_ERROR_IF(r_eff_dof_set.empty()) << "Effective DOFs set is empty. Set up the effective DOFs array first." << std::endl;

        // Check if the effective and "standard" containers are the same
        // We do it with the addresses to avoid checking the content (i.e., each DOF one-by-one)
        if (&r_eff_dof_set == &r_dof_set) {
            // Set the DOFs' effective equation global ids to match the standard ones
            // Note that these already account for the free and fixed DOFs forward and backward positioning in the system
            const std::size_t n_free_dofs = (IndexPartition<IndexType>(r_eff_dof_set.size())).template for_each<SumReduction<std::size_t>>([&](IndexType Index) {
                auto it_dof = r_eff_dof_set.begin() + Index;
                it_dof->SetEffectiveEquationId(it_dof->EquationId());
                return it_dof->IsFixed() ? 0 : 1;
            });

            // Set the problem size to the number of free effective DOFs
            this->SetProblemSize(n_free_dofs);
        } else {
            // Initialize the free and fixed effective DOFs counters
            // Note that the free DOFs ids start by 0 to position them at beginning of the system
            // while the fixed ones are positioned at the end (in opposite order) by starting the ids from the number of effective DOFs
            std::size_t free_id = 0;
            std::size_t fix_id = r_eff_dof_set.size();

            // Initialize all DOFs effective equation ids to the maximum allowable value
            // Note that this makes possible to distingish the effective DOFs from the non-effective ones
            IndexPartition<IndexType>(r_dof_set.size()).for_each([&](IndexType Index) {
                auto it_dof = r_dof_set.begin() + Index;
                it_dof->SetEffectiveEquationId(std::numeric_limits<typename Node::DofType::EquationIdType>::max());
            });

            // Set the effective DOFs' effective equation global ids
            // The free degrees of freedom are positioned at the beginning of the system, while the fixed ones are at the end (in opposite order)
            // That means that if the EquationId is greater than the fix_id (i.e., the number of free DOFs) then it means that the pointed degree of freedom is restrained
            for (auto it_eff_dof = r_eff_dof_set.begin(); it_eff_dof != r_eff_dof_set.end(); ++it_eff_dof) {
                if (it_eff_dof->IsFixed()) {
                    it_eff_dof->SetEffectiveEquationId(--fix_id);
                } else {
                    it_eff_dof->SetEffectiveEquationId(free_id++);
                }
            }

            // Set the problem size to current fix id (i.e., number of free effective DOFs)
            this->SetProblemSize(fix_id);
        }

        KRATOS_CATCH("");
    }

    void AllocateEffectiveLinearSystem(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        // Allocate the effective arrays according to the number of free effective DOFs
        const std::size_t n_free_dofs = this->GetProblemSize(); // Problem size is already set to the number of free effective DOFs
        auto p_eff_lhs = Kratos::make_shared<MatrixType>();
        auto p_eff_rhs = Kratos::make_shared<VectorType>(n_free_dofs);
        auto p_eff_dx = Kratos::make_shared<VectorType>(n_free_dofs);

        // Set the effective linear system with the effective arrays
        auto p_eff_lin_sys = Kratos::make_shared<LinearSystemType>(p_eff_lhs, p_eff_rhs, p_eff_dx, "EffectiveLinearSystem");
        rImplicitStrategyData.pSetEffectiveLinearSystem(p_eff_lin_sys);
    }

    void AllocateLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        // Check if there are master-slave constraints
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
        const std::size_t n_constraints = this->GetModelPart().NumberOfMasterSlaveConstraints();
        if (n_constraints) {
            // Fill the master-slave constraints graph
            SparseGraphType constraints_sparse_graph;
            auto& r_dof_set = *(rImplicitStrategyData.pGetDofSet());
            this->SetUpMasterSlaveConstraintsGraph(r_dof_set, r_eff_dof_set, constraints_sparse_graph);

            // Allocate the constraints arrays
            auto p_aux_q = Kratos::make_shared<VectorType>(r_dof_set.size());
            rImplicitStrategyData.pSetConstraintsQ(p_aux_q);

            auto p_aux_T = Kratos::make_shared<MatrixType>(constraints_sparse_graph);
            rImplicitStrategyData.pSetConstraintsT(p_aux_T);
        }

        // Set up Dirichlet matrix sparse graph
        // Note that the row size is the effective DOF set size as the master-slave constraints act over the already effective DOF set
        KRATOS_ERROR_IF(r_eff_dof_set.empty()) << "Effective DOF set is empty." << std::endl;
        SparseGraphType dirichlet_sparse_graph(r_eff_dof_set.size());

        // Loop the effective DOFs to add the free ones to the graph
        // Note that this graph results in a diagonal matrix with zeros in the fixed DOFs
        unsigned int aux_count = 0;
        for (IndexType i_dof = 0; i_dof < r_eff_dof_set.size(); ++i_dof) {
            // Get current DOF
            auto p_dof = *(r_eff_dof_set.ptr_begin() + i_dof);

            // Check if current DOF is free and add it to the sparse graph
            if (p_dof->IsFree()) {
                dirichlet_sparse_graph.AddEntry(p_dof->EffectiveEquationId(), aux_count);
                aux_count++;
            }
        }

        // Allocate the Dirichlet constraints relation matrix
        // Note that there is no need to allocate the Dirichlet constraints relation vector as this is never used
        // as we always solve for the solution increment (the Dirichlet values are already in the effective DOF set data)
        auto p_aux_T = Kratos::make_shared<MatrixType>(dirichlet_sparse_graph);
        mpDirichletT.swap(p_aux_T);
    }

    void ApplyLinearSystemConstraints(
        ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData,
        const bool SkipLeftHandSide = false) override
    {
        // Get effective arrays
        auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        auto& r_eff_dx = *(p_eff_lin_sys->pGetVector(DenseVectorTag::Dx));
        auto& r_eff_rhs = *(p_eff_lin_sys->pGetVector(DenseVectorTag::RHS));

        // Initialize the effective RHS
        r_eff_rhs.SetValue(0.0);

        // Initialize the effective solution vector
        r_eff_dx.SetValue(0.0);

        // Set ones in the entries of the Dirichlet constraints relation matrix
        mpDirichletT->SetValue(1.0);

        // Get the linear system to apply the constraints to
        auto p_lin_sys = rImplicitStrategyData.pGetLinearSystem();
        auto& r_lhs = *(p_lin_sys->pGetMatrix(SparseMatrixTag::LHS));
        auto& r_rhs = *(p_lin_sys->pGetVector(DenseVectorTag::RHS));

        // Check if there are master-slave constraints to do the constraints composition
        const auto& r_model_part = this->GetModelPart();
        const std::size_t n_constraints = r_model_part.NumberOfMasterSlaveConstraints();
        if (n_constraints) { //FIXME: In here we should check the number of active constraints
            // Compute the total relation matrix including master-slave and Dirichlet constraints
            auto& r_constraints_T = *rImplicitStrategyData.pGetConstraintsT();
            rImplicitStrategyData.pSetEffectiveT(AmgclCSRSpMMUtilities::SparseMultiply(r_constraints_T, *mpDirichletT));

            // Apply constraints to RHS
            rImplicitStrategyData.pGetEffectiveT()->TransposeSpMV(r_rhs, r_eff_rhs);

            // Apply constraints to LHS
            if (!SkipLeftHandSide) {
                auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(r_lhs, *rImplicitStrategyData.pGetEffectiveT());
                auto p_transT = AmgclCSRConversionUtilities::Transpose(*rImplicitStrategyData.pGetEffectiveT());
                auto p_eff_lhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
                p_eff_lin_sys->pSetMatrix(p_eff_lhs, SparseMatrixTag::LHS);
            }
        } else {
            // Assign the Dirichlet relation matrix as the effective ones since there are no other constraints
            rImplicitStrategyData.pSetEffectiveT(mpDirichletT);

            // Apply Dirichlet constraints to RHS
            rImplicitStrategyData.pGetEffectiveT()->TransposeSpMV(r_rhs, r_eff_rhs);

            // Apply Dirichlet constraints to LHS
            if (!SkipLeftHandSide) {
                auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(r_lhs, *rImplicitStrategyData.pGetEffectiveT());
                auto p_transT = AmgclCSRConversionUtilities::Transpose(*rImplicitStrategyData.pGetEffectiveT());
                auto p_eff_lhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
                p_eff_lin_sys->pSetMatrix(p_eff_lhs, SparseMatrixTag::LHS);
            }
        }
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    typename MatrixType::Pointer mpDirichletT = nullptr; // Dirichlet constraints relation matrix

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
