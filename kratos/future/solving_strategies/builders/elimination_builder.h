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
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

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

    void AllocateLinearSystem(
        const SparseGraphType& rSparseGraph,
        ImplicitStrategyDataContainer<TLinearAlgebra> &rImplicitStrategyDataContainer) override
    {
        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_dx = Kratos::make_shared<VectorType>(rSparseGraph);
        auto p_rhs = Kratos::make_shared<VectorType>(rSparseGraph);
        auto p_lhs = Kratos::make_shared<MatrixType>(rSparseGraph);

        // Set the linear system with the arrays above
        auto p_lin_sys = Kratos::make_unique<LinearSystemType>(p_lhs, p_rhs, p_dx, "LinearSystem");
        rImplicitStrategyDataContainer.mpLinearSystem = std::move(p_lin_sys); // Transfer ownership to the data container

        // Get the number of free DOFs to allocate the effective system arrays
        auto& r_eff_dof_set = *(rImplicitStrategyDataContainer.pEffectiveDofSet);
        const std::size_t n_free_dofs = IndexPartition<std::size_t>(r_eff_dof_set.size()).for_each<SumReduction<std::size_t>>([&](IndexType Index) {
            auto p_dof = r_eff_dof_set.begin() + Index;
            return p_dof->IsFixed() ? 0 : 1;
        });

        // Allocate the effective arrays according to the number of free effective DOFs
        auto p_eff_lhs = Kratos::make_shared<MatrixType>();
        auto p_eff_rhs = Kratos::make_shared<VectorType>(n_free_dofs);
        auto p_eff_dx = Kratos::make_shared<VectorType>(n_free_dofs);

        // Set the effective linear system with the effective arrays
        auto p_eff_lin_sys = Kratos::make_unique<LinearSystemType>(p_eff_lhs, p_eff_rhs, p_eff_dx, "EffectiveLinearSystem");
        rImplicitStrategyDataContainer.mpEffectiveLinearSystem = std::move(p_eff_lin_sys); // Transfer ownership to the data container
    }

    void AllocateLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer) override
    {
        // Check if there are master-slave constraints
        auto& r_eff_dof_set = *(rImplicitStrategyDataContainer.pEffectiveDofSet);
        const std::size_t n_constraints = this->GetModelPart().NumberOfMasterSlaveConstraints();
        if (n_constraints) {
            // Fill the master-slave constraints graph
            SparseGraphType constraints_sparse_graph;
            auto& r_dof_set = *(rImplicitStrategyDataContainer.pDofSet);
            this->SetUpMasterSlaveConstraintsGraph(r_dof_set, r_eff_dof_set, constraints_sparse_graph);

            // Allocate the constraints arrays (note that we are using the move assignment operator in here)
            auto p_aux_q = Kratos::make_shared<VectorType>(r_dof_set.size());
            rImplicitStrategyDataContainer.pConstraintsQ.swap(p_aux_q);

            auto p_aux_T = Kratos::make_shared<MatrixType>(constraints_sparse_graph);
            rImplicitStrategyDataContainer.pConstraintsT.swap(p_aux_T);

            // // Free the constraints sparse matrix graph memory to avoid having two graphs at the same time
            // delete constraints_sparse_graph;
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

    //FIXME: Do the RHS-only version
    void ApplyLinearSystemConstraints(ImplicitStrategyDataContainer<TLinearAlgebra>& rImplicitStrategyDataContainer) override
    {
        // Get effective arrays
        auto p_eff_lin_sys = rImplicitStrategyDataContainer.pGetEffectiveLinearSystem();
        auto p_eff_dx = p_eff_lin_sys->pGetSolution();
        auto p_eff_rhs = p_eff_lin_sys->pGetRightHandSide();

        // Initialize the effective RHS
        KRATOS_ERROR_IF(p_eff_rhs == nullptr) << "Effective RHS vector has not been initialized yet." << std::endl;
        p_eff_rhs->SetValue(0.0);

        // Initialize the effective solution vector
        KRATOS_ERROR_IF(p_eff_dx == nullptr) << "Effective solution increment vector has not been initialized yet." << std::endl;
        p_eff_dx->SetValue(0.0);

        // Set ones in the entries of the Dirichlet constraints relation matrix
        mpDirichletT->SetValue(1.0);

        // Get the linear system to apply the constraints to
        auto p_lin_sys = rImplicitStrategyDataContainer.pGetLinearSystem();
        auto& r_lhs = p_lin_sys->GetLeftHandSide();
        auto& r_rhs = p_lin_sys->GetRightHandSide();

        // Check if there are master-slave constraints to do the constraints composition
        const auto& r_model_part = this->GetModelPart();
        const std::size_t n_constraints = r_model_part.NumberOfMasterSlaveConstraints();
        if (n_constraints) { //FIXME: In here we should check the number of active constraints
            // Compute the total relation matrix including master-slave and Dirichlet constraints
            auto& r_constraints_T = *rImplicitStrategyDataContainer.pConstraintsT;
            rImplicitStrategyDataContainer.pEffectiveT = AmgclCSRSpMMUtilities::SparseMultiply(r_constraints_T, *mpDirichletT);

            // Apply constraints to RHS
            rImplicitStrategyDataContainer.pEffectiveT->TransposeSpMV(r_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(r_lhs, *rImplicitStrategyDataContainer.pEffectiveT);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rImplicitStrategyDataContainer.pEffectiveT);
            auto p_eff_lhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
            p_eff_lin_sys->SetLeftHandSide(*p_eff_lhs);
        } else {
            // Assign the Dirichlet relation matrix as the effective ones since there are no other constraints
            rImplicitStrategyDataContainer.pEffectiveT = mpDirichletT;

            // Apply Dirichlet constraints to RHS
            rImplicitStrategyDataContainer.pEffectiveT->TransposeSpMV(r_rhs, *p_eff_rhs);

            // Apply Dirichlet constraints to LHS
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(r_lhs, *rImplicitStrategyDataContainer.pEffectiveT);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rImplicitStrategyDataContainer.pEffectiveT);
            auto p_eff_lhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
            p_eff_lin_sys->SetLeftHandSide(*p_eff_lhs);
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
