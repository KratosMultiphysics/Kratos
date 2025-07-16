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
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
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

    void ResizeAndInitializeVectors(
        const TSparseGraphType& rSparseGraph,
        const typename DofsArrayType::Pointer pDofSet,
        const typename DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer) override
    {
        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_dx = Kratos::make_shared<TSystemVectorType>(rSparseGraph);
        rLinearSystemContainer.pDx.swap(p_dx);

        auto p_rhs = Kratos::make_shared<TSystemVectorType>(rSparseGraph);
        rLinearSystemContainer.pRhs.swap(p_rhs);

        auto p_lhs = Kratos::make_shared<TSparseMatrixType>(rSparseGraph);
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

    void ConstructDirichletConstraintsStructure(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Set up Dirichlet matrix sparse graph (note that the row size is the effective DOF set size from the block build)
        KRATOS_ERROR_IF(rEffectiveDofSet.empty()) << "Effective DOF set is empty." << std::endl;
        TSparseGraphType constraints_sparse_graph(rEffectiveDofSet.size());

        // Loop the effective DOFs to add the free ones to the graph
        // Note that this graph results in a diagonal matrix with zeros in the fixed DOFs
        unsigned int aux_count = 0;
        for (IndexType i_dof = 0; i_dof < rEffectiveDofSet.size(); ++i_dof) {
            // Get current DOF
            auto p_dof = *(rEffectiveDofSet.ptr_begin() + i_dof);

            // Check if current DOF is free and add it to the sparse graph
            if (p_dof->IsFree()) {
                constraints_sparse_graph.AddEntry(p_dof->EffectiveEquationId(), aux_count);
                aux_count++;
            }
        }

        // Allocate the constraints arrays (note that we are using the move assignment operator in here)
        auto p_aux_q = Kratos::make_shared<TSystemVectorType>(rEffectiveDofSet.size());
        rLinearSystemContainer.pDirichletQ.swap(p_aux_q);

        auto p_aux_T = Kratos::make_shared<TSparseMatrixType>(constraints_sparse_graph);
        rLinearSystemContainer.pDirichletT.swap(p_aux_T);
    }

    void BuildDirichletConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        // Set the BCs values in the Dirichlet constraints constant vector
        const auto dof_begin = rEffectiveDofSet.begin();
        auto& r_dirichlet_q = *rLinearSystemContainer.pDirichletQ;
        IndexPartition<std::size_t>(rEffectiveDofSet.size()).for_each([&](IndexType Index){
            auto p_dof = dof_begin + Index;
            if (p_dof->IsFixed()) {
                r_dirichlet_q[p_dof->EffectiveEquationId()] = p_dof->GetSolutionStepValue();
            } else {
                r_dirichlet_q[p_dof->EffectiveEquationId()] = 0.0;
            }
        });

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
            // Compute the total relation matrix
            auto p_effective_T = rLinearSystemContainer.pEffectiveT;
            auto& r_dirichlet_T = *rLinearSystemContainer.pDirichletT;
            auto& r_constraints_T = *rLinearSystemContainer.pConstraintsT;
            p_effective_T = AmgclCSRSpMMUtilities::SparseMultiply(r_constraints_T, r_dirichlet_T);

            // Compute the total constant vector
            auto& r_effective_q = *rLinearSystemContainer.pEffectiveQ;
            auto& r_dirichlet_q = *rLinearSystemContainer.pDirichletQ;
            auto& r_constraints_q = *rLinearSystemContainer.pConstraintsQ;
            r_effective_q = r_constraints_q;
            r_constraints_T.SpMV(r_dirichlet_q, r_effective_q);

            // Apply constraints to RHS
            auto p_rhs = rLinearSystemContainer.pRhs;
            p_effective_T->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            auto p_lhs = rLinearSystemContainer.pLhs;
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *p_effective_T);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*p_effective_T);
            rLinearSystemContainer.pEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
        } else {
            // Assign the Dirichlet constraints arrays as the effective ones since there are no other constraints
            rLinearSystemContainer.pEffectiveT = rLinearSystemContainer.pDirichletT;
            rLinearSystemContainer.pEffectiveQ = rLinearSystemContainer.pDirichletQ;

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
