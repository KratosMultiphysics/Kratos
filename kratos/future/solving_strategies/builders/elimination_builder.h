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
        //
    }

    //FIXME: Do the RHS-only version
    void ApplyLinearSystemConstraints(
        const DofsArrayType& rEffectiveDofArray,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
        //FIXME: To be implemented --> see the testbench notebook
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
