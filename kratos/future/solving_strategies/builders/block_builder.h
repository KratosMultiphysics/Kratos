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
 * @brief Utility class for handling the block build
 * @details This helper class extends the base Builder class to consider an block type build
 * The block type build keeps the Dirichlet DOFs in the the linear system of equations to be solved.
 * Hence, the Dirichlet BCs imposition is achieved by removing setting the corresponding residual
 * rows to zero and setting a diagonal contribution in the corresponding row of the Left Hand Side
 * matrix. The value of the diagonal term is computed according to the scaling type.
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class BlockBuilder : public Builder<TThreadLocalStorage, TSparseMatrixType, TSystemVectorType, TSparseGraphType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Scaling type definition
    enum class ScalingType
    {
        NoScaling,
        NormDiagonal,
        MaxDiagonal,
        PrescribedDiagonal
    };

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

        // Set scaling type
        const std::string scaling_type = Settings["scaling_type"].GetString();
            if (scaling_type == "no_scaling") {
                mScalingType = ScalingType::NoScaling;
            } else if (scaling_type == "norm_diagonal") {
                mScalingType = ScalingType::NormDiagonal;
            } else if (scaling_type == "max_diagonal") {
                mScalingType = ScalingType::MaxDiagonal;
            } else if (scaling_type == "prescribed_diagonal") {
                mScalingType = ScalingType::PrescribedDiagonal;
            } else {
                KRATOS_ERROR << "Provided \'" << scaling_type << "\' is not supported." << std::endl;
            }
    }

    virtual ~BlockBuilder() = default;

    ///@}
    ///@name Operations
    ///@{

    void AllocateLinearSystemArrays(
        const typename DofsArrayType::Pointer pDofSet,
        const typename DofsArrayType::Pointer pEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer) override
    {
        // Set up the system sparse matrix graph (note that the sparse graph will be destroyed when leaving this scope)
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

        // Set the effective arrays
        if (pDofSet == pEffectiveDofSet) {
            // If there are no constraints the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rLinearSystemContainer.pEffectiveLhs = rLinearSystemContainer.pLhs;
            rLinearSystemContainer.pEffectiveRhs = rLinearSystemContainer.pRhs;
            rLinearSystemContainer.pEffectiveDx = rLinearSystemContainer.pDx;
        } else {
            // Allocate the effective vectors according to the effective DOF set size
            auto p_eff_lhs = Kratos::make_shared<TSparseMatrixType>();
            rLinearSystemContainer.pEffectiveLhs.swap(p_eff_lhs);

            auto p_eff_rhs = Kratos::make_shared<TSystemVectorType>(pEffectiveDofSet->size());
            rLinearSystemContainer.pEffectiveRhs.swap(p_eff_rhs);

            auto p_eff_dx = Kratos::make_shared<TSystemVectorType>(pEffectiveDofSet->size());
            rLinearSystemContainer.pEffectiveDx.swap(p_eff_dx);
        }
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
        }
    }

    void BuildDirichletConstraints(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer) override
    {
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

    ScalingType mScalingType; // Diagonal scaling type

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies the master-slave constraints
     * This method applies the master-slave constraints following a block-type build
     * @param rLinearSystemContainer Auxiliary container with the linear system arrays
     */
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
            KRATOS_ERROR_IF(p_eff_dx == nullptr) << "Effective solution increment vector has not been initialized yet." << std::endl;
            p_eff_dx->SetValue(0.0);

            // Assign the master-slave relation matrix as the effective one since there are no other constraints
            rLinearSystemContainer.pEffectiveT = rLinearSystemContainer.pConstraintsT;

            // Apply constraints to RHS
            auto p_rhs = rLinearSystemContainer.pRhs;
            rLinearSystemContainer.pEffectiveT->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            auto p_lhs = rLinearSystemContainer.pLhs;
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *rLinearSystemContainer.pEffectiveT);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rLinearSystemContainer.pEffectiveT);
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

    /**
     * @brief Applies the Dirichlet conditions to the system
     * This method applies the Dirichlet conditions following a block-type build
     * @param rDofArray DOFs array from which the fixity is checked
     * @param rLHS Left Hand Side matrix
     * @param rRHS Right Hand Side vector
     */
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

    /**
     * @brief Applies the Dirichlet conditions to the system
     * This method applies the Dirichlet conditions following a block-type build
     * @param rDofArray DOFs array from which the fixity is checked
     * @param rRHS Right Hand Side vector
     */
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

    /**
     * @brief Get the Diagonal Scaling Factor
     * This method calculates the diagonal scaling value to be used in the Dirichlet BCs imposition from the provided LHS
     * @param rLHS Left Hand Side matrix
     * @return double Diagonal scaling factor
     */
    virtual double GetDiagonalScalingFactor(const TSparseMatrixType& rLHS) const
    {
        if (mScalingType == ScalingType::NoScaling) {
            return 1.0;
        } else if (mScalingType == ScalingType::NormDiagonal) {
            return rLHS.NormDiagonal();
            // return rLHS.NormDiagonal() / rLHS.size1(); //TODO: Decide which one
        } else if (mScalingType == ScalingType::MaxDiagonal) {
            return rLHS.MaxDiagonal();
            // return rLHS.MaxDiagonal() / rLHS.size1(); // TODO: Decide which one
        } else if (mScalingType == ScalingType::PrescribedDiagonal) {
            const auto& r_process_info = (this->GetModelPart()).GetProcessInfo();
            KRATOS_ERROR_IF_NOT(r_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined in ProcessInfo container. Please set 'BUILD_SCALE_FACTOR' variable." << std::endl;
            return r_process_info.GetValue(BUILD_SCALE_FACTOR);
        } else {
            KRATOS_ERROR << "Wrong scaling type." << std::endl;
        }
    }

    ///@}
}; // Class BlockBuilder

///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
