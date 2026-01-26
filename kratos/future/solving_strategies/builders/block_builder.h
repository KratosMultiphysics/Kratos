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
template<class TLinearAlgebra>
class BlockBuilder : public Builder<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BlockBuilder
    KRATOS_CLASS_POINTER_DEFINITION(BlockBuilder);

    /// Base builder type definition
    using BaseType = Builder<TLinearAlgebra>;

    /// Index type definition from sparse matrix
    using IndexType = typename TLinearAlgebra::IndexType;

    /// Matrix type definition
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Vector type definition
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Sparse graph type definition
    using SparseGraphType = typename TLinearAlgebra::SparseGraphType;

    /// DOF type definition
    using DofType = typename BaseType::DofType;

    /// DOF array type definition
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// DOF pointer vector type definition
    using DofPointerVectorType = typename BaseType::DofPointerVectorType;

    /// Linear system type definition
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

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
        mScalingType = Settings["scaling_type"].GetString();
    }

    virtual ~BlockBuilder() = default;

    ///@}
    ///@name Operations
    ///@{

    void AllocateLinearSystem(
        const SparseGraphType& rSparseGraph,
        ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData) override
    {
        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_dx = Kratos::make_shared<VectorType>(rSparseGraph);
        auto p_rhs = Kratos::make_shared<VectorType>(rSparseGraph);
        auto p_lhs = Kratos::make_shared<MatrixType>(rSparseGraph);

        // Set the linear system with the arrays above
        auto p_lin_sys = Kratos::make_shared<LinearSystemType>(p_lhs, p_rhs, p_dx, "LinearSystem");
        rImplicitStrategyData.pSetLinearSystem(p_lin_sys);

        // Set the effective arrays and corresponding linear system
        if (rImplicitStrategyData.RequiresEffectiveDofSet()) {
            // If there are no constraints the effective arrays are the same as the global ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            auto p_eff_lin_sys = Kratos::make_shared<LinearSystemType>(p_lhs, p_rhs, p_dx, "EffectiveLinearSystem");
            rImplicitStrategyData.pSetEffectiveLinearSystem(p_eff_lin_sys);
        } else {
            // Allocate the effective vectors according to the effective DOF set size
            auto p_eff_lhs = Kratos::make_shared<MatrixType>();
            auto p_eff_rhs = Kratos::make_shared<VectorType>(rImplicitStrategyData.pGetEffectiveDofSet()->size());
            auto p_eff_dx = Kratos::make_shared<VectorType>(rImplicitStrategyData.pGetEffectiveDofSet()->size());
            auto p_eff_lin_sys = Kratos::make_shared<LinearSystemType>(p_eff_lhs, p_eff_rhs, p_eff_dx, "EffectiveLinearSystem");
            rImplicitStrategyData.pSetEffectiveLinearSystem(p_eff_lin_sys);
        }
    }

    void AllocateLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        // Check if there are master-slave constraints
        //TODO: Do the MPI sum all
        const std::size_t n_constraints = this->GetModelPart().NumberOfMasterSlaveConstraints();
        if (n_constraints) {
            // Fill the master-slave constraints graph
            SparseGraphType constraints_sparse_graph;
            auto& r_dof_set = *rImplicitStrategyData.pGetDofSet();
            auto& r_eff_dof_set = *rImplicitStrategyData.pGetEffectiveDofSet();
            this->SetUpMasterSlaveConstraintsGraph(r_dof_set, r_eff_dof_set, constraints_sparse_graph);

            // Allocate the constraints arrays (note that we are using the move assignment operator in here)
            auto p_aux_q = Kratos::make_shared<VectorType>(r_dof_set.size());
            rImplicitStrategyData.pSetConstraintsQ(p_aux_q);

            auto p_aux_T = Kratos::make_shared<MatrixType>(constraints_sparse_graph);
            rImplicitStrategyData.pSetConstraintsT(p_aux_T);
        }
    }

    //FIXME: Do the RHS-only version
    void ApplyLinearSystemConstraints(ImplicitStrategyData<TLinearAlgebra>& rImplicitStrategyData) override
    {
        // Calculate the effective LHS, RHS and solution vector
        ApplyBlockBuildMasterSlaveConstraints(rImplicitStrategyData);

        // Apply the Dirichlet BCs in a block way by leveraging the CSR matrix implementation
        auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        auto& r_eff_lhs = *(p_eff_lin_sys->pGetLeftHandSide());
        auto& r_eff_rhs = *(p_eff_lin_sys->pGetRightHandSide());
        auto& r_eff_dof_set = *(rImplicitStrategyData.pGetEffectiveDofSet());
        ApplyBlockBuildDirichletConditions(r_eff_dof_set, r_eff_lhs, r_eff_rhs);
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    std::string mScalingType; // Diagonal scaling type

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Applies the master-slave constraints
     * This method applies the master-slave constraints following a block-type build
     * @param rImplicitStrategyData Auxiliary container with the linear system arrays
     */
    void ApplyBlockBuildMasterSlaveConstraints(ImplicitStrategyData<TLinearAlgebra> &rImplicitStrategyData)
    {
        const auto& r_model_part = this->GetModelPart();
        auto p_lin_sys = rImplicitStrategyData.pGetLinearSystem();
        auto p_eff_lin_sys = rImplicitStrategyData.pGetEffectiveLinearSystem();
        const std::size_t n_constraints = r_model_part.NumberOfMasterSlaveConstraints();
        if (n_constraints) { //FIXME: In here we should check the number of active constraints
            // Get effective arrays
            auto p_eff_dx = p_eff_lin_sys->pGetSolution();
            auto p_eff_rhs = p_eff_lin_sys->pGetRightHandSide();

            // Initialize the effective RHS
            KRATOS_ERROR_IF(p_eff_rhs == nullptr) << "Effective RHS vector has not been initialized yet." << std::endl;
            p_eff_rhs->SetValue(0.0);

            // Initialize the effective solution vector
            KRATOS_ERROR_IF(p_eff_dx == nullptr) << "Effective solution increment vector has not been initialized yet." << std::endl;
            p_eff_dx->SetValue(0.0);

            // Assign the master-slave relation matrix as the effective one since there are no other constraints
            auto p_const_T = rImplicitStrategyData.pGetConstraintsT();
            rImplicitStrategyData.pSetEffectiveT(p_const_T);

            // Apply constraints to RHS
            auto p_rhs = p_lin_sys->pGetRightHandSide();
            auto p_eff_T = rImplicitStrategyData.pGetEffectiveT();
            p_eff_T->TransposeSpMV(*p_rhs, *p_eff_rhs);

            // Apply constraints to LHS
            //TODO: Rethink once we figure out the LinearSystemContainer, LinearOperator and so...
            auto p_lhs = p_lin_sys->pGetLeftHandSide();
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*p_lhs, *rImplicitStrategyData.pGetEffectiveT());
            auto p_transT = AmgclCSRConversionUtilities::Transpose(*rImplicitStrategyData.pGetEffectiveT());
            auto p_eff_lhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);
            p_eff_lin_sys->SetLeftHandSide(*p_eff_lhs);

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
            auto& r_dx = p_lin_sys->GetSolution();
            auto& r_lhs = p_lin_sys->GetLeftHandSide();
            auto& r_rhs = p_lin_sys->GetRightHandSide();
            p_eff_lin_sys->SetSolution(r_dx);
            p_eff_lin_sys->SetLeftHandSide(r_lhs);
            p_eff_lin_sys->SetRightHandSide(r_rhs);
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
        MatrixType& rLHS,
        VectorType& rRHS) const
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

        //TODO: Implement this in the CSR matrix or here? @Pooyan and @Ruben think that should go to a csr_utilities.h
        // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
        // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        // Get the diagonal scaling factor
        const double diagonal_value = this->GetDiagonalScalingFactor(rLHS);

        // Apply the free DOFs (i.e., fixity) vector to the system arrays
        //TODO: Why don't we put the ApplyHomogeneousDirichlet in a csr_utilities.h to eventually have a different implementation for custom needs
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
        VectorType& rRHS) const
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
    virtual double GetDiagonalScalingFactor(const MatrixType& rLHS) const
    {
        if (mScalingType == "no_scaling") {
            return 1.0;
        } else if (mScalingType == "norm_diagonal") {
            return rLHS.NormDiagonal();
        } else if (mScalingType == "norm_diagonal_ratio") {
            return rLHS.NormDiagonal() / rLHS.size1();
        } else if (mScalingType == "max_diagonal") {
            return rLHS.MaxDiagonal();
        } else if (mScalingType == "max_diagonal_ratio") {
            return rLHS.MaxDiagonal() / rLHS.size1();
        } else if (mScalingType == "prescribed_diagonal") {
            const auto& r_process_info = (this->GetModelPart()).GetProcessInfo();
            KRATOS_ERROR_IF_NOT(r_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined in ProcessInfo container. Please set 'BUILD_SCALE_FACTOR' variable." << std::endl;
            return r_process_info.GetValue(BUILD_SCALE_FACTOR);
        } else {
            KRATOS_ERROR << "Wrong scaling type: \'" << mScalingType << "\'. Available options are:\n" <<
                "\t- \'no_scaling\'\n" <<
                "\t- \'norm_diagonal\'\n" <<
                "\t- \'norm_diagonal_ratio\'\n" <<
                "\t- \'max_diagonal\'\n" <<
                "\t- \'max_diagonal_ratio\'\n" <<
                "\t- \'prescribed_diagonal\'\n" << std::endl;
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
