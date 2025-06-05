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

//FIXME: We should place this somewhere else
/**
 * @brief Auxiliary container to store the linear system
 * This auxiliary container is intended to store all the arrays requires for the linear system setup
 * @tparam TSparseMatrixType The sparse matrix type
 * @tparam TSystemVectorType The system vector type
 */
template <class TSparseMatrixType, class TSystemVectorType>
struct LinearSystemContainer
{
    typename TSparseMatrixType::Pointer pLhs = nullptr; // Pointer to the LHS matrix

    typename TSystemVectorType::Pointer pRhs = nullptr; // Pointer to the RHS vector

    typename TSystemVectorType::Pointer pDx = nullptr; // Pointer to the solution increment vector

    typename TSparseMatrixType::Pointer pEffectiveLhs = nullptr; // Pointer to the effective LHS matrix (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveRhs = nullptr; // Pointer to the effective RHS vector (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveDx = nullptr; // Pointer to the effective solution increment vector (i.e., after applying system constraints)

    typename TSparseMatrixType::Pointer pEffectiveT = nullptr; // Linear system constraints total relation matrix

    typename TSystemVectorType::Pointer pEffectiveQ = nullptr; // Linear system constraints total constant vector

    typename TSparseMatrixType::Pointer pConstraintsT = nullptr; // Master-slave constraints relation matrix

    typename TSystemVectorType::Pointer pConstraintsQ = nullptr; // Master-slave constraints constant vector

    typename TSparseMatrixType::Pointer pDirichletT = nullptr; // Dirichlet constraints relation matrix

    typename TSystemVectorType::Pointer pDirichletQ = nullptr; // Dirichlet constraints constant vector

    typename TSparseMatrixType::Pointer pMassMatrix = nullptr; // Pointer to the mass matrix

    typename TSparseMatrixType::Pointer pDampingMatrix = nullptr; // Pointer to the damping matrix

    void Clear()
    {
        if (pLhs != nullptr) {
            pLhs->Clear();
        }
        if (pRhs != nullptr) {
            pRhs->Clear();
        }
        if (pDx != nullptr) {
            pDx->Clear();
        }
        if (pEffectiveLhs != nullptr) {
            pEffectiveLhs->Clear();
        }
        if (pEffectiveRhs != nullptr) {
            pEffectiveRhs->Clear();
        }
        if (pEffectiveDx != nullptr) {
            pEffectiveDx->Clear();
        }
        if (pEffectiveT != nullptr) {
            pEffectiveT->Clear();
        }
        if (pEffectiveQ != nullptr) {
            pEffectiveQ->Clear();
        }
        if (pConstraintsT != nullptr) {
            pConstraintsT->Clear();
        }
        if (pConstraintsQ != nullptr) {
            pConstraintsQ->Clear();
        }
        if (pDirichletT != nullptr) {
            pDirichletT->Clear();
        }
        if (pDirichletQ != nullptr) {
            pDirichletQ->Clear();
        }
        if (pMassMatrix != nullptr) {
            pMassMatrix->Clear();
        }
        if (pDampingMatrix != nullptr) {
            pDampingMatrix->Clear();
        }
    }
};

/**
 * @class Builder
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSystemVectorType, class TSparseGraphType>
class Builder
{
public:
    ///@name Type Definitions
    ///@{

    /// Build type definition
    enum class BuildType
    {
        Block,
        Elimination
    };

    /// Scaling type definition
    enum class ScalingType
    {
        NoScaling,
        NormDiagonal,
        MaxDiagonal,
        PrescribedDiagonal
    };

    /// Pointer definition of Builder
    KRATOS_CLASS_POINTER_DEFINITION(Builder);

    /// Data type definition from sparse matrix
    using DataType = typename TSparseMatrixType::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename TSparseMatrixType::IndexType;

    /// DOF type definition
    using DofType = Dof<DataType>;

    /// DOF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Effective DOFs map type definition
    using EffectiveDofsMapType = std::unordered_map<typename DofType::Pointer, IndexType>;

    /// DOF pointer vector type definition
    using DofPointerVectorType = typename MasterSlaveConstraint::DofPointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Builder() = delete;

    /// Constructor with model part
    Builder(
        const ModelPart& rModelPart,
        Parameters Settings = Parameters(R"({})"))
    : mpModelPart(&rModelPart)
    {
        Parameters default_parameters( R"({
            "name" : "builder",
            "scaling_type" : "max_diagonal",
            "echo_level" : 0
        })");
        Settings.ValidateAndAssignDefaults(default_parameters);

        // Set build type
        const std::string build_type = Settings["name"].GetString();
        if (build_type == "block_builder") {
            mBuildType = BuildType::Block;
        } else if (build_type == "elimination_builder") {
            mBuildType = BuildType::Elimination;
        } else {
            KRATOS_ERROR << "Provided 'build_type' is '" << build_type << "'. Available options are:\n"
            << "\t- 'block_builder'\n"
            << "\t- 'elimination_builder'" << std::endl;
        }

        // Set scaling type
        const std::string scaling_type = Settings["scaling_type"].GetString();
        if (mBuildType == BuildType::Block) {
            if (scaling_type == "no_scaling") {
                mScalingType = ScalingType::NoScaling;
            } else if (scaling_type == "norm_diagonal") {
                mScalingType = ScalingType::NormDiagonal;
            } else if (scaling_type == "max_diagonal") {
                mScalingType = ScalingType::MaxDiagonal;
            } else if (scaling_type == "prescribed_diagonal") {
                mScalingType = ScalingType::PrescribedDiagonal;
            }
        } else {
            mScalingType = ScalingType::NoScaling; // Note that scaling makes no sense in the elimination build
        }

        // Set verbosity level
        mEchoLevel = Settings["echo_level"].GetInt();
    }

    virtual ~Builder() = default;

    ///@}
    ///@name Operations
    ///@{

    //TODO: In the future, add the one with
    // - LHS alone
    // - MassMatrix
    // - DampingMatrix
    // TODO: To be discussed in the future. It would be great to have one with references. This will require using the move constructors
    virtual void ResizeAndInitializeVectors(
        const DofsArrayType &rDofSet,
        const DofsArrayType &rEffectiveDofSet,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer,
        const bool ReactionVector = false)
    {
        // Set up the sparse matrix graph (note that we do not need to keep it after the resizing)
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
        TSparseGraphType sparse_graph(mEquationSystemSize);
        SetUpSparseGraph(sparse_graph);

        //TODO: I think all these must be done by the system container

        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_lhs = Kratos::make_shared<TSparseMatrixType>(sparse_graph);
        rLinearSystemContainer.pLhs.swap(p_lhs);

        auto p_dx = Kratos::make_shared<TSystemVectorType>(sparse_graph);
        rLinearSystemContainer.pDx.swap(p_dx);

        auto p_rhs = Kratos::make_shared<TSystemVectorType>(sparse_graph);
        rLinearSystemContainer.pRhs.swap(p_rhs);

        //FIXME: I think we should separate these
        // Set the effective arrays
        // In a standard case we only need to allocate the effective solution update to avoid the first Predict() call to crash
        if (rDofSet == rEffectiveDofSet) {
            // If there are no constraints, the effective DOF set matches the standard one and the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rLinearSystemContainer.pEffectiveDx = rLinearSystemContainer.pDx;
        } else {
            // Allocate the effective vectors according to the effective DOF set size
            auto p_eff_lhs = Kratos::make_shared<TSparseMatrixType>();
            rLinearSystemContainer.pEffectiveLhs.swap(p_eff_lhs);

            auto p_eff_rhs = Kratos::make_shared<TSystemVectorType>(rEffectiveDofSet.size());
            rLinearSystemContainer.pEffectiveRhs.swap(p_eff_rhs);

            auto p_eff_dx = Kratos::make_shared<TSystemVectorType>(rEffectiveDofSet.size());
            rLinearSystemContainer.pEffectiveDx.swap(p_eff_dx);
        }
    }

    virtual void ConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType>& rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class ConstructMasterSlaveConstraintsStructure." << std::endl;
    }

    virtual void SetUpSparseGraph(TSparseGraphType& rSparseGraph)
    {
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Current equation system size is 0. Sparse graph cannot be set (call \'SetUpSystemIds\' first)." << std::endl;

        //TODO: This is always thrown as the mEquationSystemSize is already in the graph. I'd use the move (not implemented constructor)
        if (!rSparseGraph.IsEmpty()) {
            KRATOS_WARNING("Builder") << "Provided sparse graph is not empty and will be cleared." << std::endl;
            rSparseGraph.Clear();
        }

        //FIXME: This works for the block build but doesn't for the elimination!!!!!!!!
        Element::EquationIdVectorType eq_ids;
        for (auto& r_elem : mpModelPart->Elements()) {
            r_elem.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids); //FIXME: For the elimination one we should do the AddEntry one by one
        }
        for (auto& r_cond : mpModelPart->Conditions()) {
            r_cond.EquationIdVector(eq_ids, mpModelPart->GetProcessInfo());
            rSparseGraph.AddEntries(eq_ids); //FIXME: For the elimination one we should do the AddEntry one by one
        }
    }

    virtual std::size_t SetUpSystemIds(const DofsArrayType& rDofSet)
    {
        if (mBuildType == BuildType::Block) {
            // Set up the DOFs equation global ids
            IndexPartition<IndexType>(rDofSet.size()).for_each([&](IndexType Index) {
                auto it_dof = rDofSet.begin() + Index;
                it_dof->SetEquationId(Index);
            });

            // Save the size of the system based on the number of DOFs
            mEquationSystemSize = rDofSet.size();

        } else if (mBuildType == BuildType::Elimination) {
            // Set up the DOFs equation global ids
            // The free DOFs are positioned at the beginning of the system
            // The fixed DOFs are positioned at the end of the system (in reversed order)
            IndexType free_id = 0;
            IndexType fix_id = rDofSet.size();
            for (auto it_dof = rDofSet.begin(); it_dof != rDofSet.end(); ++it_dof) {
                if (it_dof->IsFixed()) {
                    it_dof->SetEquationId(--fix_id);
                } else {
                    it_dof->SetEquationId(free_id++);
                }
            }

            // Set the equation system size as the current fixed id
            // Note that this means that if the EquationId is greater than mEquationSystemSize the DOF is fixed
            mEquationSystemSize = fix_id;

        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }

        return mEquationSystemSize;
    }

    virtual void ApplyLinearSystemConstraints(
        const DofsArrayType &rDofArray,
        const EffectiveDofsMapType &rDofIdMap,
        LinearSystemContainer<TSparseMatrixType, TSystemVectorType> &rLinearSystemContainer)
    {
        KRATOS_ERROR << "Calling base class 'ApplyLinearSystemConstraints'." << std::endl;
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        //TODO: Do it as the other assembly functions
        if (mBuildType == BuildType::Block) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Block>(rpLhs, rpEffectiveLhs, rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else if (mBuildType == BuildType::Elimination) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Elimination>(rpLhs, rpEffectiveLhs, rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        //TODO: Do it as the other assembly functions
        if (mBuildType == BuildType::Block) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Block>(rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else if (mBuildType == BuildType::Elimination) {
            ApplyMasterSlaveConstraintsImplementation<BuildType::Elimination>(rpRhs, rEffectiveRhs, rDx, rEffectiveDx, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

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
            const auto& r_process_info = mpModelPart->GetProcessInfo();
            KRATOS_ERROR_IF_NOT(r_process_info.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined in ProcessInfo container. Please set 'BUILD_SCALE_FACTOR' variable." << std::endl;
            return r_process_info.GetValue(BUILD_SCALE_FACTOR);
        } else {
            KRATOS_ERROR << "Wrong scaling type." << std::endl;
        }
    }

    void CalculateSolutionVector(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        const EffectiveDofsMapType& rEffectiveDofIdMap,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector,
        TSystemVectorType& rSolutionVector) const
    {
        // Set an auxiliary vector containing the effective solution values
        const std::size_t n_eff_dofs = rEffectiveDofSet.size();
        TSystemVectorType y(n_eff_dofs);
        IndexPartition<IndexType>(rEffectiveDofSet.size()).for_each([&](IndexType Index) {
            // Get effective DOF
            auto p_dof = *(rEffectiveDofSet.ptr_begin() + Index);
            auto p_dof_find = rEffectiveDofIdMap.find(p_dof);
            KRATOS_ERROR_IF(p_dof_find == rEffectiveDofIdMap.end()) << "DOF cannot be found in DOF id map." << std::endl;

            // Get value from DOF and set it in the auxiliary solution values vector
            // Note that the corresponding row is retrieved from the effective DOF ids map
            y[p_dof_find->second] = p_dof->GetSolutionStepValue();
        });

        // Check solution vector size
        const std::size_t aux_size = rConstraintsConstantVector.size();
        if (rSolutionVector.size() != aux_size) {
            rSolutionVector = TSystemVectorType(aux_size);
        }

        // Initialize solution vector with the constaints constant vector values
        // Note that we deliberately avoid using the copy constructor as we dont want to overwrite the constraints constant vector
        IndexPartition<IndexType>(rSolutionVector.size()).for_each([&](IndexType i){
            rSolutionVector[i] = rConstraintsConstantVector[i];
        });

        // Compute the solution vector as x = T * y + q
        rConstraintsRelationMatrix.SpMV(1.0, y, 1.0, rSolutionVector); // Note that this performs the operation x = A*y + x
    }

    virtual void Clear()
    {
        // Clear the system info
        mEquationSystemSize = 0;
    }

    ///@}
    ///@name Input and output
    ///@{

    std::size_t GetEquationSystemSize() const
    {
        return mEquationSystemSize;
    }

    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr;

    BuildType mBuildType;

    ScalingType mScalingType;

    double mScalingValue;

    std::size_t mEchoLevel;

    std::size_t mEquationSystemSize = 0; /// Number of degrees of freedom of the problem to be solved

    ///@}
    ///@name Private Operations
    ///@{

    template <BuildType TBuildType>
    void ApplyMasterSlaveConstraintsImplementation(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const std::size_t n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSystemVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSystemVectorType(n_master));
            }
            rEffectiveDx.SetValue(0.0);

            // Apply constraints to RHS
            rConstraintsRelationMatrix.TransposeSpMV(*rpRhs, rEffectiveRhs);

            // Apply constraints to LHS
            auto p_LHS_T = AmgclCSRSpMMUtilities::SparseMultiply(*rpLhs, rConstraintsRelationMatrix);
            auto p_transT = AmgclCSRConversionUtilities::Transpose(rConstraintsRelationMatrix);
            rpEffectiveLhs = AmgclCSRSpMMUtilities::SparseMultiply(*p_transT, *p_LHS_T);

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
        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template <BuildType TBuildType>
    void ApplyMasterSlaveConstraintsImplementation(
        typename TSystemVectorType::Pointer& rpRhs,
        TSystemVectorType& rEffectiveRhs,
        TSystemVectorType& rDx,
        TSystemVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSystemVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const std::size_t n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSystemVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSystemVectorType(n_master));
            }
            rEffectiveDx.SetValue(0.0);

            // Apply constraints to RHS
            rConstraintsRelationMatrix.TransposeSpMV(*rpRhs, rEffectiveRhs);

            // // Compute the scale factor value
            // //TODO: think on how to make this user-definable
            // const double scale_factor = rpEffectiveLhs->NormDiagonal();

            // // Apply diagonal values on slave DOFs
            // IndexPartition<IndexType>(mSlaveIds.size()).for_each([&](IndexType Index){
            //     const IndexType slave_eq_id = mSlaveIds[Index];
            //     if (mInactiveSlaveDofs.find(slave_eq_id) == mInactiveSlaveDofs.end()) {
            //         rEffectiveRhs[slave_eq_id] = 0.0;
            //     }
            // });
        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    ///@}
}; // Class Builder

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
