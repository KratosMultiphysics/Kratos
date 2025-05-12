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

/**
 * @class AssemblyHelper
 * @ingroup KratosCore
 * @brief Utility class for handling the build
 * @details This class collects the methods required to
 * handle the building of the sparse sytem matrices
 * @author Ruben Zorrilla
 */
template<class TThreadLocalStorage, class TSparseMatrixType, class TSparseVectorType, class TSparseGraphType>
class AssemblyHelper
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

    /// Pointer definition of AssemblyHelper
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyHelper);

    /// Size type definition
    using SizeType = std::size_t;

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

    /// Function type for elements assembly
    using ElementAssemblyFunctionType = std::function<bool(ModelPart::ElementConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for conditions assembly
    using ConditionAssemblyFunctionType = std::function<bool(ModelPart::ConditionConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for constraints assembly
    using ConstraintAssemblyFunctionType = std::function<bool(ModelPart::MasterSlaveConstraintConstantIteratorType, const ProcessInfo&, TThreadLocalStorage&)>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyHelper() = delete;

    /// Constructor with model part
    AssemblyHelper(
        const ModelPart& rModelPart,
        Parameters AssemblySettings = Parameters(R"({})"))
    : mpModelPart(&rModelPart)
    {
        Parameters default_parameters( R"({
            "build_type" : "block",
            "scaling_type" : "max_diagonal",
            "echo_level" : 0
        })");
        AssemblySettings.ValidateAndAssignDefaults(default_parameters);

        // Set build type
        const std::string build_type = AssemblySettings["build_type"].GetString();
        if (build_type == "block") {
            mBuildType = BuildType::Block;
        } else if (build_type == "elimination") {
            mBuildType = BuildType::Elimination;
        } else {
            KRATOS_ERROR << "Provided 'build_type' is '" << build_type << "'. Available options are:\n"
            << "\t- 'block'\n"
            << "\t- 'elimination'" << std::endl;
        }

        // Set scaling type
        const std::string scaling_type = AssemblySettings["scaling_type"].GetString();
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
        mEchoLevel = AssemblySettings["echo_level"].GetInt();
    }

    virtual ~AssemblyHelper() = default;

    ///@}
    ///@name Operations
    ///@{

    //TODO: In the future, add the one with
    // - LHS alone
    // - MassMatrix
    // - DampingMatrix
    // TODO: To be discussed in the future. It would be great to have one with references. This will require using the move constructors
    virtual void ResizeAndInitializeVectors(
        const DofsArrayType& rDofSet,
        const DofsArrayType& rEffectiveDofSet,
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpRhs,
        typename TSparseVectorType::Pointer& rpEffectiveRhs,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpEffectiveDx,
        const bool ReactionVector = false)
    {
        // Set up the sparse matrix graph (note that we do not need to keep it after the resizing)
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
        TSparseGraphType sparse_graph(mEquationSystemSize);
        SetUpSparseGraph(sparse_graph);

        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_lhs = Kratos::make_shared<TSparseMatrixType>(sparse_graph);
        rpLhs.swap(p_lhs);

        auto p_dx = Kratos::make_shared<TSparseVectorType>(sparse_graph);
        rpDx.swap(p_dx);

        auto p_rhs = Kratos::make_shared<TSparseVectorType>(sparse_graph);
        rpRhs.swap(p_rhs);

        // Set the effective arrays
        // In a standard case we only need to allocate the effective solution update to avoid the first Predict() call to crash
        if (rDofSet == rEffectiveDofSet) {
            // If there are no constraints, the effective DOF set matches the standard one and the effective arrays are the same as the input ones
            // Note that we avoid duplicating the memory by making the effective pointers to point to the same object
            rpEffectiveDx = rpDx;
        } else {
            // Note that we only need to allocate the effective solution vector as the effective LHS and RHS would be allocated when applying the constraints
            auto p_eff_dx = Kratos::make_shared<TSparseVectorType>(rEffectiveDofSet.size());
            rpEffectiveDx.swap(p_eff_dx);
        }

        //TODO: Maybe we can avoid this if we rebuild the RHS for the reactions calculation --> Check it in the future
        // For the elimination build, also allocate the auxiliary reactions vector
        if (mBuildType == BuildType::Elimination && ReactionVector) {
            auto p_react = Kratos::make_shared<TSparseVectorType>(rDofSet.size() - mEquationSystemSize);
            mpReactionsVector.swap(p_react);
        }
    }

    virtual void ConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        if (mBuildType == BuildType::Block) {
            BlockConstructMasterSlaveConstraintsStructure(rModelPart, rDofSet, rEffectiveDofSet, rEffectiveDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else if (mBuildType == BuildType::Elimination) {
            EliminationConstructMasterSlaveConstraintsStructure(rModelPart, rDofSet, rEffectiveDofSet, rEffectiveDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void SetUpSparseGraph(TSparseGraphType& rSparseGraph)
    {
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Current equation system size is 0. Sparse graph cannot be set (call \'SetUpSystemIds\' first)." << std::endl;

        //TODO: This is always thrown as the mEquationSystemSize is already in the graph. I'd use the move (not implemented constructor)
        if (!rSparseGraph.IsEmpty()) {
            KRATOS_WARNING("AssemblyHelper") << "Provided sparse graph is not empty and will be cleared." << std::endl;
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

    virtual SizeType SetUpSystemIds(const DofsArrayType& rDofSet)
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

    virtual void Assemble(
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Call the implementation of the function with building type template argument
        if (mBuildType == BuildType::Block) {
            AssembleImplementation<BuildType::Block>(rLHS, rRHS, rTLS);
        } else if (mBuildType == BuildType::Elimination) {
            AssembleImplementation<BuildType::Elimination>(rLHS, rRHS, rTLS);
        } else {
            KRATOS_ERROR << "Not implemented build type." << std::endl;
        }
    }

    virtual void Assemble(
        TSparseMatrixType& rLHS,
        TThreadLocalStorage& rTLS)
    {
        // Call the implementation of the function with building type template argument
        if (mBuildType == BuildType::Block) {
            AssembleImplementation<BuildType::Block>(rLHS, rTLS);
        } else if (mBuildType == BuildType::Elimination) {
            AssembleImplementation<BuildType::Elimination>(rLHS, rTLS);
        } else {
            KRATOS_ERROR << "Not implemented build type." << std::endl;
        }
    }

    virtual void Assemble(
        TSparseVectorType& rRHS,
        TThreadLocalStorage& rTLS,
        const bool AssembleReactionVector = false)
    {
        // Call the implementation of the function with building type template argument
        if (mBuildType == BuildType::Block) {
            AssembleImplementation<BuildType::Block, false>(rRHS, rTLS);
        } else if (mBuildType == BuildType::Elimination) {
            if (AssembleReactionVector) {
                AssembleImplementation<BuildType::Elimination, true>(rRHS, rTLS);
            } else {
                AssembleImplementation<BuildType::Elimination, false>(rRHS, rTLS);
            }
        } else {
            KRATOS_ERROR << "Not implemented build type." << std::endl;
        }
    }

    virtual void AssembleMasterSlaveConstraints(
        const DofsArrayType& rDofSet,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector,
        TThreadLocalStorage& rTLS)
    {
        //TODO: Do it as the other assembly functions
        if (mBuildType == BuildType::Block) {
            AssembleMasterSlaveConstraintsImplementation<BuildType::Block>(rDofSet, rDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector, rTLS);
        } else if (mBuildType == BuildType::Elimination) {
            AssembleMasterSlaveConstraintsImplementation<BuildType::Elimination>(rDofSet, rDofIdMap, rConstraintsRelationMatrix, rConstraintsConstantVector, rTLS);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void ApplyMasterSlaveConstraints(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpRhs,
        TSparseVectorType& rEffectiveRhs,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
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
        typename TSparseVectorType::Pointer& rpRhs,
        TSparseVectorType& rEffectiveRhs,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
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

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        if (mBuildType == BuildType::Block) {
            //TODO: Implement this in the CSR matrix or here? --> Most probably we shouldn't call it here neither
            // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
            // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);
            if (rDofIdMap.empty()) {
                ApplyBlockBuildDirichletConditions(rDofArray, rLHS, rRHS);
            } else {
                ApplyBlockBuildDirichletConditions(rDofArray, rDofIdMap, rLHS, rRHS);
            }
        } else if (mBuildType == BuildType::Elimination) {
            //TODO: Implement this in the CSR matrix or here? --> Most probably we shouldn't call it here neither
            // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
            // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseVectorType& rRHS)
    {
        if (mBuildType == BuildType::Block) {
            if (rDofIdMap.empty()) {
                ApplyBlockBuildDirichletConditions(rDofArray, rRHS);
            } else {
                ApplyBlockBuildDirichletConditions(rDofArray, rDofIdMap, rRHS);
            }
        } else if (mBuildType == BuildType::Elimination) {
            return;
        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void CalculateReactionsRightHandSide(
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Initialize the provided RHS (note that this has been potentially used in the system resolution)
        rRHS.SetValue(DataType());

        if (mBuildType == BuildType::Block) {
            // Do the block RHS assembly without Dirichlet BCs
            AssembleImplementation<BuildType::Block, false>(rRHS, rTLS);

            // Set minus the RHS as reaction values
            // Note that in the block build DOFs are assumed to be numbered consecutively
            block_for_each(rDofSet, [&](DofType& rDof){
                rDof.GetSolutionStepReactionValue() = -rRHS[rDof.EquationId()];
            });
        } else if (mBuildType == BuildType::Elimination) {
            //FIXME: This is wrong
            // Do the elimination RHS assembly without Dirichlet BCs
            AssembleImplementation<BuildType::Elimination, true>(rRHS, rTLS);

            // Set minus the RHS as reaction values
            // Note that in the elimination case the fix DOFs residuals are stored in the mpReactionsVector
            auto& r_reactions_vector = *mpReactionsVector;
            block_for_each(rDofSet, [&](DofType& rDof){
                IndexType i_react = rDof.EquationId();
                if (i_react >= mEquationSystemSize) { // Check if current DOF is fixed
                    i_react -= mEquationSystemSize; // Get the corresponding row index in the reactions vector
                    rDof.GetSolutionStepReactionValue() = - r_reactions_vector[i_react];
                }
            });
            //FIXME: This is wrong
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

    virtual void Clear()
    {
        // Clear the system info
        mEquationSystemSize = 0;

        // Clear the assembly functions
        mpElementAssemblyFunction = nullptr;
        mpConditionAssemblyFunction = nullptr;
        mpConstraintAssemblyFunction = nullptr;
    }

    ///@}
    ///@name Input and output
    ///@{

    SizeType GetEquationSystemSize() const
    {
        return mEquationSystemSize;
    }

    void SetElementAssemblyFunction(ElementAssemblyFunctionType rElementAssemblyFunction)
    {
        auto p_aux = std::make_unique<ElementAssemblyFunctionType>(rElementAssemblyFunction);
        mpElementAssemblyFunction.swap(p_aux);
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Element assembly function set." << std::endl;
    }

    void SetConditionAssemblyFunction(ConditionAssemblyFunctionType rConditionAssemblyFunction)
    {
        auto p_aux = std::make_unique<ConditionAssemblyFunctionType>(rConditionAssemblyFunction);
        mpConditionAssemblyFunction.swap(p_aux);
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Condition assembly function set." << std::endl;
    }

    //TODO: For sure we'll need to modify this one
    void SetConstraintAssemblyFunction(ConstraintAssemblyFunctionType rConstraintAssemblyFunction)
    {
        auto p_aux = std::make_unique<ConstraintAssemblyFunctionType>(rConstraintAssemblyFunction);
        mpConstraintAssemblyFunction.swap(p_aux);
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Constraint assembly function set." << std::endl;
    }

    ///@}
private:
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr;

    BuildType mBuildType;

    ScalingType mScalingType;

    double mScalingValue;

    SizeType mEchoLevel;

    SizeType mEquationSystemSize = 0; /// Number of degrees of freedom of the problem to be solved

    //TODO: study performance if we get rid of the pointer to the assembly function (pass them by argument)
    std::unique_ptr<ElementAssemblyFunctionType> mpElementAssemblyFunction = nullptr; // Pointer to the function to be called in the elements assembly

    //TODO: study performance if we get rid of the pointer to the assembly function (pass them by argument)
    std::unique_ptr<ConditionAssemblyFunctionType> mpConditionAssemblyFunction = nullptr; // Pointer to the function to be called in the conditions assembly

    //TODO: study performance if we get rid of the pointer to the assembly function (pass them by argument)
    std::unique_ptr<ConstraintAssemblyFunctionType> mpConstraintAssemblyFunction = nullptr; // Pointer to the function to be called in the constraints assembly

    typename TSparseVectorType::Pointer mpReactionsVector = nullptr; // Auxiliary vector to calculate the reactions in the elimination build

    // std::vector<IndexType> mSlaveIds; /// Vector containing the equation ids of the slaves

    // std::vector<IndexType> mMasterIds; /// Vector containing the equation ids of the master

    std::unordered_set<IndexType> mInactiveSlaveDofs; /// The set containing the inactive slave DOFs (only used in the block build)

    ///@}
    ///@name Private Operations
    ///@{

    template<BuildType TBuildType>
    void AssembleImplementation(
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const SizeType n_elems = r_elems.size();
        const SizeType n_conds = r_conds.size();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.BeginAssemble();

        // Assemble entities
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, rTLS)
        {
            // Assemble elements
            if (mpElementAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    // Calculate local LHS and RHS contributions
                    auto it_elem = elems_begin + k;
                    const bool assemble = (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
                    }
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local LHS and RHS contributions
                    auto it_cond = conds_begin + k;
                    const bool assemble = (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
                    }
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.FinalizeAssemble();
        rLHS.FinalizeAssemble();
    }

    template<BuildType TBuildType>
    void AssembleImplementation(
        TSparseMatrixType& rLHS,
        TThreadLocalStorage& rTLS)
    {
        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const SizeType n_elems = r_elems.size();
        const SizeType n_conds = r_conds.size();

        // Initialize LHS assembly
        rLHS.BeginAssemble();

        // Assemble entities
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, rTLS)
        {
            // Assemble elements
            if (mpElementAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    // Calculate local LHS contributions
                    auto it_elem = elems_begin + k;
                    const bool assemble = (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS);
                    }
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local LHS contributions
                    auto it_cond = conds_begin + k;
                    const bool assemble = (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS);
                    }
                }
            }
        }

        // Finalize LHS assembly
        rLHS.FinalizeAssemble();
    }

    template<BuildType TBuildType, bool TAssembleReactionVector>
    void AssembleImplementation(
        TSparseVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Getting conditions and elements to be assembled
        const auto& r_elems = mpModelPart->Elements();
        const auto& r_conds = mpModelPart->Conditions();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting entities container data
        auto elems_begin = r_elems.begin();
        auto conds_begin = r_conds.begin();
        const SizeType n_elems = r_elems.size();
        const SizeType n_conds = r_conds.size();

        // Initialize RHS and LHS assembly
        rRHS.BeginAssemble();

        // Assemble entities
        #pragma omp parallel firstprivate(n_elems, n_conds, elems_begin, conds_begin, r_process_info, rTLS)
        {
            // Assemble elements
            if (mpElementAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    // Calculate local RHS contributions
                    auto it_elem = elems_begin + k;
                    const bool assemble = (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType, TAssembleReactionVector>(rTLS, rRHS);
                    }
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local RHS contributions
                    auto it_cond = conds_begin + k;
                    const bool assemble = (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    if (assemble) {
                        AssembleLocalContribution<TBuildType, TAssembleReactionVector>(rTLS, rRHS);
                    }
                }
            }
        }

        // Finalize RHS assembly
        rRHS.FinalizeAssemble();
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
        if constexpr (TBuildType == BuildType::Block) {
            rRHS.Assemble(GetThreadLocalStorageContainer(rRHS, rTLS), r_loc_eq_ids); // RHS contributions assembly
            rLHS.Assemble(GetThreadLocalStorageContainer(rLHS, rTLS), r_loc_eq_ids); // LHS contributions assembly

        } else if (TBuildType == BuildType::Elimination) {
            const auto& r_loc_rhs = GetThreadLocalStorageContainer(rRHS, rTLS);
            const auto& r_loc_lhs = GetThreadLocalStorageContainer(rLHS, rTLS);
            const SizeType loc_size = r_loc_rhs.size();

            for (IndexType i_loc = 0; i_loc < loc_size; ++i_loc) {
                IndexType i_glob = r_loc_eq_ids[i_loc];
                if (i_glob < mEquationSystemSize) {// Check if current row DOF is free
                    rRHS.AssembleEntry(r_loc_rhs[i_loc], i_glob); // RHS contribution assembly
                    for (IndexType j_loc = 0; j_loc < loc_size; ++j_loc) {
                        const IndexType j_glob = r_loc_eq_ids[j_loc];
                        if (j_glob < mEquationSystemSize) {// Check if current column DOF is free
                            rLHS.AssembleEntry(r_loc_lhs(i_loc, j_loc), i_glob, j_glob); // LHS contribution assembly
                        }
                    }
                }
            }

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        TSparseMatrixType& rLHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
        if constexpr (TBuildType == BuildType::Block) {
            rLHS.Assemble(GetThreadLocalStorageContainer(rLHS, rTLS), r_loc_eq_ids); // LHS contributions assembly

        } else if (TBuildType == BuildType::Elimination) {
            const auto& r_loc_lhs = GetThreadLocalStorageContainer(rLHS, rTLS);
            const SizeType loc_size = r_loc_lhs.size1();

            for (IndexType i_loc = 0; i_loc < loc_size; ++i_loc) {
                IndexType i_glob = r_loc_eq_ids[i_loc];
                if (i_glob < mEquationSystemSize) {// Check if current row DOF is free
                    for (IndexType j_loc = 0; j_loc < loc_size; ++j_loc) {
                        const IndexType j_glob = r_loc_eq_ids[j_loc];
                        if (j_glob < mEquationSystemSize) {// Check if current column DOF is free
                            rLHS.AssembleEntry(r_loc_lhs(i_loc, j_loc), i_glob, j_glob); // LHS contribution assembly
                        }
                    }
                }
            }

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template<BuildType TBuildType, bool TAssembleReactionVector>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        TSparseVectorType& rRHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
        if constexpr (TBuildType == BuildType::Block) {
            if constexpr (!TAssembleReactionVector) {
                rRHS.Assemble(GetThreadLocalStorageContainer(rRHS, rTLS), r_loc_eq_ids); // RHS contributions assembly
            } else {
                static_assert(TBuildType == BuildType::Block && TAssembleReactionVector == true, "Assemble reaction vector cannot be used with block build type.");
            }

        } else if (TBuildType == BuildType::Elimination) {
            const auto& r_loc_rhs = GetThreadLocalStorageContainer(rRHS, rTLS);
            const SizeType loc_size = r_loc_rhs.size();
            for (IndexType i_loc = 0; i_loc < loc_size; ++i_loc) {
                IndexType i_glob = r_loc_eq_ids[i_loc];
                if constexpr (!TAssembleReactionVector) {
                    if (i_glob < mEquationSystemSize) {// Check if current row DOF is free
                        rRHS.AssembleEntry(r_loc_rhs[i_loc], i_glob); // RHS contribution assembly
                    }
                } else {
                    if (i_glob < mEquationSystemSize) {// Check if current row DOF is free
                        rRHS.AssembleEntry(r_loc_rhs[i_loc], i_glob); // RHS contribution assembly
                    } else {
                        const IndexType react_vec_pos = i_glob - mEquationSystemSize;// Get the corresponding position in the reactions vector
                        mpReactionsVector->AssembleEntry(r_loc_rhs[i_loc], react_vec_pos); // RHS contribution assembly to reactions vector
                    }
                }
            }

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template<class TContainerType>
    auto& GetThreadLocalStorageContainer(
        const TContainerType& rContainer,
        const TThreadLocalStorage& rTLS)
    {
        if constexpr (std::is_same_v<TContainerType, TSparseMatrixType>) {
            return rTLS.LocalMatrix; // We can eventually do a get method in the TLS and call it in here
        } else if (std::is_same_v<TContainerType, TSparseVectorType>) {
            return rTLS.LocalVector; // We can eventually do a get method in the TLS and call it in here
        } else {
            static_assert(std::is_same_v<TContainerType, TSparseMatrixType> || std::is_same_v<TContainerType, TSparseVectorType>, "Unsupported container type.");
        }
    }

    auto& GetThreadLocalStorageEqIds(TThreadLocalStorage& rTLS)
    {
        return rTLS.LocalEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    auto& GetThreadLocalStorageEqIds(const TThreadLocalStorage& rTLS)
    {
        return rTLS.LocalEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    auto& GetThreadLocalStorageSlaveEqIds(TThreadLocalStorage& rTLS)
    {
        return rTLS.SlaveEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    auto& GetThreadLocalStorageSlaveEqIds(const TThreadLocalStorage& rTLS)
    {
        return rTLS.SlaveEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    auto& GetThreadLocalStorageMasterEqIds(TThreadLocalStorage& rTLS)
    {
        return rTLS.MasterEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    auto& GetThreadLocalStorageMasterEqIds(const TThreadLocalStorage& rTLS)
    {
        return rTLS.MasterEqIds; // We can eventually do a get method in the TLS and call it in here
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS) const
    {
        // Set the free DOFs vector (0 means fixed / 1 means free)
        // Note that we initialize to 1 so we start assuming all free
        // Also note that the type is uint_8 for the sake of efficiency
        const SizeType system_size = rLHS.size1();
        std::vector<uint8_t> free_dofs_vector(system_size, 1);

        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_ptr_begin = rDofArray.ptr_begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            const auto p_dof = *(dof_ptr_begin + Index);
            if (p_dof->IsFixed()) {
                const auto p_dof_find = rDofIdMap.find(p_dof);
                KRATOS_ERROR_IF(p_dof_find == rDofIdMap.end()) << "DOF cannot be found in DOF id map." << std::endl;
                free_dofs_vector[p_dof_find->second] = 0;
            }
        });

        //TODO: Implement this in the CSR matrix or here?
        // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
        // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        // Get the diagonal scaling factor
        const double diagonal_value = GetDiagonalScalingFactor(rLHS);

        // Apply the free DOFs (i.e., fixity) vector to the system arrays
        rLHS.ApplyHomogeneousDirichlet(free_dofs_vector, diagonal_value, rRHS);
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS) const
    {
        // Set the free DOFs vector (0 means fixed / 1 means free)
        // Note that we initialize to 1 so we start assuming all free
        // Also note that the type is uint_8 for the sake of efficiency
        const SizeType system_size = rLHS.size1();
        std::vector<uint8_t> free_dofs_vector(system_size, 1);

        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_begin = rDofArray.begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            const auto p_dof = dof_begin + Index;
            if (p_dof->IsFixed()) {
                free_dofs_vector[p_dof->EquationId()] = 0;
            }
        });

        //TODO: Implement this in the CSR matrix or here?
        // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
        // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        // Get the diagonal scaling factor
        const double diagonal_value = GetDiagonalScalingFactor(rLHS);

        // Apply the free DOFs (i.e., fixity) vector to the system arrays
        rLHS.ApplyHomogeneousDirichlet(free_dofs_vector, diagonal_value, rRHS);
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseVectorType& rRHS) const
    {
        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_ptr_begin = rDofArray.ptr_begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            auto p_dof = *(dof_ptr_begin + Index);
            if (p_dof->IsFixed()) {
                auto p_dof_find = rDofIdMap.find(p_dof);
                KRATOS_ERROR_IF(p_dof_find == rDofIdMap.end()) << "DOF cannot be found in DOF id map." << std::endl;
                rRHS[p_dof_find->second] = 0;
            }
        });
    }

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofArray,
        TSparseVectorType& rRHS) const
    {
        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_begin = rDofArray.begin();
        IndexPartition<std::size_t>(rDofArray.size()).for_each([&](IndexType Index){
            auto p_dof = dof_begin + Index;
            if (p_dof->IsFixed()) {
                rRHS[p_dof->EquationId()] = 0;
            }
        });
    }

    //FIXME: This only works in serial!!!
    void BlockConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        // Clear the provided effective DOFs map
        KRATOS_WARNING_IF("AssemblyHelper", !rEffectiveDofSet.empty()) << "Provided effective DOFs set is not empty. About to clear it." << std::endl;
        KRATOS_WARNING_IF("AssemblyHelper", !rEffectiveDofIdMap.empty()) << "Provided effective DOFs ids map is not empty. About to clear it." << std::endl;
        rEffectiveDofSet.clear();
        rEffectiveDofIdMap.clear();

        //FIXME: Do the IsActiveConstraints in here and set a flag that stays "forever"

        // Check if there are constraints to build the effective DOFs map and the corresponding arrays
        const SizeType n_constraints = rModelPart.NumberOfMasterSlaveConstraints();
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
            KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
            TSparseGraphType constraints_sparse_graph(mEquationSystemSize);

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
            rConstraintsConstantVector = std::move(TSparseVectorType(mEquationSystemSize));
            rConstraintsRelationMatrix = std::move(TSparseMatrixType(constraints_sparse_graph));

        } else {
            rEffectiveDofSet = rDofSet; // If there are no constraints the effective DOF set is the standard one
            rEffectiveDofIdMap = EffectiveDofsMapType(); // Create an empty master ids map as the standard DOF equation ids can be used
        }
    }

    void EliminationConstructMasterSlaveConstraintsStructure(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofSet,
        DofsArrayType& rEffectiveDofSet,
        EffectiveDofsMapType& rEffectiveDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector)
    {
        KRATOS_ERROR << "Not implemented yet." << std::endl;
    }

    template<BuildType TBuildType>
    void AssembleMasterSlaveConstraintsImplementation(
        const DofsArrayType& rDofSet,
        const EffectiveDofsMapType& rDofIdMap,
        TSparseMatrixType& rConstraintsRelationMatrix,
        TSparseVectorType& rConstraintsConstantVector,
        TThreadLocalStorage& rTLS)
    {
        // Getting constraints to be assembled
        const auto& r_consts = mpModelPart->MasterSlaveConstraints();
        const auto& r_process_info = mpModelPart->GetProcessInfo();

        // Getting constraints container data
        auto consts_begin = r_consts.begin();
        const SizeType n_consts = r_consts.size();

        // Initialize constraints arrays
        rConstraintsRelationMatrix.SetValue(0.0);
        rConstraintsConstantVector.SetValue(0.0);

        if constexpr (TBuildType == BuildType::Block) {
            // We clear the inactive DOFs set
            mInactiveSlaveDofs.clear();
        }

        rConstraintsRelationMatrix.BeginAssemble();
        rConstraintsConstantVector.BeginAssemble();

        #pragma omp parallel firstprivate(rDofIdMap, consts_begin, r_process_info)
        {
            // Auxiliary set to store the inactive constraints slave DOFs (required by the block build)
            std::unordered_set<IndexType> auxiliar_inactive_slave_dofs;

            // Assemble constraints
            if (mpConstraintAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_consts; ++k) {
                    // Calculate local contributions
                    auto it_const = consts_begin + k;
                    const bool assemble_const = (*mpConstraintAssemblyFunction)(it_const, r_process_info, rTLS);

                    // Set the master and slave equation ids
                    // Note that the slaves follow the system equation ids while the masters use the effective map ones
                    const auto& r_slave_dofs = it_const->GetSlaveDofsVector();
                    auto& r_slave_eq_ids = GetThreadLocalStorageSlaveEqIds(rTLS);
                    const SizeType n_slaves = r_slave_dofs.size();
                    if (r_slave_eq_ids.size() != n_slaves) {
                        r_slave_eq_ids.resize(n_slaves);
                    }
                    for (IndexType i_slave = 0; i_slave < n_slaves; ++i_slave) {
                        r_slave_eq_ids[i_slave] = (*(r_slave_dofs.begin() + i_slave))->EquationId();
                    }

                    const auto& r_master_dofs = it_const->GetMasterDofsVector();
                    auto& r_master_eq_ids = GetThreadLocalStorageMasterEqIds(rTLS);
                    const SizeType n_masters = r_master_dofs.size();
                    if (r_master_eq_ids.size() != n_masters) {
                        r_master_eq_ids.resize(n_masters);
                    }
                    for (IndexType i_master = 0; i_master < n_masters; ++i_master) {
                        auto p_master = *(r_master_dofs.begin() + i_master);
                        auto p_master_find = rDofIdMap.find(p_master);
                        KRATOS_ERROR_IF(p_master_find == rDofIdMap.end()) << "Master DOF cannot be found in DOF ids map." << std::endl;
                        r_master_eq_ids[i_master] = p_master_find->second;
                    }

                    // Assemble the constraints local contributions to the global system
                    if (assemble_const) {
                        if constexpr (TBuildType == BuildType::Block) {
                            // Assemble relation matrix contribution
                            const auto& r_loc_T = GetThreadLocalStorageContainer(rConstraintsRelationMatrix, rTLS);
                            rConstraintsRelationMatrix.Assemble(r_loc_T, r_slave_eq_ids, r_master_eq_ids);

                            // Assemble constant vector contribution
                            const auto& r_loc_v = GetThreadLocalStorageContainer(rConstraintsConstantVector, rTLS);
                            rConstraintsConstantVector.Assemble(r_loc_v, r_slave_eq_ids);
                        } else if (TBuildType == BuildType::Elimination) {

                        } else {
                            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
                        }
                    } else {
                        if constexpr (TBuildType == BuildType::Block) {
                            auxiliar_inactive_slave_dofs.insert(r_slave_eq_ids.begin(), r_slave_eq_ids.end());
                        } else if (TBuildType == BuildType::Elimination) {

                        } else {
                            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
                        }
                    }
                }

                // We merge all the sets in one thread
                if constexpr (TBuildType == BuildType::Block) {
                    #pragma omp critical
                    {
                        mInactiveSlaveDofs.insert(auxiliar_inactive_slave_dofs.begin(), auxiliar_inactive_slave_dofs.end());
                    }
                }
            }
        }

        rConstraintsRelationMatrix.FinalizeAssemble();
        rConstraintsConstantVector.FinalizeAssemble();

        if constexpr (TBuildType == BuildType::Block) {
            // Setting the missing effective but not constrain-related DOFs into the T and C system
            // For doing so we loop the standard DOF array (the one from elements and conditions)
            // We search for each DOF in the effective DOF ids map, if present it means its effective
            IndexPartition<IndexType>(rDofSet.size()).for_each([&](IndexType Index){
                const auto p_dof = *(rDofSet.ptr_begin() + Index);
                const auto p_dof_find = rDofIdMap.find(p_dof);
                if (p_dof_find != rDofIdMap.end()) {
                    rConstraintsConstantVector[p_dof->EquationId()] = 0.0;
                    rConstraintsRelationMatrix(p_dof->EquationId(), p_dof_find->second) = 1.0;
                }
            });

            // Setting inactive slave dofs in the T and C system
            //TODO: Can't this be parallel?
            for (auto eq_id : mInactiveSlaveDofs) {
                rConstraintsConstantVector[eq_id] = 0.0;
                rConstraintsRelationMatrix(eq_id, eq_id) = 1.0;
            }
        } else if (TBuildType == BuildType::Elimination) {

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unsupported build type.");
        }
    }

    template <BuildType TBuildType>
    void ApplyMasterSlaveConstraintsImplementation(
        typename TSparseMatrixType::Pointer& rpLhs,
        typename TSparseMatrixType::Pointer& rpEffectiveLhs,
        typename TSparseVectorType::Pointer& rpRhs,
        TSparseVectorType& rEffectiveRhs,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const SizeType n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSparseVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSparseVectorType(n_master));
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
        typename TSparseVectorType::Pointer& rpRhs,
        TSparseVectorType& rEffectiveRhs,
        TSparseVectorType& rDx,
        TSparseVectorType& rEffectiveDx,
        const TSparseMatrixType& rConstraintsRelationMatrix,
        const TSparseVectorType& rConstraintsConstantVector)
    {
        if constexpr (TBuildType == BuildType::Block) {
            // Get the effective size as the number of master DOFs
            // Note that this matches the number of columns in the constraints relation matrix
            const SizeType n_master = rConstraintsRelationMatrix.size2();

            // Initialize the effective RHS
            if (rEffectiveRhs.size() != n_master) {
                rEffectiveRhs = std::move(TSparseVectorType(n_master));
            }
            rEffectiveRhs.SetValue(0.0);

            // Initialize the effective solution vector
            if (rEffectiveDx.size() != n_master) {
                rEffectiveDx = std::move(TSparseVectorType(n_master));
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
}; // Class AssemblyHelper

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
