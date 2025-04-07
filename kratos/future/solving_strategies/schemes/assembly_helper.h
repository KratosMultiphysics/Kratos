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
    using DofType = ModelPart::DofType;

    /// DOF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Function type for elements assembly
    using ElementAssemblyFunctionType = std::function<void(ModelPart::ElementConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for conditions assembly
    using ConditionAssemblyFunctionType = std::function<void(ModelPart::ConditionConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for constraints assembly
    using ConstraintAssemblyFunctionType = std::function<void(ModelPart::MasterSlaveConstraintConstantIteratorType, const ProcessInfo&, TThreadLocalStorage&)>;

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
    // TODO: To be discussed in the future. It would be great to have one with references. This will require a Resize method in the sparse arrays that admits a sparse graph.
    virtual void ResizeAndInitializeVectors(
        const DofsArrayType& rDofSet,
        typename TSparseMatrixType::Pointer& rpLHS,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpRHS,
        const bool ReactionVector = false)
    {
        // Set up the sparse matrix graph (note that we do not need to keep it after the resizing)
        TSparseGraphType sparse_graph(mEquationSystemSize);
        SetUpSparseGraph(sparse_graph);

        // Set the system arrays
        // Note that the graph-based constructor does both resizing and initialization
        auto p_lhs = Kratos::make_shared<TSparseMatrixType>(sparse_graph);
        rpLHS.swap(p_lhs);

        auto p_dx = Kratos::make_shared<TSparseVectorType>(sparse_graph);
        rpDx.swap(p_dx);

        auto p_rhs = Kratos::make_shared<TSparseVectorType>(sparse_graph);
        rpRHS.swap(p_rhs);

        //TODO: Maybe we can avoid this if we rebuild the RHS for the reactions calculation --> Check it in the future
        // For the elimination build, also allocate the auxiliary reactions vector
        if (mBuildType == BuildType::Elimination && ReactionVector) {
            KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Equation system size is not set yet. Please call 'SetUpSystemIds' before this method." << std::endl;
            auto p_react = Kratos::make_shared<TSparseVectorType>(rDofSet.size() - mEquationSystemSize);
            mpReactionsVector.swap(p_react);
        }
    }

    virtual void SetUpSparseGraph(TSparseGraphType& rSparseGraph)
    {
        KRATOS_ERROR_IF(mEquationSystemSize == 0) << "Current equation system size is 0. Sparse graph cannot be set (call \'SetUpSystemIds\' first)." << std::endl;

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

    virtual void AssembleLocalSystem(
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

    virtual void AssembleLeftHandSide(
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

    virtual void AssembleRightHandSide(
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

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        if (mBuildType == BuildType::Block) {
            //TODO: Implement this in the CSR matrix or here? --> Most probably we shouldn't call it here neither
            // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
            // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

            ApplyBlockBuildDirichletConditions(rDofSet, rLHS, rRHS);

        } else if (mBuildType == BuildType::Elimination) {

            //TODO: Implement this in the CSR matrix or here? --> Most probably we shouldn't call it here neither
            // // Detect if there is a line of all zeros and set the diagonal to a certain number if this happens (1 if not scale, some norms values otherwise)
            // mScaleFactor = TSparseSpace::CheckAndCorrectZeroDiagonalValues(rModelPart.GetProcessInfo(), rA, rb, mScalingDiagonal);

        } else {
            KRATOS_ERROR << "Build type not supported." << std::endl;
        }
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS)
    {
        if (mBuildType == BuildType::Block) {
            ApplyBlockBuildDirichletConditions(rDofSet, rRHS);
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

    SizeType mEquationSystemSize = 0; /// Number of degrees of freedom of the problem to be solve

    std::unique_ptr<ElementAssemblyFunctionType> mpElementAssemblyFunction = nullptr; // Pointer to the function to be called in the elements assembly

    std::unique_ptr<ConditionAssemblyFunctionType> mpConditionAssemblyFunction = nullptr; // Pointer to the function to be called in the conditions assembly

    std::unique_ptr<ConstraintAssemblyFunctionType> mpConstraintAssemblyFunction = nullptr; // Pointer to the function to be called in the constraints assembly

    typename TSparseVectorType::Pointer mpReactionsVector = nullptr; // Auxiliary vector to calculate the reactions in the elimination build

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
                    (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local LHS and RHS contributions
                    auto it_cond = conds_begin + k;
                    (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
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
                    (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType>(rTLS, rLHS);
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local LHS contributions
                    auto it_cond = conds_begin + k;
                    (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType>(rTLS, rLHS);
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
                    (*mpElementAssemblyFunction)(it_elem, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType, TAssembleReactionVector>(rTLS, rRHS);
                }
            }

            // Assemble conditions
            if (mpConditionAssemblyFunction != nullptr) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    // Calculate local RHS contributions
                    auto it_cond = conds_begin + k;
                    (*mpConditionAssemblyFunction)(it_cond, r_process_info, rTLS);

                    // Assemble the local contributions to the global system
                    AssembleLocalContribution<TBuildType, TAssembleReactionVector>(rTLS, rRHS);
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

        if (r_loc_eq_ids.size() != 0) { // Note that inactive elements TLS is resized to zero
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
                static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
            }
        }
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        TSparseMatrixType& rLHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);

        if (r_loc_eq_ids.size() != 0) { // Note that inactive elements TLS is resized to zero
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
                static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
            }
        }
    }

    template<BuildType TBuildType, bool TAssembleReactionVector>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        TSparseVectorType& rRHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);

        if (r_loc_eq_ids.size() != 0) { // Note that inactive elements TLS is resized to zero
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
                static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
            }
        }
    }

    template<class TContainerType>
    auto& GetThreadLocalStorageContainer(
        const TContainerType& rContainer,
        const TThreadLocalStorage& rTLS)
    {
        if constexpr (std::is_same_v<TContainerType, TSparseMatrixType>) {
            return rTLS.LocalLhs; // We can eventually do a get method in the TLS and call it in here
        } else if (std::is_same_v<TContainerType, TSparseVectorType>) {
            return rTLS.LocalRhs; // We can eventually do a get method in the TLS and call it in here
        } else {
            static_assert(std::is_same_v<TContainerType, TSparseMatrixType> || std::is_same_v<TContainerType, TSparseVectorType>, "Unssupported container type.");
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

    void ApplyBlockBuildDirichletConditions(
        const DofsArrayType& rDofSet,
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
        const auto dof_begin = rDofSet.begin();
        IndexPartition<std::size_t>(rDofSet.size()).for_each([&](IndexType Index){
            auto it_dof = dof_begin + Index;
            if (it_dof->IsFixed()) {
                free_dofs_vector[Index] = 0;
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
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS) const
    {
        // Loop the DOFs to find which ones are fixed
        // Note that DOFs are assumed to be numbered consecutively in the block building
        const auto dof_begin = rDofSet.begin();
        IndexPartition<std::size_t>(rDofSet.size()).for_each([&](IndexType Index){
            auto it_dof = dof_begin + Index;
            if (it_dof->IsFixed()) {
                rRHS[Index] = 0.0;
            }
        });
    }

    ///@}
}; // Class AssemblyHelper

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
