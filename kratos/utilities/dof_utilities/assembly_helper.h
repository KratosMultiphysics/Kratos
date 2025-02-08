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
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/dof_utilities/assembly_utilities.h"

namespace Kratos
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
template<class TThreadLocalStorage, class TMatrixType=CsrMatrix<>, class TVectorType=SystemVector<>>
class AssemblyHelper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssemblyHelper
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyHelper);

    /// Size type definition
    using SizeType = std::size_t;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Data type definition from sparse matrix
    using DataType = typename SystemMatrixType::DataType;

    /// Index type definition from sparse matrix
    using IndexType = typename SystemMatrixType::IndexType;

    /// Dense space definition
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;

    /// Local system matrix type definition
    using LocalSystemMatrixType = DenseSpaceType::MatrixType;

    /// Local system vector type definition
    using LocalSystemVectorType = DenseSpaceType::VectorType;

    /// Function type for elements assembly
    using ElementAssemblyFunctionType = std::function<void(ModelPart::ElementConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for conditions assembly
    using ConditionAssemblyFunctionType = std::function<void(ModelPart::ConditionConstantIterator, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Function type for constraints assembly
    using ConstraintAssemblyFunctionType = std::function<void(ModelPart::MasterSlaveConstraintConstantIteratorType, const ProcessInfo&, TThreadLocalStorage&)>;

    /// Build enum type definition
    using BuildType = AssemblyUtilities::BuildType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyHelper() = delete;

    /// Constructor with model part
    AssemblyHelper(
        const ModelPart& rModelPart,
        const SizeType EquationSystemSize,
        Parameters AssemblySettings = Parameters(R"({})"))
    : mpModelPart(&rModelPart)
    , mEquationSystemSize(EquationSystemSize)
    {
        Parameters default_parameters( R"({
            "build_type" : "block",
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

        // Set verbosity level
        mEchoLevel = AssemblySettings["echo_level"].GetInt();
    }


    ///@}
    ///@name Operations
    ///@{

    void SetElementAssemblyFunction(const ElementAssemblyFunctionType& rElementAssemblyFunction)
    {
        mElementAssemblyFunctionIsSet = true;
        mElementAssemblyFunction = rElementAssemblyFunction;
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Element assembly function set." << std::endl;
    }

    void SetConditionAssemblyFunction(const ConditionAssemblyFunctionType& rConditionAssemblyFunction)
    {
        mConditionAssemblyFunctionIsSet = true;
        mConditionAssemblyFunction = rConditionAssemblyFunction;
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Condition assembly function set." << std::endl;
    }

    //TODO: For sure we'll need to modify this one
    void SetConstraintAssemblyFunction(const ConstraintAssemblyFunctionType& rConstraintAssemblyFunction)
    {
        mConstraintAssemblyFunctionIsSet = true;
        mConstraintAssemblyFunction = rConstraintAssemblyFunction;
        KRATOS_INFO_IF("AssemblyHelper", mEchoLevel > 1) << "Constraint assembly function set." << std::endl;
    }

    void AssembleLocalSystem(
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS,
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

    void AssembleLeftHandSide(
        SystemMatrixType& rLHS,
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

    void AssembleRighHandSide(
        SystemVectorType& rRHS,
        TThreadLocalStorage& rTLS)
    {
        // Call the implementation of the function with building type template argument
        if (mBuildType == BuildType::Block) {
            AssembleImplementation<BuildType::Block>(rRHS, rTLS);
        } else if (mBuildType == BuildType::Elimination) {
            AssembleImplementation<BuildType::Elimination>(rRHS, rTLS);
        } else {
            KRATOS_ERROR << "Not implemented build type." << std::endl;
        }
    }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{

    const ModelPart* mpModelPart = nullptr;

    SizeType mEchoLevel;

    SizeType mEquationSystemSize;

    BuildType mBuildType;

    bool mElementAssemblyFunctionIsSet = false;

    bool mConditionAssemblyFunctionIsSet = false;

    bool mConstraintAssemblyFunctionIsSet = false;

    ElementAssemblyFunctionType mElementAssemblyFunction;

    ConditionAssemblyFunctionType mConditionAssemblyFunction;

    ConstraintAssemblyFunctionType mConstraintAssemblyFunction;

    ///@}
    ///@name Private Operations
    ///@{

    template<BuildType TBuildType>
    void AssembleImplementation(
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS,
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
            if (mElementAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    auto it_elem = elems_begin + k;
                    if (it_elem->IsActive()) {
                        // Calculate local LHS and RHS contributions
                        mElementAssemblyFunction(it_elem, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_elem->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
                    }
                }
            }

            // Assemble conditions
            if (mConditionAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    auto it_cond = conds_begin + k;
                    if (it_cond->IsActive()) {
                        // Calculate local LHS and RHS contributions
                        mConditionAssemblyFunction(it_cond, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_cond->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS, rRHS);
                    }
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.FinalizeAssemble();
    }

    template<BuildType TBuildType>
    void AssembleImplementation(
        SystemMatrixType& rLHS,
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
            if (mElementAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    auto it_elem = elems_begin + k;
                    if (it_elem->IsActive()) {
                        // Calculate local LHS contributions
                        mElementAssemblyFunction(it_elem, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_elem->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS);
                    }
                }
            }

            // Assemble conditions
            if (mConditionAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    auto it_cond = conds_begin + k;
                    if (it_cond->IsActive()) {
                        // Calculate local LHS contributions
                        mConditionAssemblyFunction(it_cond, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_cond->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rLHS);
                    }
                }
            }
        }

        // Finalize LHS assembly
        rLHS.FinalizeAssemble();
    }

    template<BuildType TBuildType>
    void AssembleImplementation(
        SystemVectorType& rRHS,
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
            if (mElementAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512) nowait
                for (int k = 0; k < n_elems; ++k) {
                    auto it_elem = elems_begin + k;
                    if (it_elem->IsActive()) {
                        // Calculate local RHS contributions
                        mElementAssemblyFunction(it_elem, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_elem->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rRHS);
                    }
                }
            }

            // Assemble conditions
            if (mConditionAssemblyFunctionIsSet) {
                # pragma omp for schedule(guided, 512)
                for (int k = 0; k < n_conds; ++k) {
                    auto it_cond = conds_begin + k;
                    if (it_cond->IsActive()) {
                        // Calculate local RHS contributions
                        mConditionAssemblyFunction(it_cond, r_process_info, rTLS);

                        // Get the positions in the global system
                        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);
                        it_cond->EquationIdVector(r_loc_eq_ids, r_process_info);

                        // Assemble the local contributions to the global system
                        AssembleLocalContribution<TBuildType>(rTLS, rRHS);
                    }
                }
            }
        }

        // Finalize RHS assembly
        rRHS.BeginAssemble();
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS)
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
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
        }
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        SystemMatrixType& rLHS)
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
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
        }
    }

    template<BuildType TBuildType>
    void AssembleLocalContribution(
        const TThreadLocalStorage& rTLS,
        SystemVectorType& rRHS)
    {
        auto& r_loc_eq_ids = GetThreadLocalStorageEqIds(rTLS);

        if constexpr (TBuildType == BuildType::Block) {
            rRHS.Assemble(GetThreadLocalStorageContainer(rRHS, rTLS), r_loc_eq_ids); // RHS contributions assembly

        } else if (TBuildType == BuildType::Elimination) {
            const auto& r_loc_rhs = GetThreadLocalStorageContainer(rRHS, rTLS);
            const SizeType loc_size = r_loc_rhs.size();
            for (IndexType i_loc = 0; i_loc < loc_size; ++i_loc) {
                IndexType i_glob = r_loc_eq_ids[i_loc];
                if (i_glob < mEquationSystemSize) {// Check if current row DOF is free
                    rRHS.AssembleEntry(r_loc_rhs[i_loc], i_glob); // RHS contribution assembly
                }
            }

        } else {
            static_assert(TBuildType == BuildType::Block || TBuildType == BuildType::Elimination, "Unssupported build type.");
        }
    }

    template<class TContainerType>
    auto& GetThreadLocalStorageContainer(
        const TContainerType& rContainer,
        const TThreadLocalStorage& rTLS)
    {
        if constexpr (std::is_same_v<TContainerType, SystemMatrixType>) {
            return rTLS.LocalLhs; // We can eventually do a get method in the TLS and call it in here
        } else if (std::is_same_v<TContainerType, SystemVectorType>) {
            return rTLS.LocalRhs; // We can eventually do a get method in the TLS and call it in here
        } else {
            static_assert(std::is_same_v<TContainerType, SystemMatrixType> || std::is_same_v<TContainerType, SystemVectorType>, "Unssupported container type.");
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

    ///@}
}; // Class AssemblyHelper

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
