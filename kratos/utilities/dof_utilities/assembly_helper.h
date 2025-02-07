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
        Parameters AssemblySettings = Parameters(R"({})"))
    : mpModelPart(&rModelPart)
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

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{

    SizeType mEchoLevel;

    BuildType mBuildType;

    const ModelPart* mpModelPart = nullptr;

    bool mElementAssemblyFunctionIsSet = false;

    bool mConditionAssemblyFunctionIsSet = false;

    bool mConstraintAssemblyFunctionIsSet = false;

    ElementAssemblyFunctionType mElementAssemblyFunction;

    ConditionAssemblyFunctionType mConditionAssemblyFunction;

    ConstraintAssemblyFunctionType mConstraintAssemblyFunction;

    ///@}
    ///@name Private Operations
    ///@{

    //TODO: Do the corresponding ones for RHS and LHS only
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
                        it_elem->EquationIdVector(rTLS.LocalEqIds, r_process_info);

                        // Assemble the local contributions to the global system
                        if constexpr (TBuildType == BuildType::Block) {
                            rRHS.Assemble(rTLS.LocalRhs, rTLS.LocalEqIds);
                            rLHS.Assemble(rTLS.LocalLhs, rTLS.LocalEqIds);
                        } else if (TBuildType == BuildType::Elimination) {
                            KRATOS_ERROR << "To be implemented." << std::endl;
                        } else {
                            KRATOS_ERROR << "Not implemented build type." << std::endl;
                        }
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
                        it_cond->EquationIdVector(rTLS.LocalEqIds, r_process_info);

                        // Assemble the local contributions to the global system
                        if constexpr (TBuildType == BuildType::Block) {
                            rRHS.Assemble(rTLS.LocalRhs, rTLS.LocalEqIds);
                            rLHS.Assemble(rTLS.LocalLhs, rTLS.LocalEqIds);
                        } else if (TBuildType == BuildType::Elimination) {
                            KRATOS_ERROR << "To be implemented." << std::endl;
                        } else {
                            KRATOS_ERROR << "Not implemented build type." << std::endl;
                        }
                    }
                }
            }
        }

        // Finalize RHS and LHS assembly
        rRHS.BeginAssemble();
        rLHS.FinalizeAssemble();
    }

    ///@}
}; // Class AssemblyHelper

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
