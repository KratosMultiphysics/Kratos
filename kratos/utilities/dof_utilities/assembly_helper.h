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
template<class TThreadLocalStorage, class TDataType = double, class TIndexType=std::size_t, class TMatrixType=CsrMatrix<TDataType, TIndexType>, class TVectorType=SystemVector<TDataType, TIndexType>>
class AssemblyHelper
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssemblyHelper
    KRATOS_CLASS_POINTER_DEFINITION(AssemblyHelper);

    /// Index type definition
    using IndexType = TIndexType;

    /// Size type definition
    using SizeType = std::size_t;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Vector type definition
    using SystemVectorType = TVectorType;

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    AssemblyHelper() = delete;

    /// Constructor with model part
    AssemblyHelper(const ModelPart& rModelPart)
    : mpModelPart(&rModelPart)
    {}

    ///@}
    ///@name Operations
    ///@{

    void SetElementAssemblyFunction(const ElementAssemblyFunctionType& rElementAssemblyFunction)
    {
        mElementAssemblyFunctionIsSet = true;
        mElementAssemblyFunction = rElementAssemblyFunction;
    }

    void SetConditionAssemblyFunction(const ConditionAssemblyFunctionType& rConditionAssemblyFunction)
    {
        mConditionAssemblyFunctionIsSet = true;
        mConditionAssemblyFunction = rConditionAssemblyFunction;
    }

    void SetConstraintAssemblyFunction(const ConstraintAssemblyFunctionType& rConstraintAssemblyFunction)
    {
        mConstraintAssemblyFunctionIsSet = true;
        mConstraintAssemblyFunction = rConstraintAssemblyFunction;
    }

    void AssembleLocalSystem(
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
                        it_elem->EquationIdVector(rTLS.LocEqIds, r_process_info);

                        // Assemble the local contributions to the global system
                        rRHS.Assemble(rTLS.LocRhs, rTLS.LocEqIds);
                        rLHS.Assemble(rTLS.LocLhs, rTLS.LocEqIds);
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
                        it_cond->EquationIdVector(rTLS.LocEqIds, r_process_info);

                        // Assemble the local contributions to the global system
                        rRHS.Assemble(rTLS.LocRhs, rTLS.LocEqIds);
                        rLHS.Assemble(rTLS.LocLhs, rTLS.LocEqIds);
                    }
                }
            }
        }
    }

    // template<class TAssemblyImplementation, class TThreadLocalStorage>
    // void AssembleLocalSystemElements(
    //     const ModelPart& rModelPart,
    //     const TAssemblyImplementation& rAssemblyImplementation,
    //     SystemMatrixType& rLHS,
    //     SystemVectorType& rRHS,
    //     TThreadLocalStorage& rTLS)
    // {
    //     Assemble(rModelPart.Elements(), rModelPart.GetProcessInfo(), rAssemblyImplementation, rLHS, rRHS, rTLS);
    // }

    // template<class TAssemblyImplementation, class TThreadLocalStorage>
    // void AssembleLocalSystemConditions(
    //     const ModelPart& rModelPart,
    //     const TAssemblyImplementation& rAssemblyImplementation,
    //     SystemMatrixType& rLHS,
    //     SystemVectorType& rRHS,
    //     TThreadLocalStorage& rTLS)
    // {
    //     Assemble(rModelPart.Conditions(), rModelPart.GetProcessInfo(), rAssemblyImplementation, rLHS, rRHS, rTLS);
    // }

    ///@}
    ///@name Input and output
    ///@{


    ///@}
private:
    ///@name Member Variables
    ///@{

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

    // template<class TEntityContainer, class TAssemblyImplementation, class TThreadLocalStorage>
    // static void Assemble(
    //     const TEntityContainer& rEntities,
    //     const ProcessInfo& rProcessInfo,
    //     const TAssemblyImplementation& rAssemblyImplementation,
    //     SystemMatrixType& rLHS,
    //     SystemVectorType& rRHS,
    //     TThreadLocalStorage& rTLS)
    // {
    //     // Getting entities container data
    //     auto ent_begin = rEntities.begin();
    //     const SizeType n_ent = rEntities.size();

    //     // Assemble entities
    //     #pragma omp parallel firstprivate(n_ent, rTLS, rProcessInfo)
    //     {
    //         // Assemble elements
    //         # pragma omp for schedule(guided, 512) nowait
    //         for (int k = 0; k < n_ent; ++k) {
    //             auto it_ent = ent_begin + k;
    //             if (it_ent->IsActive()) {
    //                 // Calculate local LHS and RHS contributions
    //                 rAssemblyImplementation(it_ent, rProcessInfo, rTLS);

    //                 // Get the positions in the global system
    //                 it_ent->EquationIdVector(rTLS.LocEqIds, rProcessInfo);

    //                 // Assemble the local contributions to the global system
    //                 rRHS.Assemble(rTLS.LocRhs, rTLS.LocEqIds);
    //                 rLHS.Assemble(rTLS.LocLhs, rTLS.LocEqIds);
    //             }
    //         }
    //     }
    // }

    ///@}
}; // Class AssemblyHelper

///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.
