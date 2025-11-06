
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
#include "linear_operator.h"
#include "includes/model_part.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Auxiliary container to store the linear system
 * This auxiliary container is intended to store all the arrays requires for the linear system setup
 * @tparam TSparseMatrixType The sparse matrix type
 * @tparam TSystemVectorType The system vector type
 * TODO: to be changed to LinearSystem and have the matrices and vectors stored by their name in a dictionary, so we can reuse them if needed in the strategy
 */
template <class TSparseMatrixType, class TSystemVectorType>
struct LinearSystemContainer
{

    /// Pointer definition of LinearSystemContainer
    KRATOS_CLASS_POINTER_DEFINITION(LinearSystemContainer);

    typename TSparseMatrixType::Pointer pLhs = nullptr; // Pointer to the LHS matrix

    typename TSystemVectorType::Pointer pRhs = nullptr; // Pointer to the RHS vector

    typename TSystemVectorType::Pointer pDx = nullptr; // Pointer to the solution increment vector

    typename TSparseMatrixType::Pointer pEffectiveLhs = nullptr; // Pointer to the effective LHS matrix (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveRhs = nullptr; // Pointer to the effective RHS vector (i.e., after applying system constraints)

    typename TSystemVectorType::Pointer pEffectiveDx = nullptr; // Pointer to the effective solution increment vector (i.e., after applying system constraints)

    typename TSparseMatrixType::Pointer pEffectiveT = nullptr; // Linear system constraints total relation matrix

    typename TSparseMatrixType::Pointer pConstraintsT = nullptr; // Master-slave constraints relation matrix

    typename TSystemVectorType::Pointer pConstraintsQ = nullptr; // Master-slave constraints constant vector

    typename ModelPart::DofsArrayType::Pointer pDofSet = Kratos::make_shared<ModelPart::DofsArrayType>(); // The PVS containing the DOFs of the system

    typename ModelPart::DofsArrayType::Pointer pEffectiveDofSet = Kratos::make_shared<ModelPart::DofsArrayType>(); /// The PVS containing the effective DOFs of the system

    void Clear()
    {
        if (pLhs) {
            pLhs->Clear();
        }
        if (pRhs) {
            pRhs->Clear();
        }
        if (pDx) {
            pDx->Clear();
        }
        if (pEffectiveLhs) {
            pEffectiveLhs->Clear();
        }
        if (pEffectiveRhs) {
            pEffectiveRhs->Clear();
        }
        if (pEffectiveDx) {
            pEffectiveDx->Clear();
        }
        if (pEffectiveT) {
            pEffectiveT->Clear();
        }
        if (pConstraintsT) {
            pConstraintsT->Clear();
        }
        if (pConstraintsQ) {
            pConstraintsQ->Clear();
        }
        if (pDofSet) {
            pDofSet->clear();
        }
        if (pEffectiveDofSet) {
            pEffectiveDofSet->clear();
        }
    }

    bool RequiresEffectiveDofSet() const
    {
        return pDofSet == pEffectiveDofSet;
    }

}; // Class LinearSolverContainer

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.