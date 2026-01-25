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
#include "includes/model_part.h"
#include "future/containers/linear_system.h"
#include "future/linear_operators/linear_operator.h"

namespace Kratos::Future
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Auxiliary container to store the linear system
 * This auxiliary container is intended to store all the arrays requires for the linear system setup
 * @tparam TLinearAlgebra The linear algebra type
 * TODO: to be changed to LinearSystem and have the matrices and vectors stored by their name in a dictionary, so we can reuse them if needed in the strategy
 */
template <class TLinearAlgebra>
struct ImplicitStrategyDataContainer
{

    // ImplicitStrategyDataContainer()
    // {
    //     mpLinearSystem = LinearSystem<TLinearAlgebra>(pLhs, pDx, pRhs);
    //     mpEffectiveLinearSystem = LinearSystem<TLinearAlgebra>(pEffectiveLhs, pEffectiveDx, pEffectiveRhs);
    // }

    /// Pointer definition of ImplicitStrategyDataContainer
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitStrategyDataContainer);

    using MatrixType = typename TLinearAlgebra::MatrixType;

    using VectorType = typename TLinearAlgebra::VectorType;

    using LinearSystemType = LinearSystem<TLinearAlgebra>;

    // //TODO: All these must dissapear from here
    // typename MatrixType::Pointer pLhs = nullptr; // Pointer to the LHS matrix

    // //TODO: All these must dissapear from here
    // typename VectorType::Pointer pRhs = nullptr; // Pointer to the RHS vector

    // //TODO: All these must dissapear from here
    // typename VectorType::Pointer pDx = nullptr; // Pointer to the solution increment vector

    // //TODO: All these must dissapear from here
    // typename MatrixType::Pointer pEffectiveLhs = nullptr; // Pointer to the effective LHS matrix (i.e., after applying system constraints)

    // //TODO: All these must dissapear from here
    // typename VectorType::Pointer pEffectiveRhs = nullptr; // Pointer to the effective RHS vector (i.e., after applying system constraints)

    // //TODO: All these must dissapear from here
    // typename VectorType::Pointer pEffectiveDx = nullptr; // Pointer to the effective solution increment vector (i.e., after applying system constraints)

    typename MatrixType::Pointer pEffectiveT = nullptr; // Linear system constraints total relation matrix (combination of MPCs and eventual Dirichlet BCs)

    typename MatrixType::Pointer pConstraintsT = nullptr; // Master-slave constraints relation matrix

    typename VectorType::Pointer pConstraintsQ = nullptr; // Master-slave constraints constant vector

    typename ModelPart::DofsArrayType::Pointer pDofSet = Kratos::make_shared<ModelPart::DofsArrayType>(); // The PVS containing the DOFs of the system

    typename ModelPart::DofsArrayType::Pointer pEffectiveDofSet = Kratos::make_shared<ModelPart::DofsArrayType>(); /// The PVS containing the effective DOFs of the system

    typename LinearSystemType::Pointer mpLinearSystem = nullptr;

    typename LinearSystemType::Pointer mpEffectiveLinearSystem = nullptr;

    /**
     * @brief Returns the pointer to the linear system
     * @return const LinearSystemType::Pointer The pointer to the linear system
     */
    typename LinearSystemType::Pointer pGetLinearSystem()
    {
        return mpLinearSystem;
    }

    /**
     * @brief Returns the pointer to the linear system
     * @return const LinearSystemType::Pointer The pointer to the linear system
     */
    const typename LinearSystemType::Pointer pGetLinearSystem() const
    {
        return mpLinearSystem;
    }

    /**
     * @brief Sets the pointer to the linear system
     * @param pLinearSystem The pointer to the linear system
     */
    void pSetLinearSystem(typename LinearSystemType::Pointer pLinearSystem)
    {
        mpLinearSystem = pLinearSystem;
    }

    /**
     * @brief Returns the pointer to the effective linear system
     * @return const LinearSystemType::Pointer The pointer to the effective linear system
     */
    typename LinearSystemType::Pointer pGetEffectiveLinearSystem()
    {
        return mpEffectiveLinearSystem;
    }

    /**
     * @brief Returns the pointer to the effective linear system
     * @return const LinearSystemType::Pointer The pointer to the effective linear system
     */
    const typename LinearSystemType::Pointer pGetEffectiveLinearSystem() const
    {
        return mpEffectiveLinearSystem;
    }

    /**
     * @brief Sets the pointer to the effective linear system
     * @param pEffectiveLinearSystem The pointer to the effective linear system
     */
    void pSetEffectiveLinearSystem(typename LinearSystemType::Pointer pEffectiveLinearSystem)
    {
        mpEffectiveLinearSystem = pEffectiveLinearSystem;
    }

    void Clear()
    {
        // if (pLhs) {
        //     pLhs->Clear();
        // }
        // if (pRhs) {
        //     pRhs->Clear();
        // }
        // if (pDx) {
        //     pDx->Clear();
        // }
        // if (pEffectiveLhs) {
        //     pEffectiveLhs->Clear();
        // }
        // if (pEffectiveRhs) {
        //     pEffectiveRhs->Clear();
        // }
        // if (pEffectiveDx) {
        //     pEffectiveDx->Clear();
        // }
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

}; // Class ImplicitStrategyDataContainer

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.