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
 */
template <class TLinearAlgebra>
class ImplicitStrategyData
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImplicitStrategyData
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitStrategyData);

    /// Type alias for the matrix type
    using MatrixType = typename TLinearAlgebra::MatrixType;

    /// Type alias for the vector type
    using VectorType = typename TLinearAlgebra::VectorType;

    /// Type alias for the linear system type
    using LinearSystemType = LinearSystem<TLinearAlgebra>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ImplicitStrategyData() = default;

    /// Copy constructor
    ImplicitStrategyData(const ImplicitStrategyData& other) = delete;

    /// Move constructor
    ImplicitStrategyData(ImplicitStrategyData&& other) = delete;

    /// Destructor
    ~ImplicitStrategyData() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy assignment operator
    ImplicitStrategyData& operator=(const ImplicitStrategyData& other) = delete;

    /// Move assignment operator
    ImplicitStrategyData& operator=(ImplicitStrategyData&& other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clear all the pointers
     * @details This function is intended to be called to leave the object in a clean state
     * Note that the clear function does not deallocate the memory of the arrays as this
     * container is not the owner of these. Conversely, it clears the DOF sets.
     */
    void Clear()
    {
        // Clear the linear system and constraints pointers
        mpEffectiveT = nullptr;
        mpConstraintsT = nullptr;
        mpConstraintsQ = nullptr;
        mpLinearSystem = nullptr;
        mpEffectiveLinearSystem = nullptr;

        // Clear the DOF set PVS
        if (mpDofSet) {
            mpDofSet->clear();
        }
        if (mpEffectiveDofSet) {
            mpEffectiveDofSet->clear();
        }
    }

    ///@}
    ///@name Access
    ///@{

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

        /**
     * @brief Returns the pointer to the effective relation matrix
     * @return const MatrixType::Pointer The pointer to the effective relation matrix
     */
    typename MatrixType::Pointer pGetEffectiveT()
    {
        return mpEffectiveT;
    }

    /**
     * @brief Returns the pointer to the effective relation matrix
     * @return const MatrixType::Pointer The pointer to the effective relation matrix
     */
    const typename MatrixType::Pointer pGetEffectiveT() const
    {
        return mpEffectiveT;
    }

    /**
     * @brief Sets the pointer to the effective relation matrix
     * @param pEffectiveT The pointer to the effective relation matrix
     */
    void pSetEffectiveT(typename MatrixType::Pointer pEffectiveT)
    {
        mpEffectiveT = pEffectiveT;
    }

    /**
     * @brief Returns the pointer to the constraints relation matrix
     * @return const MatrixType::Pointer The pointer to the constraints relation matrix
     */
    typename MatrixType::Pointer pGetConstraintsT()
    {
        return mpConstraintsT;
    }

    /**
     * @brief Returns the pointer to the constraints relation matrix
     * @return const MatrixType::Pointer The pointer to the constraints relation matrix
     */
    const typename MatrixType::Pointer pGetConstraintsT() const
    {
        return mpConstraintsT;
    }

    /**
     * @brief Sets the pointer to the constraints relation matrix
     * @param pConstraintsT The pointer to the constraints relation matrix
     */
    void pSetConstraintsT(typename MatrixType::Pointer pConstraintsT)
    {
        mpConstraintsT = pConstraintsT;
    }

    /**
     * @brief Returns the pointer to the constraints constant vector
     * @return const VectorType::Pointer The pointer to the constraints constant vector
     */
    typename VectorType::Pointer pGetConstraintsQ()
    {
        return mpConstraintsQ;
    }

    /**
     * @brief Returns the pointer to the constraints constant vector
     * @return const VectorType::Pointer The pointer to the constraints constant vector
     */
    const typename VectorType::Pointer pGetConstraintsQ() const
    {
        return mpConstraintsQ;
    }

    /**
     * @brief Sets the pointer to the constraints constant vector
     * @param pConstraintsQ The pointer to the constraints constant vector
     */
    void pSetConstraintsQ(typename VectorType::Pointer pConstraintsQ)
    {
        mpConstraintsQ = pConstraintsQ;
    }

    /**
     * @brief Returns the pointer to the DOF set
     * @return const ModelPart::DofsArrayType::Pointer The pointer to the DOF set
     */
    typename ModelPart::DofsArrayType::Pointer pGetDofSet()
    {
        return mpDofSet;
    }

    /**
     * @brief Returns the pointer to the DOF set
     * @return const ModelPart::DofsArrayType::Pointer The pointer to the DOF set
     */
    const typename ModelPart::DofsArrayType::Pointer pGetDofSet() const
    {
        return mpDofSet;
    }

    /**
     * @brief Sets the pointer to the DOF set
     * @param pDofSet The pointer to the DOF set
     */
    void pSetDofSet(typename ModelPart::DofsArrayType::Pointer pDofSet)
    {
        mpDofSet = pDofSet;
    }

    /**
     * @brief Returns the pointer to the effective DOF set
     * @return const ModelPart::DofsArrayType::Pointer The pointer to the effective DOF set
     */
    typename ModelPart::DofsArrayType::Pointer pGetEffectiveDofSet()
    {
        return mpEffectiveDofSet;
    }

    /**
     * @brief Returns the pointer to the effective DOF set
     * @return const ModelPart::DofsArrayType::Pointer The pointer to the effective DOF set
     */
    const typename ModelPart::DofsArrayType::Pointer pGetEffectiveDofSet() const
    {
        return mpEffectiveDofSet;
    }

    /**
     * @brief Sets the pointer to the effective DOF set
     * @param pEffectiveDofSet The pointer to the effective DOF set
     */
    void pSetEffectiveDofSet(typename ModelPart::DofsArrayType::Pointer pEffectiveDofSet)
    {
        mpEffectiveDofSet = pEffectiveDofSet;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if the DOF set and the effective DOF set are the same
     * @return bool True if the DOF set and the effective DOF set are the same, false otherwise
     */
    bool RequiresEffectiveDofSet() const
    {
        return mpDofSet == mpEffectiveDofSet;
    }

    ///@}
private:

    ///@name Member Variables
    ///@{

    /// Linear system constraints total relation matrix (combination of MPCs and eventual Dirichlet BCs)
    typename MatrixType::Pointer mpEffectiveT = nullptr;

    /// Master-slave constraints relation matrix
    typename MatrixType::Pointer mpConstraintsT = nullptr;

    /// Master-slave constraints constant vector
    typename VectorType::Pointer mpConstraintsQ = nullptr;

    /// Linear system (i.e., linear system before applying constraints)
    typename LinearSystemType::Pointer mpLinearSystem = nullptr;

    /// Effective linear system (i.e., linear system after applying constraints)
    typename LinearSystemType::Pointer mpEffectiveLinearSystem = nullptr;

    /// The PVS containing the DOFs of the system
    typename ModelPart::DofsArrayType::Pointer mpDofSet = Kratos::make_shared<ModelPart::DofsArrayType>();

    /// The PVS containing the effective DOFs of the system
    typename ModelPart::DofsArrayType::Pointer mpEffectiveDofSet = Kratos::make_shared<ModelPart::DofsArrayType>();


    ///@}
}; // Class ImplicitStrategyData

///@}
///@name Input and output
///@{

///@}
///@} addtogroup block

} // namespace Kratos.