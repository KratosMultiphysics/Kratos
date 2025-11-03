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
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class DofArrayUtilities
 * @ingroup KratosCore
 * @brief Utility class for handling the build of DOFs
 * @details This static class collects some functions to handle the DOFs.
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) DofArrayUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DofArrayUtilities
    KRATOS_CLASS_POINTER_DEFINITION(DofArrayUtilities);

    /// Definition of the index type
    using IndexType = std::size_t;

    /// DOFs array type definition from ModelPart (note that this is PointerVectorSet<DofType>)
    using DofsArrayType = ModelPart::DofsArrayType;

    /// DOFs vector type definition from ModelPart (note that this is std::vector<DofType::Pointer>)s
    using DofsVectorType = ModelPart::DofsVectorType;

    /// Auxilary DOF set type definition
    using AuxiliaryDofsSetType = std::unordered_set<Node::DofType::Pointer, DofPointerHasher>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DofArrayUtilities() = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set the Up Dof Array object
     * This function fills a sorted DOFs array with the DOFs of the provided model part (elements, conditions and constraints)
     * @param rModelPart The model part containing the entities from which the DOFs are to be taken
     * @param rDofArray The DOFs array to be filled
     * @param EchoLevel The echo level (i.e., verbosity level)
     */
    static void SetUpDofArray(
        const ModelPart& rModelPart,
        DofsArrayType& rDofArray,
        const unsigned int EchoLevel = 0);

    /**
     * @brief Set the Up Effective Dof Array and Slave To Master Dofs objects
     * This function fills a map relating each slave DOF with a vector containing its master DOF(s)
     * If there are constraints, the effective (i.e., non-slave) DOF array is set based on this map
     * If there are no constraints, the effective DOF array is the standard DOF array
     * @param rModelPart The model part containing the constraints
     * @param rDofArray The already filled and sorted DOF array
     * @param rEffectiveDofArray The effective DOF array to be filled
     * @param EchoLevel The echo level (i.e., verbosity level)
     */
    static void SetUpEffectiveDofArray(
        const ModelPart& rModelPart,
        const DofsArrayType& rDofArray,
        DofsArrayType& rEffectiveDofArray,
        const unsigned int EchoLevel = 0);

    /**
     * @brief Set the Dof Equation Ids
     * This function sets the equation id for each DOF in the given DOFs array
     * @param rDofArray The already filled and sorted DOF array
     */
    static void SetDofEquationIds(const DofsArrayType& rDofArray);

    /**
     * @brief Set the Effective Dof Equation Ids
     * This function sets the effective equation ids in the effective DOF array
     * If the effective DOF array matches the DOF array (i.e., all DOFs are effective), the effective equation ids are set as the equation ids.
     * If there are non-effective DOFs, the effective equation id of the non-effective DOFs is initialized to the maximum value.
     * This effectively makes possible to distinguish a non-effective DOF by checking its effective equation id
     * @param rDofArray The already filled and sorted DOF array
     * @param rEffectiveDofArray The effective DOFs array in which the effective DOF id are set
     */
    static void SetEffectiveDofEquationIds(
        const DofsArrayType& rDofArray,
        DofsArrayType& rEffectiveDofArray);

    ///@}
}; // Class DofArrayUtilities

///@} addtogroup block

}  // namespace Kratos.
