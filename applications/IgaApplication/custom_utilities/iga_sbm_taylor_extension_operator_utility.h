//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ricky Aristio
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// Utilities for a classic Taylor-expansion operator for patch coupling
class KRATOS_API(IGA_APPLICATION) IgaSbmTaylorExtensionOperatorUtility
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @param rInterfaceConditions CouplingSbmTaylorInterface6pCondition instances on true boundary
    * @param rMasterSurrogateConditions Master patch's own surrogate boundary conditions
    * @param rSlaveSurrogateConditions Slave patch's own surrogate boundary conditions
    */
    static void PrecomputeAndStoreDualCouplingTaylorData(
        const std::vector<Condition::Pointer>& rInterfaceConditions,
        const std::vector<Condition::Pointer>& rMasterSurrogateConditions,
        const std::vector<Condition::Pointer>& rSlaveSurrogateConditions);

    /**
    * @brief Replaces every condition currently in rCouplingModelPart 
    * @param rCouplingModelPart Model part holding the true-interface conditions to replace
    * @param StartId First id to use for the new conditions
    * @return The newly created Condition instances, in the same order as original conditions
    */
    static std::vector<Condition::Pointer> ReplaceWithInterfaceConditions(
        ModelPart& rCouplingModelPart,
        IndexType StartId);

    ///@}

}; // Class IgaSbmTaylorExtensionOperatorUtility

///@}

}  // namespace Kratos.
