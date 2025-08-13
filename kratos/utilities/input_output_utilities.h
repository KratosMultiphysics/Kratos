//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore

///@name Kratos Classes
///@{
/**
 * @class InputOutputUtilities
 * @ingroup KratosCore
 * @brief Utilities for common input/output operations
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) InputOutputUtilities
{
public:
    /**
     * @brief Checks if the given entity can be skipped during output
     * @param rEntity The entity to be checked
     * @param rWarningLabel The warning label to be used
     * @return true if the entity can be skipped, false otherwise
     * @tparam TEntityType the type of the entity
     */
    template<typename TEntityType>
    static bool SkippableEntity(
        const TEntityType& rEntity,
        const std::string& rWarningLabel
        );

    /**
     * @brief This method generates a string label for the step with a given step and the number of spaces considered.
     * @details For example Step = 3 and Spaces = 6, would be "000003".
     * @param Step The numerical value of the step.
     * @param Spaces The total desired width of the output string, padded with leading zeros.
     * @return A zero-padded string representation of the step.
     */
    static std::string GenerateStepLabel(
        const std::size_t Step,
        const std::size_t Spaces
        );
};

///@}

} // namespace Kratos
