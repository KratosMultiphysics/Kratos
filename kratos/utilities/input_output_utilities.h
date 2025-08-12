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
};

///@}

} // namespace Kratos
