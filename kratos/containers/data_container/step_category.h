//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

#pragma once

// System includes

// External includes

// Project includes

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name  Enum's
///@{

/**
 * @brief Categories of solution steps tracked by the data container history policies.
 * @details Used by DataHistoryPolicyBase and its derived classes to decide which stored step ring buffer a given operation refers to, and by DataAccessor / DataContainer to request the data of a specific step. Cloning the step data of a DataChunk only advances the internal step index of chunks whose history policy matches the requested category (see HistoricalDataPolicy::CloneStepIndex).
 */
enum class StepCategory {
    AnyStep,            /// Not bound to a specific step category (e.g. non-historical data)
    TimeStep,           /// Time step
    IterationStep,      /// Non-linear iteration step
    SubStep,            /// Sub-step within a time step
    NumberOfCategories  /// Sentinel value holding the number of categories (for array sizing)
};

///@}

///@} addtogroup block

} // namespace Kratos
