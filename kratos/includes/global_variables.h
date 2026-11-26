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
//                   Riccardo Rossi
//

#pragma once

// System includes
#include <numbers>

// External includes

// Project includes

namespace Kratos::Globals
{
    ///@name Kratos Globals
    ///@{

    // Using the STL to constants variables
	constexpr double Pi = std::numbers::pi;
	constexpr double pi = std::numbers::pi; // NOTE: Consistent lowercase naming is often preferred
    constexpr double e = std::numbers::e;
    constexpr double sqrt2 = std::numbers::sqrt2;
    constexpr double sqrt3 = std::numbers::sqrt3;
    constexpr double log2e = std::numbers::log2e;
    constexpr double log10e = std::numbers::log10e;
    constexpr double phi = std::numbers::phi;

    static constexpr int MaxAllowedThreads = 128; //we assume that no more than MaxAllowedThreads is used in SMP

    ///@}
    ///@name Enums
    ///@{

    /**
     * @brief Enum for Initial and Current configurations
     */
    enum class Configuration
    {
        Initial = 0,
        Current = 1
    };

    /**
     * @brief Enum for location of data
     */
    enum class DataLocation
    {
        NodeHistorical,
        NodeNonHistorical,
        Element,
        Condition,
        ModelPart,
        ProcessInfo,
        NumberOfDataLocations // This needs to be last always, as this is used as a counter for number of DataLocations
    };

    ///@}
}  // namespace Kratos::Globals
