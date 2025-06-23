//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#pragma once

// System includes

// External includes

// Project includes

namespace Kratos::Globals
{
	constexpr double Pi   = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651L;

    static constexpr int MaxAllowedThreads = 128; //we assume that no more than MaxAllowedThreads is used in SMP

    ///@name Kratos Globals
    ///@{
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
        ProcessInfo
    };

    ///@}
    ///@}
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}

}  // namespace Kratos::Globals
