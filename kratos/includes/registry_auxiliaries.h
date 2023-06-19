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
//                   Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/registry.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

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
///@name Kratos Classes
///@{

/**
 * @brief Kratos registry auxiliaries
 * This static class collect all the auxiliary functions to be used to register c++ items
 */
class RegistryAuxiliaries final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Registry
    KRATOS_CLASS_POINTER_DEFINITION(RegistryAuxiliaries);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RegistryAuxiliaries() = delete;

    /// Copy constructor.
    RegistryAuxiliaries(RegistryAuxiliaries const& rOther) = delete;

    /// Destructor.
    ~RegistryAuxiliaries() = delete;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<typename TPrototypeType>
    static void RegisterProcessWithPrototype(
        const std::string ModuleName,
        const std::string ProcessName,
        TPrototypeType rProcessPrototype)
    {
        const std::string all_path = std::string("Processes.All.") + ProcessName;
        RegisterPrototype(all_path, rProcessPrototype);
        const std::string module_path = std::string("Processes.") + ModuleName + std::string(".") + ProcessName;
        RegisterPrototype(module_path, rProcessPrototype);
    }

    template<typename TPrototypeType>
    static void RegisterOperationWithPrototype(
        const std::string ModuleName,
        const std::string OperationName,
        TPrototypeType rOperationPrototype)
    {
        const std::string all_path = std::string("Operations.All.") + OperationName;
        RegisterPrototype(all_path, rOperationPrototype);
        const std::string module_path = std::string("Operations.") + ModuleName + std::string(".") + OperationName;
        RegisterPrototype(module_path, rOperationPrototype);
    }

    template<typename TPrototypeType>
    static void RegisterControllerWithPrototype(
        const std::string ModuleName,
        const std::string ControllerName,
        TPrototypeType rControllerPrototype)
    {
        const std::string all_path = std::string("Controllers.All.") + ControllerName;
        RegisterPrototype(all_path, rControllerPrototype);
        const std::string module_path = std::string("Controllers.") + ModuleName + std::string(".") + ControllerName;
        RegisterPrototype(module_path, rControllerPrototype);
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    template<typename TPrototypeType>
    static void RegisterPrototype(
        const std::string RegistryEntryName,
        TPrototypeType rPrototype)
    {
        if (!Registry::HasItem(RegistryEntryName)) {
            auto& r_item = Registry::AddItem<RegistryItem>(RegistryEntryName);
            r_item.AddItem<TPrototypeType>("Prototype", rPrototype);
        } else {
            KRATOS_ERROR << "'" << RegistryEntryName << "' is already registered." << std::endl;
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class RegistryAuxiliaries

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
