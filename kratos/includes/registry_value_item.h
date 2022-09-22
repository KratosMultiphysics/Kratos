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
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/registry_item.h"

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

/// The registry value item is a registry item which stores a copy of a value.
/** This class is derived from RegistryItem and the difference is that it creates a
 *  copy of the value in initialization and delete that copy at destruction.
 *  This make it suitable for storing values that we want to have one copy of them
 *  in the registry but not for the values that are stored in other places (like static variables)
 *  Value type should be copy constructable.
*/
template<typename TValueType>
class RegistryValueItem : public RegistryItem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RegistryValueItem
    KRATOS_CLASS_POINTER_DEFINITION(RegistryValueItem);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    RegistryValueItem() = delete;

    /// Constructor with the name
    RegistryValueItem(std::string Name) : RegistryItem(Name){}

    /// Constructor with the name and value
    template<class... TArgumentsList >
    RegistryValueItem(
        std::string const& Name,
        TArgumentsList&&... Arguments)
        : RegistryItem(Name, *(new TValueType(std::forward<TArgumentsList>(Arguments)...)))
    {}

    /// Destructor.
    virtual ~RegistryValueItem()
    {
        delete static_cast<TValueType*>(const_cast<void*>(mpValue));
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return Name() + " RegistryValueItem ";
    }

    std::string ToJson(std::string const& Indentation = "") const override
    {
        std::stringstream buffer;

        if(HasValue()){
            buffer << Indentation << "\"" << mName << "\" : " << ValueToJson() << std::endl;
        }
        else{
            buffer << Indentation << "\"" << mName << "\" : \"\""  << std::endl;

        }

        return buffer.str();
    }

    std::string ValueToJson() const
    {
        std::stringstream buffer;
        buffer << "\"" << *static_cast<const TValueType*>(mpValue) << "\"";

        return buffer.str();
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    RegistryValueItem& operator=(RegistryValueItem const& rOther);

    /// Copy constructor.
    RegistryValueItem(RegistryValueItem const& rOther);


    ///@}

}; // Class RegistryValueItem

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.
