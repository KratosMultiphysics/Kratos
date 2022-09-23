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
#include <unordered_map>

// External includes


// Project includes
#include "includes/define.h"

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

/// The registry item to be stored by Registry class. It is the base class for some more specific ones.
/** RegistryItem has a tree node structure and stores its name, an optional
 *  value, and an unorder_set of its sub data.
 *  This structure let us to have registry of the elements and then different
 *  registries for each elements inside it.
 *  Please note that RegistryItem stores a pointer to the value.
 *  To have a copy of the value you may use the derived RegistryValueItem
 *  which crates a copy in construction and delete it in its destructor
 *  to make the memory management easier.
*/
class RegistryItem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RegistryItem
    KRATOS_CLASS_POINTER_DEFINITION(RegistryItem);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    RegistryItem() = delete;

    /// Constructor with the name
    RegistryItem(std::string Name) : mName(Name), mpValue(nullptr){}

    /// Constructor with the name and value
    template<typename TValueType>
    RegistryItem(std::string Name, TValueType const& Value) : mName(Name), mpValue(&Value){}

    /// Destructor.
    virtual ~RegistryItem(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Items
    ///@{

    template< typename TItemType, class... TArgumentsList >
    RegistryItem& AddItem(
        std::string const& ItemName,
        TArgumentsList&&... Arguments)
    {
        KRATOS_ERROR_IF(this->HasItem(ItemName)) << "The RegistryItem " << this->Name() << " already has an item with name " << ItemName << std::endl;
        auto insert_result = mSubRegistryItem.emplace(std::make_pair(ItemName, Kratos::make_unique<TItemType>(ItemName, std::forward<TArgumentsList>(Arguments)...)));
        KRATOS_ERROR_IF_NOT(insert_result.second) << "Error in inserting " << ItemName << " in registry item with name " << this->Name() << std::endl;
        return *insert_result.first->second;
    }



    ///@}
    ///@name Access
    ///@{

    const std::string& Name() const
    {
        return mName;
    }

    RegistryItem const& GetItem(std::string const& rItemName) const;

    RegistryItem& GetItem(std::string const& rItemName);

    void RemoveItem(std::string const& rItemName);

    ///@}
    ///@name Inquiry
    ///@{

    bool HasValue() const
    {
        return (mpValue != nullptr);
    }

    bool HasItems() const
    {
        return (!mSubRegistryItem.empty());
    }

    bool HasItem(std::string const& rItemName) const
    {
        return (mSubRegistryItem.find(rItemName) != mSubRegistryItem.end());
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

    virtual std::string ToJson(std::string const& rIndetation = "") const;

    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    std::string mName;
    const void* mpValue;
    std::unordered_map<std::string, Kratos::unique_ptr<RegistryItem>> mSubRegistryItem;

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
    RegistryItem& operator=(RegistryItem const& rOther);

    /// Copy constructor.
    RegistryItem(RegistryItem const& rOther);

    ///@}
}; // Class RegistryItem

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    RegistryItem& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const RegistryItem& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.
