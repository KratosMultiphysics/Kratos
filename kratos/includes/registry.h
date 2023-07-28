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
//                   Ruben Zorrilla
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/registry_item.h"
#include "utilities/parallel_utilities.h"
#include "utilities/string_utilities.h"

namespace Kratos
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @brief Kratos base registry
 * This class is intended to act as global registry
 * Each time the AddItem method is called a pair of name and prototype is called
 */
class KRATOS_API(KRATOS_CORE) Registry final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Registry
    KRATOS_CLASS_POINTER_DEFINITION(Registry);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Registry(){}

    /// Destructor.
    ~Registry(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template<typename TItemType, class... TArgumentsList >
    static RegistryItem& AddItem(
        std::string const& rItemFullName,
        TArgumentsList&&... Arguments)
    {
        KRATOS_TRY

        const std::lock_guard<LockObject> scope_lock(ParallelUtilities::GetGlobalLock());

        auto item_path = StringUtilities::SplitStringByDelimiter(rItemFullName, '.');
        KRATOS_ERROR_IF(item_path.empty()) << "The item full name is empty" << std::endl;

        RegistryItem* p_current_item = &GetRootRegistryItem();

        for(std::size_t i = 0 ; i < item_path.size() - 1 ; i++){
            auto& r_item_name = item_path[i];
            if(p_current_item->HasItem(r_item_name)){
                p_current_item = &p_current_item->GetItem(r_item_name);
            }
            else{
                p_current_item = &p_current_item->AddItem<RegistryItem>(r_item_name);
            }
        }

        // I am doing the last one out of the loop to create it with the given type and argument
        auto& r_item_name = item_path.back();
        if(p_current_item->HasItem(r_item_name)){
            KRATOS_ERROR << "The item \"" << rItemFullName << "\" is already registered." << std::endl;
        }
        else{
            p_current_item = &p_current_item->AddItem<TItemType>(r_item_name, std::forward<TArgumentsList>(Arguments)...);
        }

        return *p_current_item;

        KRATOS_CATCH("")
    }

    template<typename... Types>
    static std::string RegistryTemplateToString(Types&&... args) {
        std::string f_name = (... += ("," + std::to_string(args)));
        f_name.erase(0,1);
        return f_name;
    }

    ///@}
    ///@name Access
    ///@{

    static auto begin()
    {
        return mspRootRegistryItem->begin();
    }

    static auto cbegin()
    {
        return mspRootRegistryItem->cbegin();
    }

    static auto end()
    {
        return mspRootRegistryItem->end();
    }

    static auto const cend()
    {
        return mspRootRegistryItem->cend();
    }

    static RegistryItem& GetItem(std::string const& rItemFullName);

    template<typename TDataType>
    static TDataType const& GetValue(std::string const& rItemFullName)
    {
        return GetItem(rItemFullName).GetValue<TDataType>();
    }

    static void RemoveItem(std::string const& ItemName);

    ///@}
    ///@name Inquiry
    ///@{

    static std::size_t size();

    static bool HasItem(std::string const& rItemFullName);

    static bool HasValue(std::string const& rItemFullName);

    static bool HasItems(std::string const& rItemFullName);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    std::string ToJson(std::string const& Indentation) const;

    ///@}
    ///@name Friends
    ///@{


    ///@}
private:
    ///@name Static Member Variables
    ///@{

        static RegistryItem* mspRootRegistryItem;

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

        static RegistryItem& GetRootRegistryItem();

        static std::vector<std::string> SplitFullName(std::string const& FullName);

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    Registry& operator=(Registry const& rOther);

    /// Copy constructor.
    Registry(Registry const& rOther);

    ///@}
}; // Class Registry

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    Registry& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const Registry& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
///@} addtogroup block

}  // namespace Kratos.
