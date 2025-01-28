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
 * @class Registry
 * @ingroup KratosCore
 * @brief Kratos base registry class
 * This class serves as a global registry to store hierarchical data or object prototypes.
 * Items in the registry are organized as a tree structure, where each node can have sub-nodes (items) or values.
 * The `Registry` class supports adding, removing, and retrieving items, and is optimized for thread safety.
 * @details
 * - The registry is a singleton, ensuring a unique global instance.
 * - Each registry item is represented by the `RegistryItem` class.
 * - Thread safety is guaranteed via locks.
 * @author Pooyan Dadvand
 *         Ruben Zorrilla
 * @see RegistryItem
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

    /**
     * @brief Default constructor
     */
    Registry() {}

    /**
     * @brief Destructor
     */
    ~Registry() {}

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Adds a new item to the registry
     * @tparam TItemType The type of the item to add
     * @tparam TArgumentsList Variadic template arguments to construct the item
     * @param rItemFullName Full hierarchical path of the item (e.g., "Physics.Gravity")
     * @param Arguments Arguments used to construct the new item
     * @return Reference to the newly created registry item
     * @throws An error if the item already exists
     * @note If the full path does not exist, intermediate nodes are created as `RegistryItem` instances.
     */
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

        for(std::size_t i = 0 ; i < item_path.size() - 1 ; i++) {
            auto& r_item_name = item_path[i];
            if(p_current_item->HasItem(r_item_name)){
                p_current_item = &p_current_item->GetItem(r_item_name);
            } else {
                p_current_item = &p_current_item->AddItem<RegistryItem>(r_item_name);
            }
        }

        // I am doing the last one out of the loop to create it with the given type and argument
        auto& r_item_name = item_path.back();
        if (p_current_item->HasItem(r_item_name)) {
            KRATOS_ERROR << "The item \"" << rItemFullName << "\" is already registered." << std::endl;
        } else {
            p_current_item = &p_current_item->AddItem<TItemType>(r_item_name, std::forward<TArgumentsList>(Arguments)...);
        }

        return *p_current_item;

        KRATOS_CATCH("")
    }

    /**
     * @brief Converts template arguments to a string representation
     * @tparam Types Variadic template types
     * @param args Arguments to convert to string
     * @return A comma-separated string representation of the arguments
     */
    template<typename... Types>
    static std::string RegistryTemplateToString(Types&&... args) {
        std::string f_name = (... += ("," + std::to_string(args)));
        f_name.erase(0,1);
        return f_name;
    }

    ///@}
    ///@name Iterators
    ///@{

    /**
     * @brief The begin iterator of the registry
     * @return An iterator to the first item in the registry
     */
    static auto begin()
    {
        return mspRootRegistryItem->begin();
    }

    /**
     * @brief The begin iterator of the registry (const version)
     * @return A const iterator to the first item in the registry
     */
    static auto cbegin()
    {
        return mspRootRegistryItem->cbegin();
    }

    /**
     * @brief The end iterator of the registry
     * @return An iterator to the end of the registry
     */
    static auto end()
    {
        return mspRootRegistryItem->end();
    }

    /**
     * @brief The end iterator of the registry (const version)
     * @return A const iterator to the end of the registry
     */
    static auto cend()
    {
        return mspRootRegistryItem->cend();
    }

    ///@}
    ///@name Access
    ///@{


    /**
     * @brief Retrieves a registry item by its full hierarchical path
     * @param rItemFullName Full hierarchical path of the item
     * @return Reference to the registry item
     */
    static RegistryItem& GetItem(std::string const& rItemFullName);

    /**
     * @brief Retrieves the value stored in a registry item
     * @tparam TDataType Type of the value
     * @param rItemFullName Full hierarchical path of the item
     * @return Const reference to the value stored in the item
     * @throws An error if the item does not exist or does not have a value
     */
    template<typename TDataType>
    static TDataType const& GetValue(std::string const& rItemFullName)
    {
        return GetItem(rItemFullName).GetValue<TDataType>();
    }

    /**
     * @brief Retrieves a value casted to a derived type
     * @tparam TDataType Base type of the value
     * @tparam TCastType Derived type to cast to
     * @param rItemFullName Full hierarchical path of the item
     * @return Const reference to the value casted to the derived type
     */
    template<typename TDataType, typename TCastType>
    static typename std::enable_if<std::is_base_of<TDataType, TCastType>::value, TCastType>::type const GetValueAs(std::string const& rItemFullName)
    {
        return GetItem(rItemFullName).GetValueAs<TDataType, TCastType>();
    }

    /**
     * @brief Removes an item from the registry
     * @param rItemFullName Full hierarchical path of the item
     * @throws An error if the item does not exist
     */
    static void RemoveItem(std::string const& ItemName);

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Returns the number of items in the registry
     * @return Total number of items
     */
    static std::size_t size();

    /**
     * @brief Checks if an item exists in the registry
     * @param rItemFullName Full hierarchical path of the item
     * @return True if the item exists, false otherwise
     */
    static bool HasItem(std::string const& rItemFullName);

    /**
     * @brief Checks if an item has an associated value
     * @param rItemFullName Full hierarchical path of the item
     * @return True if the item has a value, false otherwise
     */
    static bool HasValue(std::string const& rItemFullName);

    /**
     * @brief Checks if an item has sub-items
     * @param rItemFullName Full hierarchical path of the item
     * @return True if the item has sub-items, false otherwise
     */
    static bool HasItems(std::string const& rItemFullName);

    ///@}
    ///@name Input and output
    ///@{


    /**
     * @brief Turn back information as a string
     * @return Information about this object
     */
    std::string Info() const;


    /**
     * @brief Prints information about this object
     * @param rOStream Output stream
     */
    void PrintInfo(std::ostream& rOStream) const;

    /**
     * @brief Prints data stored in this object
     * @param rOStream Output stream
     */
    void PrintData(std::ostream& rOStream) const;

    /**
     * @brief Converts the registry to a JSON string
     * @param Indentation Indentation string for formatting
     * @return JSON representation of the registry
     */
    std::string ToJson(std::string const& Indentation) const;

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    /// Root registry item
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

    /**
     * @brief This method throws an error message for a not found item.
     * @param rFullName The full name of the item.
     * @param rItemName The name of the item.
     * @param pCurrentItem The current item.
     */
    static void NotFoundError(
        const std::string& rFullName,
        const std::string& rItemName,
        RegistryItem* pCurrentItem
        );

    ///@}
    ///@name Private  Access
    ///@{

    /**
     * @brief Returns the root registry item
     * @return Reference to the root registry item
     */
    static RegistryItem& GetRootRegistryItem();

    /**
     * @brief Splits a full item name into a vector of strings
     * @param FullName Full hierarchical path of the item
     * @return Vector of strings representing the item path
     */
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
