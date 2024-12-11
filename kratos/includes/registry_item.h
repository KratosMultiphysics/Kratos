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
//                   Carlos Roig
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <unordered_map>

#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
    #include <boost/any.hpp>
#else
    #include <any>
#endif

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{
///@name Kratos Classes
///@{

/**
 * @class RegistryItem
 * @ingroup KratosCore
 * @brief Represents a node in the registry hierarchy.
 * A `RegistryItem` can store a value or act as a container for sub-items, forming a tree structure.
 * It is the base class for more specific registry item types.
 * @details
 * - Each `RegistryItem` has a unique name and optionally holds a value.
 * - Sub-items are stored in an unordered map for quick access.
 * - Supports serialization to JSON for debugging and inspection.
 * @author Pooyan Dadvand
 *         Ruben Zorrilla
 *         Carlos Roig
 * @see Registry
 */
class KRATOS_API(KRATOS_CORE) RegistryItem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RegistryItem
    KRATOS_CLASS_POINTER_DEFINITION(RegistryItem);

    /// Subregistry item type definition
    using SubRegistryItemType = std::unordered_map<std::string, Kratos::shared_ptr<RegistryItem>>;

    /// Pointer definition of SubRegistryItemType
    using SubRegistryItemPointerType = Kratos::shared_ptr<SubRegistryItemType>;

    /**
     * @class KeyReturnConstIterator
     * @brief Custom iterator for registry items returning keys.
     * @details This iterator is designed for Python bindings, allowing easy iteration over item keys.
     */
    class KeyReturnConstIterator
    {
    public:
        ///@name Type Definitions
        ///@{

        using BaseIterator      = SubRegistryItemType::const_iterator;
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = BaseIterator::difference_type;
        using value_type        = SubRegistryItemType::key_type;
        using const_pointer     = const value_type*;
        using const_reference   = const value_type&;

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor
        KeyReturnConstIterator() {}

        /// Constructor with base iterator
        KeyReturnConstIterator(const BaseIterator Iterator)
            : mIterator(Iterator) {}

        /// Copy constructor
        KeyReturnConstIterator(const KeyReturnConstIterator& rIterator)
            : mIterator(rIterator.mIterator) {}

        ///@}
        ///@name Operators
        ///@{

        /// Assignment operator
        KeyReturnConstIterator& operator=(const KeyReturnConstIterator& rIterator)
        {
            this->mIterator = rIterator.mIterator;
            return *this;
        }

        /// Dereference operator
        const_reference operator*() const
        {
            return mIterator->first;
        }

        /// Pointer dereference operator
        const_pointer operator->() const
        {
            return &(mIterator->first);
        }

        /// Pre-increment operator
        KeyReturnConstIterator& operator++()
        {
            ++mIterator;
            return *this;
        }

        /// Post-increment operator
        KeyReturnConstIterator operator++(int)
        {
            KeyReturnConstIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        /// Equality comparison operator
        bool operator==(const KeyReturnConstIterator& rIterator) const
        {
            return this->mIterator == rIterator.mIterator;
        }

        /// Inequality comparison operator
        bool operator!=(const KeyReturnConstIterator& rIterator) const
        {
            return this->mIterator != rIterator.mIterator;
        }

        ///@}
    private:
        ///@name Member Variables
        ///@{

        BaseIterator mIterator; /// Internal iterator

        ///@}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor is deleted
    RegistryItem() = delete;

    /**
     * @brief Constructor with name
     * @param rName The name of the registry item
     */
    RegistryItem(const std::string& rName)
        : mName(rName),
          mpValue(Kratos::make_shared<SubRegistryItemType>()),
          mGetValueStringMethod(&RegistryItem::GetRegistryItemType) {}

    /**
     * @brief Constructor with name and lambda function for value
     * @tparam TItemType Type of the item
     * @tparam TArgs Argument types for the lambda
     * @param rName The name of the registry item
     * @param rValue Lambda function to initialize the value
     */
    template <typename TItemType, typename... TArgs>
    RegistryItem(
        const std::string &rName,
        const std::function<std::shared_ptr<TItemType>(TArgs...)> &rValue)
        : mName(rName),
          mpValue(rValue()),
          mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>) {}

    /**
     * @brief Constructor with name and value
     * @tparam TItemType Type of the value
     * @param rName The name of the registry item
     * @param rValue The value to store
     */
    template<class TItemType>
    RegistryItem(
        const std::string& rName,
        const TItemType& rValue)
        : mName(rName),
          mpValue(Kratos::make_shared<TItemType>(rValue)),
          mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>) {}

    /**
     * @brief Constructor with name and shared pointer to value
     * @tparam TItemType Type of the value
     * @param rName The name of the registry item
     * @param pValue Shared pointer to the value
     */
    template<class TItemType>
    RegistryItem(
        const std::string& rName,
        const shared_ptr<TItemType>& pValue)
        : mName(rName),
          mpValue(pValue),
          mGetValueStringMethod(&RegistryItem::GetItemString<TItemType>) {}

    /// Copy constructor is deleted
    RegistryItem(RegistryItem const& rOther) = delete;

    /// Destructor
    ~RegistryItem() = default;

    /// Assignment operator is deleted
    RegistryItem& operator=(RegistryItem& rOther) = delete;

    ///@}
    ///@name Items
    ///@{

    /**
     * @brief Adds a sub-item to the registry item
     * @tparam TItemType Type of the sub-item
     * @tparam TArgumentsList Variadic template arguments for the sub-item constructor
     * @param ItemName Name of the sub-item
     * @param Arguments Arguments for the sub-item constructor
     * @return Reference to the added sub-item
     */
    template<typename TItemType, class... TArgumentsList>
    RegistryItem& AddItem(
        std::string const& ItemName,
        TArgumentsList&&... Arguments)
    {
        KRATOS_ERROR_IF(this->HasItem(ItemName))
            << "The RegistryItem '" << this->Name() << "' already has an item with name "
            << ItemName << "." << std::endl;

        using ValueType = typename std::conditional<std::is_same<TItemType, RegistryItem>::value, SubRegistryItemFunctor, SubValueItemFunctor<TItemType>>::type;

        auto insert_result = GetSubRegistryItemMap().emplace(
            std::make_pair(
                ItemName,
                ValueType::Create(ItemName, std::forward<TArgumentsList>(Arguments)...)
                )
            );

        KRATOS_ERROR_IF_NOT(insert_result.second)
            << "Error in inserting '" << ItemName
            << "' in registry item with name '" << this->Name() << "'." << std::endl;

        return *insert_result.first->second;
    }

    ///@}
    ///@name Iterators
    ///@{

    /**
     * @brief Returns an iterator to the beginning of sub-items
     * @return Iterator pointing to the first sub-item
     * @throws An error if the current item holds a value instead of sub-items
     */
    SubRegistryItemType::iterator begin();

    /**
     * @brief Returns a constant iterator to the beginning of sub-items
     * @return Constant iterator pointing to the first sub-item
     * @throws An error if the current item holds a value instead of sub-items
     */
    SubRegistryItemType::const_iterator cbegin() const;

    /**
     * @brief Returns an iterator to the end of sub-items
     * @return Iterator pointing past the last sub-item
     * @throws An error if the current item holds a value instead of sub-items
     */
    SubRegistryItemType::iterator end();

    /**
     * @brief Returns a constant iterator to the end of sub-items
     * @return Constant iterator pointing past the last sub-item
     * @throws An error if the current item holds a value instead of sub-items
     */
    SubRegistryItemType::const_iterator cend() const;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Returns an iterator to the beginning of sub-item keys
     * @return Iterator for sub-item keys
     */
    KeyReturnConstIterator KeyConstBegin() const;

    /**
     * @brief Returns an iterator to the end of sub-item keys
     * @return Iterator past the last sub-item key
     */
    KeyReturnConstIterator KeyConstEnd() const;

    /**
     * @brief Returns the name of the current registry item
     * @return Name of the registry item
     */
    const std::string& Name() const { return mName; }

    /**
     * @brief Retrieves a sub-item by name
     * @param rItemName Name of the sub-item to retrieve
     * @return Constant reference to the requested sub-item
     * @throws An error if the sub-item does not exist
     */
    RegistryItem const& GetItem(std::string const& rItemName) const;

    /**
     * @brief Retrieves a sub-item by name
     * @param rItemName Name of the sub-item to retrieve
     * @return Reference to the requested sub-item
     * @throws An error if the sub-item does not exist
     */
    RegistryItem& GetItem(std::string const& rItemName);

    /**
     * @brief Retrieves the value stored in this registry item
     * @tparam TDataType Type of the value to retrieve
     * @return Constant reference to the value
     * @throws An error if the item does not store a value or if the type does not match
     */
    template<typename TDataType>
    TDataType const& GetValue() const
    {
        KRATOS_TRY

#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
        return *(boost::any_cast<std::shared_ptr<TDataType>>(mpValue));
#else
        return *(std::any_cast<std::shared_ptr<TDataType>>(mpValue));
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Retrieves a value casted to a derived type
     * @tparam TDataType Base type of the stored value
     * @tparam TCastType Derived type to cast the value to
     * @return Constant reference to the value casted to the derived type
     * @throws An error if the cast fails or if the type does not match
     */
    template<typename TDataType, typename TCastType>
    TCastType const& GetValueAs() const
    {
        KRATOS_TRY

#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
        return *std::dynamic_pointer_cast<TCastType>(boost::any_cast<std::shared_ptr<TDataType>>(mpValue));
#else
        return *std::dynamic_pointer_cast<TCastType>(std::any_cast<std::shared_ptr<TDataType>>(mpValue));
#endif

        KRATOS_CATCH("");
    }

    /**
     * @brief Removes a sub-item by name
     * @param rItemName Name of the sub-item to remove
     * @throws An error if the sub-item does not exist
     */
    void RemoveItem(std::string const& rItemName);

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Returns the number of sub-items
     * @return Number of sub-items
     */
    std::size_t size();

    /**
     * @brief Checks if the item holds a value
     * @return True if the item has a value, false otherwise
     */
    bool HasValue() const;

    /**
     * @brief Checks if the item has sub-items
     * @return True if the item has sub-items, false otherwise
     */
    bool HasItems() const;

    /**
     * @brief Checks if a sub-item exists
     * @param rItemName Name of the sub-item to check
     * @return True if the sub-item exists, false otherwise
     */
    bool HasItem(std::string const& rItemName) const;

    /**
     * @brief Checks if the type of the stored value matches a given type
     * @tparam OtherType Type to compare against
     * @param rOther A dummy variable of the type to compare against
     * @return True if the types match, false otherwise
     */
    template<class OtherType>
    bool IsSameType(const OtherType& rOther) const {
        try {
            GetValue<OtherType>();
        } catch (...) {
            return false;
        }
        return true;
    }

    ///@}
    ///@name Input and Output
    ///@{

    /**
     * @brief Returns a string representation of the item
     * @return A string describing the registry item
     */
    std::string Info() const;

    /**
     * @brief Prints detailed information about the item
     * @param rOStream Output stream to write to
     */
    void PrintInfo(std::ostream& rOStream) const;

    /**
     * @brief Prints the item's data
     * @param rOStream Output stream to write to
     */
    void PrintData(std::ostream& rOStream) const;

    /**
     * @brief Serializes the item and its sub-items to JSON
     * @param rTabSpacing String to use for indentation
     * @param Level Depth of the current item in the hierarchy
     * @return A JSON string representing the item
     */
    std::string ToJson(std::string const& rTabSpacing = "", const std::size_t Level = 0) const;

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    std::string mName; /// Name of the registry item

#if defined (__GNUC__) && __GNUC__ <= 8 && __GNUC_MINOR__ <= 3
    boost::any mpValue; /// Value stored in the item (if any)
#else
    std::any mpValue; /// Value stored in the item (if any)
#endif

    std::string (RegistryItem::*mGetValueStringMethod)() const; /// Pointer to a function to get the string representation of the value

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the type of the value as a string
     * @return String representation of the type of the stored value
     */
    std::string GetRegistryItemType() const
    {
        return mpValue.type().name();
    }

    /**
     * @brief Converts the stored value to a string
     * @tparam TItemType Type of the stored value
     * @return String representation of the value
     */
    template<class TItemType>
    std::string GetItemString() const
    {
        std::stringstream buffer;
        buffer << this->GetValue<TItemType>();
        return buffer.str();
    }

    ///@}

    ///@name Private Classes
    ///@{

    /**
     * @class SubRegistryItemFunctor
     * @brief Functor for creating sub-items with no values
     */
    class SubRegistryItemFunctor
    {
    public:
        template<class... TArgumentsList>
        static inline RegistryItem::Pointer Create(
            std::string const& ItemName,
            TArgumentsList&&... Arguments)
        {
            return Kratos::make_shared<RegistryItem>(ItemName);
        }
    };

    /**
     * @class SubValueItemFunctor
     * @brief Functor for creating sub-items with values
     * @tparam TItemType Type of the value
     */
    template<typename TItemType>
    class SubValueItemFunctor
    {
    public:
        template<class... TArgumentsList, class TFunctionType = std::function<std::shared_ptr<TItemType>(TArgumentsList...)>>
        static inline RegistryItem::Pointer Create(
            std::string const& ItemName,
            TFunctionType && Function)
        {
            return Kratos::make_shared<RegistryItem>(ItemName, std::forward<TFunctionType>(Function));
        }

        template<class... TArgumentsList>
        static inline RegistryItem::Pointer Create(
            std::string const& ItemName,
            TArgumentsList&&... Arguments)
        {
            return Kratos::make_shared<RegistryItem>(ItemName, Kratos::make_shared<TItemType>(std::forward<TArgumentsList>(Arguments)...));
        }

    };

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * @brief Returns the string representation of the stored value
     * @details This method converts the stored value to a string, regardless of its type. It uses the internal
     * function pointer `mGetValueStringMethod` to handle the conversion.
     * @return String representation of the value
     * @throws An error if no value is stored or if the conversion fails
     */
    std::string GetValueString() const;

    /**
     * @brief Returns a reference to the map of sub-items
     * @details This method provides access to the unordered map that stores the sub-items of this registry item.
     * It ensures that the current registry item is a container (i.e., does not hold a value).
     * @return Reference to the unordered map of sub-items
     * @throws An error if the current registry item holds a value
     */
    SubRegistryItemType& GetSubRegistryItemMap();

    /**
     * @brief Returns a constant reference to the map of sub-items
     * @details This method provides read-only access to the unordered map that stores the sub-items of this registry item.
     * It ensures that the current registry item is a container (i.e., does not hold a value).
     * @return Constant reference to the unordered map of sub-items
     * @throws An error if the current registry item holds a value
     */
    SubRegistryItemType& GetSubRegistryItemMap() const;

    ///@}

}; // Class RegistryItem

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
