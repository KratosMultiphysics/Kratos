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
#include <any>

// External includes


// Project includes
#include "includes/define.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{
///@name Kratos Classes
///@{

/// The registry item to be stored by Registry class. It is the base class for some more specific ones.
/** RegistryItem has a tree node structure and stores its name, an optional
 *  value, and an unorder_set of its sub data.
 *  This structure let us to have registry of the elements and then different
 *  registries for each elements inside it.
 *  Please note that RegistryItem stores a shared pointer to the value.
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

    /// Custom iterator with key as return type to be used in the Python export
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

        KeyReturnConstIterator()
        {}

        KeyReturnConstIterator(const BaseIterator Iterator)
            : mIterator(Iterator)
        {}

        KeyReturnConstIterator(const KeyReturnConstIterator& rIterator)
            : mIterator(rIterator.mIterator)
        {}

        ///@}
        ///@name Operators
        ///@{

        KeyReturnConstIterator& operator=(const KeyReturnConstIterator& rIterator)
        {
            this->mIterator = rIterator.mIterator;
            return *this;
        }

        const_reference operator*() const
        {
            return mIterator->first;
        }

        const_pointer operator->() const
        {
            return &(mIterator->first);
        }

        KeyReturnConstIterator& operator++()
        {
            ++mIterator;
            return *this;
        }

        KeyReturnConstIterator operator++(int)
        {
            KeyReturnConstIterator tmp(*this);
            ++(*this);
            return tmp;
        }

        bool operator==(const KeyReturnConstIterator& rIterator) const
        {
            return this->mIterator == rIterator.mIterator;
        }

        bool operator!=(const KeyReturnConstIterator& rIterator) const
        {
            return this->mIterator != rIterator.mIterator;
        }

        ///@}
    private:
        ///@name Member Variables
        ///@{

        BaseIterator mIterator;

        ///@}
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor deleted.
    RegistryItem() = delete;

    /// Constructor with the name
    RegistryItem(const std::string &rName)
        : mName(rName),
          mpValue(Kratos::make_shared<SubRegistryItemType>()),
          mValueName(GetRegistryItemType()) {}

    /// Constructor with the name and lambda
    template <typename TItemType, typename... TArgs>
    RegistryItem(
        const std::string &rName,
        const std::function<TItemType(TArgs...)> &rValue)
        : mName(rName),
          mpValue(rValue),
          mValueName(typeid(TItemType).name()) {}

    /// Constructor with the name and value
    template<class TItemType>
    RegistryItem(
        const std::string&  rName,
        const TItemType& rValue)
        : mName(rName),
          mpValue((std::function<std::shared_ptr<TItemType>()>)[=](){return Kratos::make_shared<TItemType>(rValue);})
    {
        std::stringstream buffer;
        buffer << rValue;
        mValueName = buffer.str();
    }

    /// Constructor with the name and shared ptr
    template<class TItemType>
    RegistryItem(
        const std::string&  rName,
        const shared_ptr<TItemType>& pValue)
        : mName(rName),
          mpValue((std::function<std::shared_ptr<TItemType>()>)[=](){return pValue;})
    {
        std::stringstream buffer;
        buffer << *pValue;
        mValueName = buffer.str();
    }

    // Copy constructor deleted
    RegistryItem(RegistryItem const& rOther) = delete;

    /// Destructor.
    ~RegistryItem() = default;

    /// Assignment operator deleted.
    RegistryItem& operator=(RegistryItem& rOther) = delete;

    ///@}
    ///@name Items
    ///@{

    // Adding sub value item
    template <typename TItemType, class... TArgumentsList>
    RegistryItem &AddItem(
        std::string const &ItemName,
        TArgumentsList &&...Arguments)
    {
        KRATOS_ERROR_IF(this->HasItem(ItemName))
            << "The RegistryItem '" << this->Name() << "' already has an item with name "
            << ItemName << "." << std::endl;

        if constexpr (std::is_same<TItemType, RegistryItem>::value == true) {
            using ValueType = typename std::conditional<std::is_same<TItemType, RegistryItem>::value, SubRegistryItemFunctor, SubValueItemFunctor<TItemType>>::type;

            auto insert_result = GetSubRegistryItemMap().emplace(
                std::make_pair(
                    ItemName,
                    ValueType::Create(ItemName, std::forward<TArgumentsList>(Arguments)...)));

            KRATOS_ERROR_IF_NOT(insert_result.second)
                << "Error in inserting '" << ItemName
                << "' in registry item with name '" << this->Name() << "'." << std::endl;

            return *insert_result.first->second;
        }
        else
        {
            std::function<std::shared_ptr<TItemType>()> aux_lambda = [=]() -> std::shared_ptr<TItemType> {
                return std::make_shared<TItemType>((Arguments)...);
            };

            using CallableType = typename std::conditional<std::is_same<std::function<std::shared_ptr<TItemType>()>, RegistryItem>::value, SubRegistryItemFunctor, SubValueItemFunctor<std::function<std::shared_ptr<TItemType>()>>>::type;

            auto insert_result = GetSubRegistryItemMap().emplace(
                std::make_pair(
                    ItemName,
                    std::make_shared<RegistryItem>(ItemName, aux_lambda)
                )
            );

            KRATOS_ERROR_IF_NOT(insert_result.second)
                << "Error in inserting '" << ItemName
                << "' in registry item with name '" << this->Name() << "'." << std::endl;

            // Force the evaluation of the lambda in order to set mpCallable and mpValueName
            this->EvaluateValue<TItemType>();

            return *insert_result.first->second;
        }
    }

    // Adding sub value callable item
    template <typename TItemType, class TFunctionType>
    RegistryItem &AddCallable(
        std::string const &ItemName,
        TFunctionType &&Function)
    {
        KRATOS_ERROR_IF(this->HasItem(ItemName))
            << "The RegistryItem '" << this->Name() << "' already has an item with name "
            << ItemName << "." << std::endl;

        auto insert_result = GetSubRegistryItemMap().emplace(
            std::make_pair(
                ItemName,
                std::make_shared<RegistryItem>(ItemName, Function)
            )
        );

        KRATOS_ERROR_IF_NOT(insert_result.second)
            << "Error in inserting '" << ItemName
            << "' in registry item with name '" << this->Name() << "'." << std::endl;

        return *insert_result.first->second;
    }

    ///@}
    ///@name Access
    ///@{

    SubRegistryItemType::iterator begin();

    SubRegistryItemType::const_iterator cbegin() const;

    SubRegistryItemType::iterator end();

    SubRegistryItemType::const_iterator cend() const;

    KeyReturnConstIterator KeyConstBegin() const;

    KeyReturnConstIterator KeyConstEnd() const;

    const std::string& Name() const  { return mName; }

    RegistryItem const& GetItem(std::string const& rItemName) const;

    RegistryItem& GetItem(std::string const& rItemName);

    template <typename TDataType>
    void EvaluateValue()
    {
         // Assign callable value
        using TFunctionType = std::function<std::shared_ptr<TDataType>()>;
        TFunctionType func = std::any_cast<TFunctionType>(mpValue);
        mpCallable = func();

        // Set value name
        std::stringstream buffer;
        buffer << *(std::any_cast<std::shared_ptr<TDataType>>(mpCallable));
        mValueName = buffer.str();
    }

    template <typename TDataType>
    TDataType const &GetValue()
    {
        KRATOS_TRY

        // This is executed the first time we access the GetValue for this item
        using TFunctionType = std::function<std::shared_ptr<TDataType>()>;
        if (std::any_cast<std::shared_ptr<TDataType>>(&mpCallable) == nullptr) {
            this->EvaluateValue<TDataType>();
        }

        return *(std::any_cast<std::shared_ptr<TDataType>>(mpCallable));

        KRATOS_CATCH("");
    }

    void RemoveItem(std::string const& rItemName);

    ///@}
    ///@name Inquiry
    ///@{

    std::size_t size();

    bool HasValue() const;

    bool HasItems() const;

    bool HasItem(std::string const& rItemName) const;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    std::string ToJson(std::string const& rTabSpacing = "", const std::size_t Level = 0) const;

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    std::string mName;
    std::any mpValue;
    std::any mpCallable = nullptr;
    std::string mValueName;

    ///@}
    ///@name Private operations
    ///@{

    std::string GetRegistryItemType() const
    {
        return mpValue.type().name();
    }

    ///@}
    ///@name Private classes
    ///@{

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

    template<typename TItemType>
    class SubValueItemFunctor
    {
    public:
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

    std::string GetValueString() const;

    SubRegistryItemType& GetSubRegistryItemMap();

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
