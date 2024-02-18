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
//

#pragma once

// System includes
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <utility>

// External includes
#include <boost/iterator/indirect_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "containers/set_identity_function.h"

namespace Kratos
{

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
 * @ingroup KratosCore
 * @class PointerVectorSet
 * @brief A sorted associative container similar to an STL set, but uses a vector to store pointers to its data.
 * @details The `PointerVectorSet` is a sorted associative container that behaves like an STL set but employs a vector
 * to store pointers to its data elements. Many of its methods are inspired by the Boost ptr_container library,
 * with modifications to support shared pointers.
 * @note Unlike the Boost counterpart, this container does not manage memory deallocation automatically. It relies on
 * the use of smart pointers or manual memory management for resource cleanup.
 * @author Pooyan Dadvand
 */
template<class TDataType,
         class TGetKeyType = SetIdentityFunction<TDataType>,
         class TCompareType = std::less<decltype(std::declval<TGetKeyType>()(std::declval<TDataType>()))>,
         class TEqualType = std::equal_to<decltype(std::declval<TGetKeyType>()(std::declval<TDataType>()))>,
         class TPointerType = typename TDataType::Pointer,
         class TContainerType = std::vector<TPointerType>>
class PointerVectorSet final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointerVectorSet
    KRATOS_CLASS_POINTER_DEFINITION(PointerVectorSet);

    /// Key type for searching in this container.
    using key_type = typename std::remove_reference<decltype(std::declval<TGetKeyType>()(std::declval<TDataType>()))>::type;

    // Data type stored in this container.
    using data_type = TDataType;
    using value_type = TDataType;
    using key_compare = TCompareType;
    using pointer = TPointerType;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    using ContainerType = TContainerType;

    /// @}
    /// @name Iterators
    /// @{
    using iterator = boost::indirect_iterator<typename TContainerType::iterator>;
    using const_iterator = boost::indirect_iterator<typename TContainerType::const_iterator>;
    using reverse_iterator = boost::indirect_iterator<typename TContainerType::reverse_iterator>;
    using const_reverse_iterator = boost::indirect_iterator<typename TContainerType::const_reverse_iterator>;

    /// @}
    /// @name Other definitions
    /// @{
    using size_type = typename TContainerType::size_type;
    using ptr_iterator = typename TContainerType::iterator;
    using ptr_const_iterator = typename TContainerType::const_iterator;
    using ptr_reverse_iterator = typename TContainerType::reverse_iterator;
    using ptr_const_reverse_iterator = typename TContainerType::const_reverse_iterator;
    using difference_type = typename TContainerType::difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointerVectorSet()
        : mData()
    {}

    /**
    * @brief Constructs a PointerVectorSet from a range of elements.
    * @details This constructor initializes a PointerVectorSet with elements in the range [First, Last).
    * @tparam TInputIteratorType The type of the input iterator.
    * @param First An input iterator pointing to the beginning of the range.
    * @param Last An input iterator pointing to the end of the range.
    * @param NewMaxBufferSize The maximum buffer size (default is 1).
    */
    template <class TInputIteratorType>
    PointerVectorSet(
        TInputIteratorType First,
        TInputIteratorType Last)
    {
        mData.reserve(std::distance(First, Last));
        insert(First, Last);
    }

    /**
     * @brief Copy constructor for PointerVectorSet.
     * @param rOther The PointerVectorSet to copy from.
     */
    PointerVectorSet(const PointerVectorSet& rOther)
        : mData(rOther.mData)
    {}

    /**
     * @brief Constructs a PointerVectorSet from a container.
     * @details This constructor initializes a PointerVectorSet with elements from a container.
     * @param rContainer The container to copy elements from.
     */
    explicit PointerVectorSet(const TContainerType& rContainer)
        : mData(rContainer)
    {
        Conform();
    }

    /// Destructor.
    ~PointerVectorSet() {}

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment operator for PointerVectorSet.
     * @details Assigns the contents of another PointerVectorSet to this one.
     * @param rOther The PointerVectorSet to assign from.
     * @return A reference to the modified PointerVectorSet.
     */
    PointerVectorSet& operator=(const PointerVectorSet& rOther)
    {
        mData = rOther.mData;
        return *this;
    }

    /**
     * @brief Accesses an element by key and returns a reference.
     * @details This operator allows you to access an element in the PointerVectorSet by its key. If the key is not found,
     * a new element with the provided key is inserted into the set.
     * @param Key The key of the element to access or insert.
     * @return A reference to the accessed or newly inserted element.
     */
    TDataType& operator[](const key_type& Key)
    {
        auto itr_pos = std::lower_bound(mData.begin(), mData.end(), Key, CompareKey());
        if (itr_pos == mData.end()) {
            // insert a new element with id.
            mData.push_back(TPointerType(new TDataType(Key)));
            return back();
        } else if (EqualKeyTo(Key)(*itr_pos)) {
            // already found existing element with the same key, hence returning the existing element.
            return **itr_pos;
        } else {
            // insert the new value before the itr_pos.
            return **(mData.insert(itr_pos, TPointerType(new TDataType(Key))));
        }
    }

    /**
     * @brief Function for inserting or retrieving a pointer associated with a key.
     * @details This function allows you to insert or retrieve a pointer associated with a key in the set.
     * If the key already exists in the set, it returns the corresponding pointer. If not, it inserts
     * the key and a new pointer into the set and returns the newly inserted pointer.
     * @param Key The key for which you want to insert or retrieve a pointer.
     * @return A reference to the pointer associated with the given key.
     */
    pointer& operator()(const key_type& Key)
    {
        auto itr_pos = std::lower_bound(mData.begin(), mData.end(), Key, CompareKey());
        if (itr_pos == mData.end()) {
            // insert a new element with id.
            mData.push_back(TPointerType(new TDataType(Key)));
            return mData.back();
        } else if (EqualKeyTo(Key)(*itr_pos)) {
            // already found existing element with the same key, hence returning the existing element.
            return *itr_pos;
        } else {
            // insert the new value before the itr_pos.
            return *(mData.insert(itr_pos, TPointerType(new TDataType(Key))));
        }
    }

    /**
     * @brief Equality comparison operator for two PointerVectorSet objects.
     * @details This operator checks if two PointerVectorSet objects are equal by comparing their sizes
     * and the equality of their elements using the EqualKeyTo comparison function.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @param r The PointerVectorSet to compare with.
     * @return true if the two sets are equal, false otherwise.
     */
    bool operator==(const PointerVectorSet& r) const noexcept
    {
        assert( !empty() );
        if (size() != r.size())
            return false;
        else
            return std::equal(mData.begin(), mData.end(), r.mData.begin(), EqualKeyTo());
    }

    /**
     * @brief Less than comparison operator for two PointerVectorSet objects.
     * @details This operator checks if one PointerVectorSet is less than another by comparing their
     * elements using the CompareKey comparison function in a lexicographical order.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @param r The PointerVectorSet to compare with.
     * @return true if this set is less than r, false otherwise.
     */
    bool operator<(const PointerVectorSet& r) const noexcept
    {
        assert( !empty() );
        return std::lexicographical_compare(mData.begin(), mData.end(), r.mData.begin(), r.mData.end(), CompareKey());
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Returns an iterator pointing to the beginning of the container.
     * @return An iterator pointing to the beginning of the container.
     */
    iterator begin()
    {
        return iterator( mData.begin() );
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the container.
     * @return A constant iterator pointing to the beginning of the container.
     */
    const_iterator begin() const
    {
        return const_iterator( mData.begin() );
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the container.
     * @return A constant iterator pointing to the beginning of the container.
     */
    const_iterator cbegin()
    {
        return const_iterator(mData.begin());
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the container.
     * @return A constant iterator pointing to the beginning of the container.
     */
    const_iterator cbegin() const
    {
        return const_iterator(mData.begin());
    }

    /**
     * @brief Returns an iterator pointing to the end of the container.
     * @return An iterator pointing to the end of the container.
     */
    iterator end()
    {
        return iterator( mData.end() );
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the container.
     * @return A constant iterator pointing to the end of the container.
     */
    const_iterator end() const
    {
        return const_iterator( mData.end() );
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the container.
     * @return A constant iterator pointing to the end of the container.
     */
    const_iterator cend()
    {
        return const_iterator(mData.end());
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the container.
     * @return A constant iterator pointing to the end of the container.
     */
    const_iterator cend() const
    {
        return const_iterator(mData.end());
    }

    /**
     * @brief Returns a reverse iterator pointing to the beginning of the container.
     * @return A reverse iterator pointing to the beginning of the container.
     */
    reverse_iterator rbegin()
    {
        return reverse_iterator( mData.rbegin() );
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the beginning of the container.
     * @return A constant reverse iterator pointing to the beginning of the container.
     */
    const_reverse_iterator rbegin() const
    {
        return const_reverse_iterator( mData.rbegin() );
    }

    /**
     * @brief Returns a reverse iterator pointing to the end of the container.
     * @return A reverse iterator pointing to the end of the container.
     */
    reverse_iterator rend()
    {
        return reverse_iterator( mData.rend() );
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the end of the container.
     * @return A constant reverse iterator pointing to the end of the container.
     */
    const_reverse_iterator rend() const
    {
        return const_reverse_iterator( mData.rend() );
    }

    /**
     * @brief Returns an iterator pointing to the beginning of the underlying data container.
     * @return An iterator pointing to the beginning of the underlying data container.
     */
    ptr_iterator ptr_begin()
    {
        return mData.begin();
    }

    /**
     * @brief Returns a constant iterator pointing to the beginning of the underlying data container.
     * @return A constant iterator pointing to the beginning of the underlying data container.
     */
    ptr_const_iterator ptr_begin() const
    {
        return mData.begin();
    }

    /**
     * @brief Returns an iterator pointing to the end of the underlying data container.
     * @return An iterator pointing to the end of the underlying data container.
     */
    ptr_iterator ptr_end()
    {
        return mData.end();
    }

    /**
     * @brief Returns a constant iterator pointing to the end of the underlying data container.
     * @return A constant iterator pointing to the end of the underlying data container.
     */
    ptr_const_iterator ptr_end() const
    {
        return mData.end();
    }

    /**
     * @brief Returns a reverse iterator pointing to the beginning of the underlying data container.
     * @return A reverse iterator pointing to the beginning of the underlying data container.
     */
    ptr_reverse_iterator ptr_rbegin()
    {
        return mData.rbegin();
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the beginning of the underlying data container.
     * @return A constant reverse iterator pointing to the beginning of the underlying data container.
     */
    ptr_const_reverse_iterator ptr_rbegin() const
    {
        return mData.rbegin();
    }

    /**
     * @brief Returns a reverse iterator pointing to the end of the underlying data container.
     * @return A reverse iterator pointing to the end of the underlying data container.
     */
    ptr_reverse_iterator ptr_rend()
    {
        return mData.rend();
    }

    /**
     * @brief Returns a constant reverse iterator pointing to the end of the underlying data container.
     * @return A constant reverse iterator pointing to the end of the underlying data container.
     */
    ptr_const_reverse_iterator ptr_rend() const
    {
        return mData.rend();
    }

    /**
     * @brief Returns a reference to the first element in the container.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @return A reference to the first element in the container.
     */
    reference front() noexcept
    {
        assert( !empty() );
        return *(mData.front());
    }

    /**
     * @brief Returns a constant reference to the first element in the container.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @return A constant reference to the first element in the container.
     */
    const_reference front() const noexcept
    {
        assert( !empty() );
        return *(mData.front());
    }

    /**
     * @brief Returns a reference to the last element in the container.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @return A reference to the last element in the container.
     */
    reference back() noexcept
    {
        assert( !empty() );
        return *(mData.back());
    }

    /**
     * @brief Returns a constant reference to the last element in the container.
     * @note This function is marked as "noexcept," and it asserts that the container is not empty.
     * @return A constant reference to the last element in the container.
     */
    const_reference back() const noexcept
    {
        assert( !empty() );
        return *(mData.back());
    }

    /**
     * @brief Returns the number of elements in the container.
     * @return The number of elements in the container.
     */
    size_type size() const
    {
        return mData.size();
    }

    /**
     * @brief Returns the maximum possible number of elements the container can hold.
     * @return The maximum possible number of elements the container can hold.
     */
    size_type max_size() const
    {
        return mData.max_size();
    }

    /**
     * @brief Returns the key comparison function used for ordering elements in the container.
     * @details This function returns an instance of the key comparison function (TCompareType) used for ordering elements in the container.
     * @return The key comparison function.
     */
    key_compare key_comp() const
    {
        return TCompareType();
    }

    /**
     * @brief Swaps the contents of this PointerVectorSet with another.
     * @details This function swaps the contents of this PointerVectorSet with another set, including
     * @param rOther The other PointerVectorSet to swap with.
     */
    void swap(PointerVectorSet& rOther)
    {
        mData.swap(rOther.mData);
    }

    /**
     * @brief Adds a pointer to the end of the set.
     * @details This function appends a given pointer to the end of the set.
     * @param x The pointer to be added to the end of the set.
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use insert instead)") void push_back(TPointerType x)
    {
        insert(x);
    }

    /**
     * @brief Removes the last element from the set.
     * @details This function removes the last element (pointer) from the set.
     */
    void pop_back()
    {
        mData.pop_back();
    }

    /**
     * @brief Inserts a pointer.
     * @details This function inserts a given pointer such that the resulting PointerVectorSet
     * is kept sorted. If there exists already a pointer with a key same as the key of the value, then
     * this will return iterator of that existing pointer (The value will not be inserted.)
     * @param value The pointer to be inserted.
     * @return An iterator pointing to the inserted element.
     */
    iterator insert(const TPointerType& value)
    {
        auto itr_pos = std::lower_bound(mData.begin(), mData.end(), KeyOf(*value), CompareKey());
        if (itr_pos == mData.end()) {
            // the position to insert is at the end.
            mData.push_back(value);
            return iterator(mData.end() - 1);
        } else if (EqualKeyTo(KeyOf(*value))(*itr_pos)) {
            // already found existing element with the same key, hence returning the existing element.
            return iterator(itr_pos);
        } else {
            // insert the new value before the itr_pos.
            return mData.insert(itr_pos, value);
        }
    }

    /**
     * @brief Inserts a pointer at the specified position.
     * @details This function inserts a given pointer. If the given position_hint is valid, then
     * it uses that to insert the value, otherwise the position_hint is discarded to maintain the dataset
     * sorted. If there is an existing element with the same key as in the value, then an iterator for
     * the existing element is returned.
     * @param position_hint An iterator pointing to the position where the pointer may be inserted.
     * @param value The pointer to be inserted.
     * @return An iterator pointing to the inserted element.
     */
    iterator insert(const_iterator position_hint, const TPointerType& value)
    {
        if (empty()) {
            // the dataset is empty. So use push back.
            mData.push_back(value);
            return iterator(mData.end() - 1);
        } else if (position_hint == cend()) {
            // trying to insert at the end.
            if (KeyOf(position_hint - 1) < KeyOf(*value)) {
                // key at the position hint is less than the value of key. Hence position hint
                // is valid. So using the push back.
                mData.push_back(value);
                return iterator(mData.end() - 1);
            } else {
                // given position is invalid. Hence, discarding the hint.
                return insert(value);
            }
        } else if (EqualKeyTo(KeyOf(position_hint))(*ptr_begin())) {
            // trying to insert at the front.
            if (KeyOf(*value) < KeyOf(position_hint)) {
                // key at the position hint is greater than the value of key. Hence position hint
                // is valid. So using the push back.
                return mData.insert(mData.begin(), value);
            } else {
                // given position is invalid. Hence, discarding the hint.
                return insert(value);
            }
        } else {
            // trying to insert at an arbitrary position.
            if (KeyOf(*value) < KeyOf(position_hint) && (KeyOf(position_hint - 1) < KeyOf(*value))) {
                return mData.insert(mData.begin() + (position_hint - cbegin()), value);
            } else {
                // given position is invalid. Hence, discarding the hint.
                return insert(value);
            }
        }
    }

    /**
     * @brief Insert elements from a range of iterators.
     * @details This function inserts element pointers from a range defined by the iterators `first` and `last`
     * into the set. This will not insert any elements in the range, if there exists an element with a key
     * which is equal to an element's key in the input range.
     * @param first An input iterator pointing to the beginning of the range to insert.
     * @param last An input iterator pointing to the end of the range to insert.
     */
    template <class InputIterator>
    void insert(InputIterator first, InputIterator last)
    {
        // first sorts the input iterators and make the input unique.
        std::sort(first, last, CompareKey());
        auto new_last = std::unique(first, last, EqualKeyTo());

        if (empty()) {
            for (; first != new_last; ++first) {
                mData.push_back(*first);
            }
        } else {
            auto p_current_itr = mData.begin();
            // now add the new elements
            for (; first != new_last; ++first) {
                // find the lower bound element.
                p_current_itr = std::lower_bound(p_current_itr, mData.end(), KeyOf(**first), CompareKey());
                if (!EqualKeyTo(KeyOf(**first))(*p_current_itr)) {
                    p_current_itr = mData.insert(p_current_itr, *first);
                }
            }
        }
    }

    /**
     * @brief Erase an element at the specified position.
     * @details This function erases the element at the specified position. If the provided position is equal to `end()`,
     * it returns `end()`.
     * @param pos An iterator pointing to the position of the element to erase.
     * @return An iterator pointing to the element following the erased element, or `end()` if the
     *         provided position was equal to `end()`.
     */
    iterator erase(iterator pos)
    {
        if (pos.base() == mData.end())
            return mData.end();
        iterator new_end = iterator(mData.erase(pos.base()));
        return new_end;
    }

    /**
     * @brief Erase a range of elements defined by iterators.
     * @details This function erases a range of elements defined by the iterators `first` and `last`.
     * @param first An iterator pointing to the beginning of the range to erase.
     * @param last An iterator pointing to the end of the range to erase.
     * @return An iterator pointing to the element following the last erased element.
     */
    iterator erase(iterator first, iterator last)
    {
        iterator new_end = iterator(mData.erase(first.base(), last.base()));
        return new_end;
    }

    /**
     * @brief Erase an element with the specified key.
     * @details This function erases an element with the specified key by first finding the element using
     * the `find` function and then erasing it using the `erase(iterator)` function.
     * @param k The key of the element to erase.
     * @return An iterator pointing to the element following the erased element, or `end()` if the
     *         element with the specified key was not found.
     */
    iterator erase(const key_type& k)
    {
        return erase(find(k));
    }

    /**
     * @brief Clear the set, removing all elements.
     * @details This function clears the set by removing all elements
    */
    void clear()
    {
        mData.clear();
    }

    /**
     * @brief Find an element with the specified key.
     * @details This function searches for an element with the specified key in the set. If the element is found,
     * it returns an iterator to the found element. If the element is not found, it returns an iterator
     * pointing to the end of the set.
     * @param Key The key to search for.
     * @return An iterator pointing to the found element or the end of the set if not found.
     */
    iterator find(const key_type& Key)
    {
        auto itr = std::lower_bound(mData.begin(), mData.end(), Key, CompareKey());
        if (itr != mData.end() && !EqualKeyTo(Key)(*itr)) {
            return iterator(mData.end());
        } else {
            return iterator(itr);
        }
    }

    /**
     * @brief Find an element with the specified key (const version).
     * @details This function is a const version of find() and searches for an element with the specified key
     * in the set. If the element is found, it returns a const_iterator to the found element. If the
     * element is not found, it returns a const_iterator pointing to the end of the set.
     * @param Key The key to search for.
     * @return A const_iterator pointing to the found element or the end of the set if not found.
     */
    const_iterator find(const key_type& Key) const
    {
        auto itr = std::lower_bound(mData.begin(), mData.end(), Key, CompareKey());
        if (itr != mData.end() && !EqualKeyTo(Key)(*itr)) {
            return const_iterator(mData.end());
        } else {
            return const_iterator(itr);
        }
    }

    /**
     * @brief Reserves memory for a specified number of elements.
     * @details This function reserves memory in the underlying data container for a specified number of elements.
     * @param reservedsize The number of elements to reserve memory for.
     */
    void reserve(int reservedsize)
    {
        mData.reserve(reservedsize);
    }

    /**
     * @brief Get the current capacity of the underlying data container.
     * @details This function returns the current capacity of the underlying data container.
     * @return The current capacity of the data container.
     */
    int capacity()
    {
        return mData.capacity();
    }

    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use Conform instead)") void Sort()
    {
        Conform();
    }

    KRATOS_DEPRECATED_MESSAGE("This is legacy version (use Conform instead)") void Unique()
    {
        Conform();
    }

    void Conform()
    {
        std::sort(mData.begin(), mData.end(), CompareKey());
        auto new_last = std::unique(mData.begin(), mData.end(), EqualKeyTo());
        mData.erase(new_last, mData.end());
    }

    ///@}
    ///@name Access
    ///@{

    /** Gives a reference to underly normal container. */
    TContainerType& GetContainer()
    {
        return mData;
    }

    /** Gives a constant reference to underly normal container. */
    const TContainerType& GetContainer() const
    {
        return mData;
    }

    /**
     * @brief Get the maximum size of buffer used in the container
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version. This is a blank call") size_type GetMaxBufferSize() const
    {
        return 0;
    }

    /**
     * @brief Set the maximum size of buffer used in the container.
     * @details This container uses a buffer which keep data unsorted. After buffer size arrived to the MaxBufferSize it will sort all container and empties buffer.
     * @param NewSize Is the new buffer maximum size.
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version. This is a blank call") void SetMaxBufferSize(const size_type NewSize)
    {
    }

    /**
     * @brief Get the sorted part size of buffer used in the container.
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version. This is a blank call") size_type GetSortedPartSize() const
    {
        return 0;
    }

    /**
     * @brief Set the sorted part size of buffer used in the container.
     * @param NewSize Is the new buffer maximum size.
     */
    KRATOS_DEPRECATED_MESSAGE("This is legacy version. This is a blank call") void SetSortedPartSize(const size_type NewSize)
    {
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Check if the data container is empty.
     * @details This function checks if the data container, represented by the member variable mData, is empty.
     * @return True if the data container is empty, false otherwise.
     */
    bool empty() const
    {
        return mData.empty();
    }

    bool IsSorted() const
    {
        return  true;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pointer vector set (size = " << size() << ") : ";

        return buffer.str();
    }


    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        std::copy(begin(), end(), std::ostream_iterator<TDataType>(rOStream, "\n "));
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@{

    /**
     * @class CompareKey
     * @brief A class providing comparison operators for keys in a custom map.
     * @details This class defines comparison operators for keys in a custom map. It allows you to compare keys using the specified comparison type and key extraction functions.
     */
    class CompareKey
    {
    public:
        /**
         * @brief Compare a key with a pointer to an object.
         * @details This function compares a key of type `key_type` with a pointer to an object of type `TPointerType` using the specified comparison function and key extraction function.
         * @param a The key of type `key_type`.
         * @param b The pointer to an object of type `TPointerType`.
         * @return True if the key `a` is less than the extracted key from `b`, false otherwise.
         */
        bool operator()(key_type a, TPointerType b) const
        {
            return TCompareType()(a, TGetKeyType()(*b));
        }

        /**
         * @brief Compare a pointer to an object with a key.
         * @details This function compares a pointer to an object of type `TPointerType` with a key of type `key_type` using the specified comparison function and key extraction function.
         * @param a The pointer to an object of type `TPointerType`.
         * @param b The key of type `key_type`.
         * @return True if the extracted key from `a` is less than the key `b`, false otherwise.
         */
        bool operator()(TPointerType a, key_type b) const
        {
            return TCompareType()(TGetKeyType()(*a), b);
        }

        /**
         * @brief Compare two pointers to objects.
         * @details This function compares two pointers to objects of type `TPointerType` using the specified comparison function and key extraction function.
         * @param a The pointer to the first object of type `TPointerType`.
         * @param b The pointer to the second object of type `TPointerType`.
         * @return True if the extracted key from `a` is less than the extracted key from `b`, false otherwise.
         */
        bool operator()(TPointerType a, TPointerType b) const
        {
            return TCompareType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };

    /**
    * @class EqualKeyTo
    * @brief A class providing equality comparison operators for keys in a custom map.
    * @details This class defines equality comparison operators for keys in a custom map. It allows you to check if a key is equal to the specified key using the specified equality function and key extraction function.
    */
    class EqualKeyTo
    {
        key_type mKey;  /// The key to compare against.
    public:
        /**
        * @brief Default constructor.
        * @details Initializes the `EqualKeyTo` object with a default-constructed key.
        */
        EqualKeyTo() : mKey() {}

        /**
        * @brief Explicit constructor with a specified key.
        * @details Initializes the `EqualKeyTo` object with the specified key.
        * @param Key The key of type `key_type` to compare against.
        */
        explicit EqualKeyTo(key_type Key) : mKey(Key) {}

        /**
        * @brief Compare a pointer to an object with the stored key.
        * @details This function checks if the key stored in this `EqualKeyTo` object is equal to the extracted key from the pointer to an object of type `TPointerType`.
        * @param a The pointer to an object of type `TPointerType`.
        * @return True if the stored key is equal to the extracted key from `a`, false otherwise.
        */
        bool operator()(TPointerType a) const
        {
            return TEqualType()(mKey, TGetKeyType()(*a));
        }

        /**
        * @brief Compare two pointers to objects with each other.
        * @details This function checks if the extracted keys from two pointers to objects of type `TPointerType` are equal using the specified equality function and key extraction function.
        * @param a The pointer to the first object of type `TPointerType`.
        * @param b The pointer to the second object of type `TPointerType`.
        * @return True if the extracted key from `a` is equal to the extracted key from `b`, false otherwise.
        */
        bool operator()(TPointerType a, TPointerType b) const
        {
            return TEqualType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };

    ///@}
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /// The data container holding the elements.
    TContainerType mData;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Extract the key from an iterator and apply a key extraction function.
     * @details This function extracts the key from an iterator and applies a key extraction function of type `TGetKeyType` to it.
     * @param i The iterator from which the key is extracted.
     * @return The key extracted from the iterator after applying the key extraction function.
     */
    key_type KeyOf(const_iterator i)
    {
        return TGetKeyType()(*i);
    }

    /**
     * @brief Extract the key from a pointer iterator and apply a key extraction function.
     * @details This function extracts the key from a pointer iterator and applies a key extraction function of type `TGetKeyType` to it.
     * @param i The pointer iterator from which the key is extracted.
     * @return The key extracted from the pointer iterator after applying the key extraction function.
     */
    key_type KeyOf(ptr_const_iterator i)
    {
        return TGetKeyType()(**i);
    }

    /**
    * @brief Extract the key from a data element and apply a key extraction function.
    * @details This function extracts the key from a data element of type `TDataType` and applies a key extraction function of type `TGetKeyType` to it.
    * @param i The data element from which the key is extracted.
    * @return The key extracted from the data element after applying the key extraction function.
    */
    key_type KeyOf(const TDataType &i)
    {
        return TGetKeyType()(i);
    }

    ///@}
    ///@name Serialization
    ///@{

    /**
     * @class Serializer
     * @brief A fried class responsible for handling the serialization process.
     */
    friend class Serializer;

    /**
     * @brief Extract the object's state and uses the Serializer to store it.
     * @param rSerializer Serializer instance to be used for saving.
     */
    virtual void save(Serializer& rSerializer) const
    {
        size_type local_size = mData.size();

        rSerializer.save("size", local_size);

        for(size_type i = 0 ; i < local_size ; i++)
            rSerializer.save("E", mData[i]);
    }

    /**
     * @brief Set the object's state based on data provided by the Serializer.
     * @param rSerializer Serializer instance to be used for loading.
     */
    virtual void load(Serializer& rSerializer)
    {
        size_type local_size;

        rSerializer.load("size", local_size);

        mData.resize(local_size);

        for(size_type i = 0 ; i < local_size ; i++)
            rSerializer.load("E", mData[i]);
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

}; // Class PointerVectorSet

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template<class TDataType,
         class TGetKeyType,
         class TCompareType,
         class TEqualType,
         class TPointerType,
         class TContainerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  PointerVectorSet<TDataType, TGetKeyType, TCompareType, TEqualType, TPointerType, TContainerType>& rThis);

/// output stream function
template<class TDataType,
         class TGetKeyType,
         class TCompareType,
         class TEqualType,
         class TPointerType,
         class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PointerVectorSet<TDataType, TGetKeyType, TCompareType, TEqualType, TPointerType, TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
