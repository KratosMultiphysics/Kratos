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
    typedef typename std::remove_reference<decltype(std::declval<TGetKeyType>()(std::declval<TDataType>()))>::type key_type;

    // Data type stored in this container.
    using data_type = TDataType;
    using value_type = TDataType;
    using key_compare = TCompareType;
    using pointer = TPointerType;
    using reference = TDataType&;
    using const_reference = const TDataType&;
    using ContainerType = TContainerType;

    using iterator = boost::indirect_iterator<typename TContainerType::iterator>;
    using const_iterator = boost::indirect_iterator<typename TContainerType::const_iterator>;
    using reverse_iterator = boost::indirect_iterator<typename TContainerType::reverse_iterator>;
    using const_reverse_iterator = boost::indirect_iterator<typename TContainerType::const_reverse_iterator>;

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
    PointerVectorSet() : mData(), mSortedPartSize(size_type()), mMaxBufferSize(1) {}

    /**
    * @brief Constructs a PointerVectorSet from a range of elements.
    * @details This constructor initializes a PointerVectorSet with elements in the range [First, Last).
    * @tparam TInputIteratorType The type of the input iterator.
    * @param First An input iterator pointing to the beginning of the range.
    * @param Last An input iterator pointing to the end of the range.
    * @param NewMaxBufferSize The maximum buffer size (default is 1).
    */
    template <class TInputIteratorType>
    PointerVectorSet(TInputIteratorType First, TInputIteratorType Last, size_type NewMaxBufferSize = 1)
    : mSortedPartSize(size_type()), mMaxBufferSize(NewMaxBufferSize)
    {
    for (; First != Last; ++First)
        insert(begin(), *First);
    }

    /**
     * @brief Copy constructor for PointerVectorSet.
     * @param rOther The PointerVectorSet to copy from.
     */
    PointerVectorSet(const PointerVectorSet& rOther)
        :  mData(rOther.mData), mSortedPartSize(rOther.mSortedPartSize), mMaxBufferSize(rOther.mMaxBufferSize) {}
    
    /**
     * @brief Constructs a PointerVectorSet from a container.
     * @details This constructor initializes a PointerVectorSet with elements from a container.
     * @param rContainer The container to copy elements from.
     */
    explicit PointerVectorSet(const TContainerType& rContainer) :  mData(rContainer), mSortedPartSize(size_type()), mMaxBufferSize(1)
    {
        Sort();
        std::unique(mData.begin(), mData.end(), EqualKeyTo());
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
        mSortedPartSize = rOther.mSortedPartSize;

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
        ptr_iterator sorted_part_end;

        if (mData.size() - mSortedPartSize >= mMaxBufferSize) {
            Sort();
            sorted_part_end = mData.end();
        } else {
            sorted_part_end = mData.begin() + mSortedPartSize;
        }

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end) {
            mSortedPartSize++;
            return **mData.insert(sorted_part_end, TPointerType(new TDataType(Key)));
        }

        if (!EqualKeyTo(Key)(*i)) {
            if ((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end()) {
                mData.push_back(TPointerType(new TDataType(Key)));
                return **(mData.end() - 1);
            }
        }

        return **i;
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
        ptr_iterator sorted_part_end;

        if (mData.size() - mSortedPartSize >= mMaxBufferSize) {
            Sort();
            sorted_part_end = mData.end();
        } else
            sorted_part_end = mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end) {
            mSortedPartSize++;
            return *mData.insert(sorted_part_end, TPointerType(new TDataType(Key)));
        }

        if (!EqualKeyTo(Key)(*i))
            if ((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end()) {
                mData.push_back(TPointerType(new TDataType(Key)));
                return *(mData.end() - 1);
            }

        return *i;
    }

    /**
     * @brief Equality comparison operator for two PointerVectorSet objects.
     * @details This operator checks if two PointerVectorSet objects are equal by comparing their sizes
     * and the equality of their elements using the EqualKeyTo comparison function.
     * @param r The PointerVectorSet to compare with.
     * @return true if the two sets are equal, false otherwise.
     */
    bool operator==(const PointerVectorSet& r) const // nothrow
    {
        if (size() != r.size())
            return false;
        else
            return std::equal(mData.begin(), mData.end(), r.mData.begin(), EqualKeyTo());
    }

    /**
     * @brief Less than comparison operator for two PointerVectorSet objects.
     * @details This operator checks if one PointerVectorSet is less than another by comparing their
     * elements using the CompareKey comparison function in a lexicographical order.
     * @param r The PointerVectorSet to compare with.
     * @return true if this set is less than r, false otherwise.
     */
    bool operator<(const PointerVectorSet& r) const // nothrow
    {
        return std::lexicographical_compare(mData.begin(), mData.end(), r.mData.begin(), r.mData.end(), CompareKey());
    }

    ///@}
    ///@name Operations
    ///@{

    iterator                   begin()
    {
        return iterator( mData.begin() );
    }
    const_iterator             begin() const
    {
        return const_iterator( mData.begin() );
    }

    const_iterator cbegin()
    {
        return const_iterator(mData.begin());
    }

    const_iterator cbegin() const
    {
        return const_iterator(mData.begin());
    }

    iterator                   end()
    {
        return iterator( mData.end() );
    }
    const_iterator             end() const
    {
        return const_iterator( mData.end() );
    }

    const_iterator cend()
    {
        return const_iterator(mData.end());
    }

    const_iterator cend() const
    {
        return const_iterator(mData.end());
    }

    reverse_iterator           rbegin()
    {
        return reverse_iterator( mData.rbegin() );
    }
    const_reverse_iterator     rbegin() const
    {
        return const_reverse_iterator( mData.rbegin() );
    }
    reverse_iterator           rend()
    {
        return reverse_iterator( mData.rend() );
    }
    const_reverse_iterator     rend() const
    {
        return const_reverse_iterator( mData.rend() );
    }
    ptr_iterator               ptr_begin()
    {
        return mData.begin();
    }
    ptr_const_iterator         ptr_begin() const
    {
        return mData.begin();
    }
    ptr_iterator               ptr_end()
    {
        return mData.end();
    }
    ptr_const_iterator         ptr_end() const
    {
        return mData.end();
    }
    ptr_reverse_iterator       ptr_rbegin()
    {
        return mData.rbegin();
    }
    ptr_const_reverse_iterator ptr_rbegin() const
    {
        return mData.rbegin();
    }
    ptr_reverse_iterator       ptr_rend()
    {
        return mData.rend();
    }
    ptr_const_reverse_iterator ptr_rend() const
    {
        return mData.rend();
    }

    reference        front()       /* nothrow */
    {
        assert( !empty() );
        return *(mData.front());
    }
    const_reference  front() const /* nothrow */
    {
        assert( !empty() );
        return *(mData.front());
    }
    reference        back()        /* nothrow */
    {
        assert( !empty() );
        return *(mData.back());
    }
    const_reference  back() const  /* nothrow */
    {
        assert( !empty() );
        return *(mData.back());
    }

    size_type size() const
    {
        return mData.size();
    }

    size_type max_size() const
    {
        return mData.max_size();
    }

    key_compare key_comp() const
    {
        return TCompareType();
    }

    /**
     * @brief Swaps the contents of this PointerVectorSet with another.
     * @details This function swaps the contents of this PointerVectorSet with another set, including
     * mSortedPartSize, mMaxBufferSize, and mData.
     * @param rOther The other PointerVectorSet to swap with.
     */
    void swap(PointerVectorSet& rOther)
    {
        std::swap(mSortedPartSize, rOther.mSortedPartSize);
        std::swap(mMaxBufferSize, rOther.mMaxBufferSize);
        mData.swap(rOther.mData);
    }

    /**
     * @brief Adds a pointer to the end of the set.
     * @details This function appends a given pointer to the end of the set.
     * @param x The pointer to be added to the end of the set.
     */
    void push_back(TPointerType x)
    {
        mData.push_back(x);
    }

    /**
     * @brief Removes the last element from the set.
     * @details This function removes the last element (pointer) from the set and updates mSortedPartSize
     * if necessary.
     */
    void pop_back()
    {
        mData.pop_back();
        if (mSortedPartSize > mData.size())
            mSortedPartSize = mData.size();
    }

    /**
     * @brief Inserts a pointer at the specified position.
     * @details This function inserts a given pointer at the specified position in the set. It also maintains
     * the sorting order and updates mSortedPartSize if necessary.
     * @param Position An iterator pointing to the position where the pointer should be inserted.
     * @param pData The pointer to be inserted.
     * @return An iterator pointing to the inserted element.
     */
    iterator insert(iterator Position, const TPointerType pData)
    {
        ptr_iterator sorted_part_end;

        key_type key = KeyOf(*pData);

        if (mData.size() - mSortedPartSize >= mMaxBufferSize) {
            Sort();
            sorted_part_end = mData.end();
        } else
            sorted_part_end = mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, key, CompareKey()));
        if (i == sorted_part_end) {
            mSortedPartSize++;
            return mData.insert(sorted_part_end, pData);
        }

        if (!EqualKeyTo(key)(*i))
            if ((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(key))) == mData.end()) {
                mData.push_back(pData);
                return iterator(mData.end() - 1);
            }

        *i = pData;
        return i;
    }

    /**
     * @brief Insert elements from a range of iterators.
     * @details This function inserts elements from a range defined by the iterators `First` and `Last`
     * into the set. It uses the `insert` function to insert each element.
     * @param First An input iterator pointing to the beginning of the range to insert.
     * @param Last An input iterator pointing to the end of the range to insert.
     */
    template <class InputIterator>
    void insert(InputIterator First, InputIterator Last)
    {
        for (; First != Last; ++First)
            insert(begin(), *First);
    }

    /**
     * @brief Erase an element at the specified position.
     * @details This function erases the element at the specified position and updates `mSortedPartSize`
     * to match the size of the data container. If the provided position is equal to `end()`,
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
        mSortedPartSize = mData.size();
        return new_end;
    }

    /**
     * @brief Erase a range of elements defined by iterators.
     * @details This function erases a range of elements defined by the iterators `first` and `last`
     * and updates `mSortedPartSize` to match the size of the data container.
     * @param first An iterator pointing to the beginning of the range to erase.
     * @param last An iterator pointing to the end of the range to erase.
     * @return An iterator pointing to the element following the last erased element.
     */
    iterator erase(iterator first, iterator last)
    {
        iterator new_end = iterator(mData.erase(first.base(), last.base()));
        // TODO: Sorted part size must change
        mSortedPartSize = mData.size();
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
     * @details This function clears the set by removing all elements, resetting `mSortedPartSize` to zero,
     * and setting `mMaxBufferSize` to 1.
     */
    void clear()
    {
        mData.clear();
        mSortedPartSize = size_type();
        mMaxBufferSize = 1;
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
        ptr_iterator sorted_part_end;

        if (mData.size() - mSortedPartSize >= mMaxBufferSize) {
            Sort();
            sorted_part_end = mData.end();
        } else
            sorted_part_end = mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end || (!EqualKeyTo(Key)(*i)))
            if ((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
                return mData.end();

        return i;
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
        ptr_const_iterator sorted_part_end(mData.begin() + mSortedPartSize);

        ptr_const_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end || (!EqualKeyTo(Key)(*i)))
            if ((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
                return mData.end();

        return const_iterator(i);
    }

    /**
     * @brief Count the number of elements with the specified key.
     * @details This function counts the number of elements with the specified key in the set. It returns 1
     * if the element is found and 0 if it's not found.
     * @param Key The key to count.
     * @return The number of elements with the specified key (0 or 1).
     */
    size_type count(const key_type& Key)
    {
        return find(Key) == mData.end() ? 0 : 1;
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

    /**
     * @brief Sort the elements in the set.
     * @details This function sorts the elements in the set using the CompareKey comparison function. After sorting,
     * it updates mSortedPartSize to match the size of the data container.
     */
    void Sort()
    {
        std::sort(mData.begin(), mData.end(), CompareKey());
        mSortedPartSize = mData.size();
    }

    /**
     * @brief Remove duplicate elements from the set.
     * @details This function removes duplicate elements from the set using the EqualKeyTo comparison function. After
     * removing duplicates, it updates mSortedPartSize to match the size of the data container.
     */
    void Unique()
    {
        typename TContainerType::iterator end_it = mData.end();
        std::sort(mData.begin(), mData.end(), CompareKey());
        typename TContainerType::iterator new_end_it = std::unique(mData.begin(), mData.end(), EqualKeyTo());
        mData.erase(new_end_it, end_it);
        mSortedPartSize = mData.size();
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
    size_type GetMaxBufferSize() const 
    {
        return mMaxBufferSize;
    }

    /** 
     * @brief Set the maximum size of buffer used in the container.
     * @details This container uses a buffer which keep data unsorted. After buffer size arrived to the MaxBufferSize it will sort all container and empties buffer.
     * @param NewSize Is the new buffer maximum size. 
     */
    void SetMaxBufferSize(const size_type NewSize)
    {
        mMaxBufferSize = NewSize;
    }

    /**
     * @brief Get the sorted part size of buffer used in the container. 
     */
    size_type GetSortedPartSize() const 
    {
        return mSortedPartSize;
    }

    /** 
     * @brief Set the sorted part size of buffer used in the container.
     * @param NewSize Is the new buffer maximum size. 
     */
    void SetSortedPartSize(const size_type NewSize)
    {
        mSortedPartSize = NewSize;
    }

    ///@}
    ///@name Inquiry
    ///@{

    bool empty() const
    {
        return mData.empty();
    }

    bool IsSorted() const
    {
        return (mSortedPartSize == mData.size());
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

    class CompareKey
    {
    public:
        bool operator()(key_type a, TPointerType b) const
        {
            return TCompareType()(a, TGetKeyType()(*b));
        }
        bool operator()(TPointerType a, key_type b) const
        {
            return TCompareType()(TGetKeyType()(*a), b);
        }
        bool operator()(TPointerType a, TPointerType b) const
        {
            return TCompareType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };

    class EqualKeyTo
    {
        key_type mKey;
    public:
        EqualKeyTo() : mKey() {}
        explicit EqualKeyTo(key_type Key) : mKey(Key) {}
        bool operator()(TPointerType a) const
        {
            return TEqualType()(mKey, TGetKeyType()(*a));
        }
        bool operator()(TPointerType a, TPointerType b) const
        {
            return TEqualType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    TContainerType mData;
    size_type  mSortedPartSize;
    size_type mMaxBufferSize;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    key_type KeyOf(iterator i)
    {
        return TGetKeyType()(*i);
    }

    key_type KeyOf(ptr_iterator i)
    {
        return TGetKeyType()(**i);
    }

    key_type KeyOf(const TDataType & i)
    {
        return TGetKeyType()(i);
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        size_type local_size = mData.size();

        rSerializer.save("size", local_size);

        for(size_type i = 0 ; i < local_size ; i++)
            rSerializer.save("E", mData[i]);

        rSerializer.save("Sorted Part Size",mSortedPartSize);
        rSerializer.save("Max Buffer Size",mMaxBufferSize);
    }

    virtual void load(Serializer& rSerializer)
    {
        size_type local_size;

        rSerializer.load("size", local_size);

        mData.resize(local_size);

        for(size_type i = 0 ; i < local_size ; i++)
            rSerializer.load("E", mData[i]);

        rSerializer.load("Sorted Part Size",mSortedPartSize);
        rSerializer.load("Max Buffer Size",mMaxBufferSize);
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
