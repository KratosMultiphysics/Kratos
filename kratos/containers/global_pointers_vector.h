//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <vector>
#include <iostream>

// External includes
#include <boost/iterator/indirect_iterator.hpp>

// Project includes
#include "includes/define.h"
#include "includes/global_pointer.h"
#include "includes/serializer.h"

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
 * @class GlobalPointersVector
 * @ingroup KratosCore
 * @brief This class is a vector which stores global pointers
 * @details Uses boost::indirect_iterator
 * @tparam TDataType The type of data stored
 * @author Riccardo Rossi
 */
template< class TDataType >
class GlobalPointersVector final
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GlobalPointersVector
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointersVector);

    using TContainerType = std::vector<GlobalPointer<TDataType>>;
    using TPointerType = GlobalPointer<TDataType>;
    using data_type = TDataType;
    using value_type = TPointerType;
    using pointer = TPointerType;
    using const_pointer = const TPointerType;
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
    GlobalPointersVector() : mData() {}

    GlobalPointersVector(const std::initializer_list<GlobalPointer<TDataType>>& l)
        : mData(l)
    {}

    /// Destructor.
    ~GlobalPointersVector() {}

    /**
    * @brief Fill the container from another container.
    * @tparam TContainerType The type of the source container.
    * @param rContainer The source container.
    */
    template <class TContainerType>
    void FillFromContainer(TContainerType& rContainer)
    {
        this->reserve(rContainer.size());
        for (auto& r_item : rContainer) {
            this->push_back(GlobalPointer<TDataType>(&r_item));
        }
    }

    /**
     * @brief Sort the elements in the container.
     */
    void Sort()
    {
        std::sort(mData.begin(), mData.end(), GlobalPointerCompare<TDataType>());
    }

    /**
     * @brief Remove duplicate elements from the container.
     * @details This function first sorts the elements and then removes duplicates.
     */
    void Unique()
    {
        Sort();
        auto end_it = mData.end();
        auto new_end_it = std::unique(mData.begin(), mData.end(), GlobalPointerComparor<TDataType>());
        this->erase(new_end_it, end_it);
    }

    /**
    * @brief Reduce the capacity of the container to fit its size.
    */
    void shrink_to_fit()
    {
        mData.shrink_to_fit();
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Assignment operator to copy the contents of another GlobalPointersVector.
     * @param rOther The GlobalPointersVector to copy from.
     * @return A reference to the modified GlobalPointersVector.
     */
    GlobalPointersVector& operator=(const GlobalPointersVector& rOther)
    {
        mData = rOther.mData;
        return *this;
    }

    /**
     * @brief Access an element in the container by index.
     * @param i The index of the element to access.
     * @return A reference to the element.
     */
    TDataType& operator[](const size_type& i)
    {
        return *(mData[i]);
    }

    /**
     * @brief Access a constant element in the container by index.
     * @param i The index of the element to access.
     * @return A constant reference to the element.
     */
    TDataType const& operator[](const size_type& i) const
    {
        return *(mData[i]);
    }

    /**
     * @brief Access an element in the container by index.
     * @param i The index of the element to access.
     * @return A reference to the element.
     */
    pointer& operator()(const size_type& i)
    {
        return mData[i];
    }

    /**
     * @brief Access a constant element in the container by index.
     * @param i The index of the element to access.
     * @return A constant reference to the element.
     */
    const_pointer& operator()(const size_type& i) const
    {
        return mData[i];
    }

    /**
     * @brief Equality comparison operator to check if two GlobalPointersVector objects are equal.
     * @details This function checks if the sizes are equal and then compares the elements for equality
     * using the EqualKeyTo() function.
     * @param r The GlobalPointersVector to compare with.
     * @return True if the containers are equal, false otherwise.
     */
    bool operator==(const GlobalPointersVector& r) const // nothrow
    {
        if (size() != r.size())
            return false;
        else
            return std::equal(mData.begin(), mData.end(), r.mData.begin(), this->EqualKeyTo());
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
    iterator                   end()
    {
        return iterator( mData.end() );
    }
    const_iterator             end() const
    {
        return const_iterator( mData.end() );
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
        assert( !mData.empty() );
        return *(mData.front());
    }
    const_reference  front() const /* nothrow */
    {
        assert( !mData.empty() );
        return *(mData.front());
    }
    reference        back()        /* nothrow */
    {
        assert( !mData.empty() );
        return *(mData.back());
    }
    const_reference  back() const  /* nothrow */
    {
        assert( !mData.empty() );
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

    void swap(GlobalPointersVector& rOther)
    {
        mData.swap(rOther.mData);
    }

    void push_back(TPointerType x)
    {
        mData.push_back(x);
    }

    // template<class TOtherDataType>
    // void push_back(TOtherDataType const& x)
    // {
    //     push_back(TPointerType(new TOtherDataType(x)));
    // }

    // template<class TOtherDataType>
    // iterator insert(iterator Position, const TOtherDataType& rData)
    // {
    //     return iterator(mData.insert(Position, TPointerType(new TOtherDataType(rData))));
    // }

    iterator insert(iterator Position, const TPointerType pData)
    {
        return iterator(mData.insert(Position, pData));
    }

    template <class InputIterator>
    void insert(InputIterator First, InputIterator Last)
    {
        for(; First != Last; ++First)
            insert(*First);
    }

    iterator erase(iterator pos)
    {
        return iterator(mData.erase(pos.base()));
    }

    iterator erase( iterator first, iterator last )
    {
        return iterator( mData.erase( first.base(), last.base() ) );
    }

    void clear()
    {
        mData.clear();
    }

    void resize(size_type new_dim) const
    {
        return mData.resize(new_dim);
    }
    void resize(size_type new_dim)
    {
        mData.resize(new_dim);
    }

    void reserve(int dim)
    {
        mData.reserve(dim);
    }

    int capacity()
    {
        return mData.capacity();
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

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointersVector" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointersVector";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

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

    TContainerType mData; /// The data stored in the class

    ///@}
    ///@name Private Operators
    ///@{

    friend class Serializer;

    /**
    * @brief Serialize the GlobalPointersVector for saving.
    * @details This function saves the size of the container and each element in the container using the provided serializer.
    * @param rSerializer The serializer used for saving.
    */
    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Size", this->size());

        for (std::size_t i = 0; i < this->size(); i++) {
            rSerializer.save("Data", mData[i]);
        }
    }

    /**
    * @brief Deserialize the GlobalPointersVector for loading.
    * @details This function loads the size of the container and each element from the serializer and populates the container.
    * @param rSerializer The serializer used for loading.
    */
    void load(Serializer& rSerializer)
    {
        std::size_t size;

        rSerializer.load("Size", size);

        for (std::size_t i = 0; i < size; i++) {
            GlobalPointer<TDataType> p(nullptr);
            rSerializer.load("Data", p);
            this->push_back(p);
        }
    }

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

    ///@}
}; // Class GlobalPointersVector

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template< class TDataType >
inline std::istream& operator >> (std::istream& rIStream,
                                  GlobalPointersVector<TDataType>& rThis)
{
    return rIStream;
}

/// output stream function
template< class TDataType >
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GlobalPointersVector<TDataType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup KratosCore

}  // namespace Kratos.
