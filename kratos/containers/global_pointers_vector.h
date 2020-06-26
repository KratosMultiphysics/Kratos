//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED )
#define  KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED

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
class GlobalPointersVector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GlobalPointersVector
    KRATOS_CLASS_POINTER_DEFINITION(GlobalPointersVector);

    typedef std::vector< GlobalPointer<TDataType> > TContainerType;
    typedef GlobalPointer<TDataType> TPointerType;
    typedef TDataType data_type;
    typedef TPointerType value_type;
    typedef TPointerType pointer;
    typedef const TPointerType const_pointer;
    typedef TDataType& reference;
    typedef const TDataType& const_reference;
    typedef TContainerType ContainerType;


    typedef boost::indirect_iterator<typename TContainerType::iterator>                iterator;
    typedef boost::indirect_iterator<typename TContainerType::const_iterator>          const_iterator;
    typedef boost::indirect_iterator<typename TContainerType::reverse_iterator>        reverse_iterator;
    typedef boost::indirect_iterator<typename TContainerType::const_reverse_iterator>  const_reverse_iterator;

    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::iterator ptr_iterator;
    typedef typename TContainerType::const_iterator ptr_const_iterator;
    typedef typename TContainerType::reverse_iterator ptr_reverse_iterator;
    typedef typename TContainerType::const_reverse_iterator ptr_const_reverse_iterator;
    typedef typename TContainerType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GlobalPointersVector() : mData() {}

    GlobalPointersVector(const std::initializer_list<GlobalPointer<TDataType>>& l)
        : mData(l)
    {}

    /// Destructor.
    virtual ~GlobalPointersVector() {}

    template < class TContainerType >
    void FillFromContainer( TContainerType& container)
    {
        this->reserve(container.size());
        for(auto& item : container)
        {
            this->push_back(GlobalPointer<TDataType>(&item));
        }
    }

    void Sort()
    {
        std::sort(mData.begin(), mData.end(), GlobalPointerCompare<TDataType>());
    }

    void Unique()
    {
        Sort();
        auto end_it = mData.end();
        auto new_end_it = std::unique(mData.begin(), mData.end(), GlobalPointerComparor<TDataType>());
        this->erase(new_end_it, end_it);
    }

    void shrink_to_fit()
    {
        mData.shrink_to_fit();
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "GlobalPointersVector" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "GlobalPointersVector";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    GlobalPointersVector& operator=(const GlobalPointersVector& rOther)
    {
        mData = rOther.mData;
        return *this;
    }

    TDataType& operator[](const size_type& i)
    {
        return *(mData[i]);
    }

    TDataType const& operator[](const size_type& i) const
    {
        return *(mData[i]);
    }

    pointer& operator()(const size_type& i)
    {
        return mData[i];
    }

    const_pointer& operator()(const size_type& i) const
    {
        return mData[i];
    }

    bool operator==( const GlobalPointersVector& r ) const // nothrow
    {
        if( size() != r.size() )
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

//     template<class TOtherDataType>
//     void push_back(TOtherDataType const& x)
//     {
//         push_back(TPointerType(new TOtherDataType(x)));
//     }
/*
    template<class TOtherDataType>
    iterator insert(iterator Position, const TOtherDataType& rData)
    {
        return iterator(mData.insert(Position, TPointerType(new TOtherDataType(rData))));
    }*/

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
    TContainerType mData;


    ///@}
    ///@name Private Operators
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Size", this->size());

        for(std::size_t i=0; i<this->size(); i++) {
            rSerializer.save("Data", mData[i]);
        }
    }

    void load(Serializer& rSerializer)
    {
        std::size_t size;

        rSerializer.load("Size", size);

        for(std::size_t i = 0; i < size; i++) {
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

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GLOBAL_POINTER_VECTOR_H_INCLUDED  defined
