//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#if !defined(KRATOS_WEAK_POINTER_VECTOR_H_INCLUDED )
#define  KRATOS_WEAK_POINTER_VECTOR_H_INCLUDED



// System includes
#include <functional>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// External includes


// Project includes
#include "includes/define.h"
#include "containers/weak_pointer_vector_iterator.h"

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

/// WeakPointerVector is a  container like stl vector but using a vector to store pointers to its data.
/** WeakPointerVector is a container like stl vector
    but using a vector to store pointers its data. Many methods are
    copied from boost ptr_container library. There is modification
    to make it capable to work with shared pointers.

    This Container unlike the boost one does not free the memory by
    itself and relies on using of counted pointers or manual
    deleting.
 */
template<class TDataType,
         class TPointerType = Kratos::weak_ptr<TDataType>,
         class TContainerType = std::vector<TPointerType> >
class WeakPointerVector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WeakPointerVector
    KRATOS_CLASS_POINTER_DEFINITION(WeakPointerVector);

    /// data type stores in this container.
    typedef TDataType data_type;
    typedef TPointerType value_type;
    typedef TPointerType pointer;
    typedef const TPointerType const_pointer;
    typedef TDataType& reference;
    typedef const TDataType& const_reference;
    typedef TContainerType ContainerType;

    typedef WeakPointerVectorIterator<typename TContainerType::iterator, TDataType>                iterator;
    typedef WeakPointerVectorIterator<typename TContainerType::const_iterator, TDataType>          const_iterator;
    typedef WeakPointerVectorIterator<typename TContainerType::reverse_iterator, TDataType>        reverse_iterator;
    typedef WeakPointerVectorIterator<typename TContainerType::const_reverse_iterator, TDataType>  const_reverse_iterator;

    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::iterator ptr_iterator;
    typedef typename TContainerType::const_iterator ptr_const_iterator;
    typedef typename TContainerType::reverse_iterator ptr_reverse_iterator;
    typedef typename TContainerType::const_reverse_iterator ptr_const_reverse_iterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WeakPointerVector() : mData() {}

    template <class TInputIteratorType>
    WeakPointerVector(TInputIteratorType First, TInputIteratorType Last)
        : mData(First, Last)
    {
    }

    WeakPointerVector(const WeakPointerVector& rOther) :  mData(rOther.mData) {}

    WeakPointerVector(const TContainerType& rContainer) :  mData(rContainer)
    {
    }

    /// Destructor.
    virtual ~WeakPointerVector() {}


    ///@}
    ///@name Operators
    ///@{

    WeakPointerVector& operator=(const WeakPointerVector& rOther)
    {
        mData = rOther.mData;
        return *this;
    }

    TDataType& operator[](const size_type& i)
    {
        return *(mData[i].lock());
    }

    TDataType const& operator[](const size_type& i) const
    {
        return *(mData[i].lock());
    }

    pointer& operator()(const size_type& i)
    {
        return mData[i];
    }

    const_pointer& operator()(const size_type& i) const
    {
        return mData[i];
    }

    bool operator==( const WeakPointerVector& r ) const // nothrow
    {
        if( size() != r.size() )
            return false;
        else
            return std::equal(mData.begin(), mData.end(), r.mData.begin());
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
        assert( !empty() );
        return *((mData.front()).lock());
    }
    const_reference  front() const /* nothrow */
    {
        assert( !empty() );
        return *((mData.front()).lock());
    }
    reference        back()        /* nothrow */
    {
        assert( !empty() );
        return *((mData.back()).lock());
    }
    const_reference  back() const  /* nothrow */
    {
        assert( !empty() );
        return *((mData.back()).lock());
    }

    size_type size() const
    {
        return mData.size();
    }
    void resize(size_type new_dim) const
    {
        return mData.resize(new_dim);
    }
    void resize(size_type new_dim)
    {
        mData.resize(new_dim);
    }

    size_type max_size() const
    {
        return mData.max_size();
    }

    void swap(WeakPointerVector& rOther)
    {
        mData.swap(rOther.mData);
    }

    void push_back(TPointerType x)
    {
        mData.push_back(x);
    }

    void push_back(const_reference x)
    {
        push_back(TPointerType(new TDataType(x)));
    }

    iterator insert(iterator Position, const TDataType& rData)
    {
        return iterator(mData.insert(Position, TPointerType(new TDataType(rData))));
    }

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
    ///@name Inquiry
    ///@{

    bool empty() const
    {
        return mData.empty();
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "WeakPointerVector (size = " << size() << ") : ";

        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        std::copy(begin(), end(), std::ostream_iterator<TDataType>(rOStream, "\t "));
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

}; // Class WeakPointerVector

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TDataType,
         class TPointerType,
         class TContainerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  WeakPointerVector<TDataType, TPointerType, TContainerType>& rThis);

/// output stream function
template<class TDataType,
         class TPointerType,
         class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const WeakPointerVector<TDataType, TPointerType, TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_WEAK_POINTER_VECTOR_H_INCLUDED  defined
