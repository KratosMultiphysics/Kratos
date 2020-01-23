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


#if !defined(KRATOS_POINTER_VECTOR_H_INCLUDED )
#define  KRATOS_POINTER_VECTOR_H_INCLUDED



// System includes
#include <functional>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// External includes
#include <boost/iterator/indirect_iterator.hpp>


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"


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

/// PointerVector is a  container like stl vector but using a vector to store pointers to its data.
/** PointerVector is a container like stl vector
    but using a vector to store pointers its data. Many methods are
    copied from boost ptr_container library. There is modification
    to make it capable to work with shared pointers.

    This Container unlike the boost one does not free the memory by
    itself and relies on using of counted pointers or manual
    deleting.
 */
template<class TDataType,
         class TPointerType = typename TDataType::Pointer,
         class TContainerType = std::vector<TPointerType> >
class PointerVector
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointerVector
    KRATOS_CLASS_POINTER_DEFINITION(PointerVector);

    /// data type stores in this container.
    typedef TDataType data_type;
    typedef TDataType value_type;
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
    PointerVector() : mData() {}

    template <class TInputIteratorType>
    PointerVector(TInputIteratorType First, TInputIteratorType Last)
        : mData(First, Last)
    {
    }

    PointerVector(const PointerVector& rOther) :  mData(rOther.mData) {}

    explicit PointerVector(const TContainerType& rContainer) :  mData(rContainer)
    {
    }

    explicit PointerVector(std::size_t NewSize) :  mData(NewSize)
    {
    }
/*
    template<class TOtherDataType>
    PointerVector(std::size_t NewSize, TOtherDataType const& Value) :  mData(NewSize)
    {
        for(size_type i = 0 ; i < NewSize ; i++)
            mData[i] = pointer(new TOtherDataType(Value));
    }*/

    /// Destructor.
    virtual ~PointerVector() {}


    ///@}
    ///@name Operators
    ///@{

    PointerVector& operator=(const PointerVector& rOther)
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

    bool operator==( const PointerVector& r ) const // nothrow
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

    void swap(PointerVector& rOther)
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
        buffer << "PointerVector (size = " << size() << ") : ";

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        size_type local_size = mData.size();

        rSerializer.save("size", local_size);

        for(size_type i = 0 ; i < local_size ; i++)
            rSerializer.save("E", mData[i]);
    }

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

}; // Class PointerVector

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
                                  PointerVector<TDataType, TPointerType, TContainerType>& rThis);

/// output stream function
template<class TDataType,
         class TPointerType,
         class TContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PointerVector<TDataType, TPointerType, TContainerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_POINTER_VECTOR_SET_H_INCLUDED  defined 
