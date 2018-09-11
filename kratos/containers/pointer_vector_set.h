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
//                    
//


#if !defined(KRATOS_POINTER_VECTOR_SET_H_INCLUDED )
#define  KRATOS_POINTER_VECTOR_SET_H_INCLUDED



// System includes
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes
#include "boost/smart_ptr.hpp"
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

/// PointerVectorSet is a sorted associative container like stl set but using a vector to store pointers to its data.
/** PointerVectorSet is a sorted associative container like stl set
    but using a vector to store pointers its data. Many methods are
    copied from boost ptr_container library. There is modification
    to make it capable to work with shared pointers.

    This Container unlike the boost one does not free the memory by
    itself and relies on using of counted pointers or manual
    deleting.
 */
template<class TDataType,
         class TGetKeyType = SetIdentityFunction<TDataType>,
         class TCompareType = std::less<typename TGetKeyType::result_type>,
         class TEqualType = std::equal_to<typename TGetKeyType::result_type>,
         class TPointerType = Kratos::shared_ptr<TDataType>,
         class TContainerType = std::vector<TPointerType> >
class PointerVectorSet
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointerVectorSet
    KRATOS_CLASS_POINTER_DEFINITION(PointerVectorSet);

    /// Key type for searching in this container.
    typedef typename TGetKeyType::result_type key_type;

    /// data type stores in this container.
    typedef TDataType data_type;
    typedef TDataType value_type;
    typedef TCompareType key_compare;
    typedef TPointerType pointer;
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
    PointerVectorSet() : mData(), mSortedPartSize(size_type()), mMaxBufferSize(1) {}

    template <class TInputIteratorType>
    PointerVectorSet(TInputIteratorType First, TInputIteratorType Last, size_type NewMaxBufferSize = 1)
        : mSortedPartSize(size_type()), mMaxBufferSize(NewMaxBufferSize)
    {
        for(; First != Last; ++First)
            insert(begin(), *First);
    }

    PointerVectorSet(const PointerVectorSet& rOther)
        :  mData(rOther.mData), mSortedPartSize(rOther.mSortedPartSize), mMaxBufferSize(rOther.mMaxBufferSize) {}

    explicit PointerVectorSet(const TContainerType& rContainer) :  mData(rContainer), mSortedPartSize(size_type()), mMaxBufferSize(1)
    {
        Sort();
        std::unique(mData.begin(), mData.end(), EqualKeyTo());
    }

    /// Destructor.
    virtual ~PointerVectorSet() {}


    ///@}
    ///@name Operators
    ///@{

    PointerVectorSet& operator=(const PointerVectorSet& rOther)
    {
        mData = rOther.mData;
        mSortedPartSize = rOther.mSortedPartSize;

        return *this;
    }

    TDataType& operator[](const key_type& Key)
    {
        ptr_iterator sorted_part_end;

        if(mData.size() - mSortedPartSize >= mMaxBufferSize)
        {
            Sort();
            sorted_part_end = mData.end();
        }
        else
            sorted_part_end	= mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end)
        {
            mSortedPartSize++;
            return **mData.insert(sorted_part_end, TPointerType(new TDataType(Key)));
        }

        if (!EqualKeyTo(Key)(*i))
            if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
            {
                mData.push_back(TPointerType(new TDataType(Key)));
                return **(mData.end()-1);
// 	  return **(--mData.end());
            }

        return **i;
    }

    pointer& operator()(const key_type& Key)
    {
        ptr_iterator sorted_part_end;

        if(mData.size() - mSortedPartSize >= mMaxBufferSize)
        {
            Sort();
            sorted_part_end = mData.end();
        }
        else
            sorted_part_end	= mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end)
        {
            mSortedPartSize++;
            return *mData.insert(sorted_part_end, TPointerType(new TDataType(Key)));
        }

        if (!EqualKeyTo(Key)(*i))
            if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
            {
                mData.push_back(TPointerType(new TDataType(Key)));
                return *(mData.end()-1);
// 	  return *(--mData.end());
            }

        return *i;
    }

    bool operator==( const PointerVectorSet& r ) const // nothrow
    {
        if( size() != r.size() )
            return false;
        else
            return std::equal(mData.begin(), mData.end(), r.mData.begin(), EqualKeyTo());
    }

    bool operator<( const PointerVectorSet& r ) const // nothrow
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

    key_compare key_comp() const
    {
        return TCompareType();
    }

    void swap(PointerVectorSet& rOther)
    {
        mData.swap(rOther.mData);
    }

    void push_back(TPointerType x)
    {
        mData.push_back(x);
    }
    void pop_back()
    {
//KRATOS_WATCH("before mData.pop_back")
//KRATOS_WATCH(**((mData.end()-1).base()))
        mData.pop_back();
//KRATOS_WATCH("after mData.pop_back")
        if(mSortedPartSize>mData.size())
            mSortedPartSize = mData.size();
//KRATOS_WATCH("finished pop_back")
    }

//     template<class TOtherDataType>
//     void push_back(TOtherDataType const& x)
//     {
//         push_back(TPointerType(new TOtherDataType(x)));
//     }

//     template<class TOtherDataType>
//     iterator insert(iterator Position, const TOtherDataType& rData)
//     {
//         ptr_iterator sorted_part_end;
// 
//         key_type key = KeyOf(rData);
// 
//         if(mData.size() - mSortedPartSize >= mMaxBufferSize)
//         {
//             Sort();
//             sorted_part_end = mData.end();
//         }
//         else
//             sorted_part_end	= mData.begin() + mSortedPartSize;
// 
//         ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, key, CompareKey()));
//         if (i == sorted_part_end)
//         {
//             mSortedPartSize++;
//             return mData.insert(sorted_part_end, TPointerType(new TOtherDataType(rData)));
//         }
// 
//         if (!EqualKeyTo(key)(*i))
//             if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(key))) == mData.end())
//             {
//                 mData.push_back(TPointerType(new TOtherDataType(rData)));
//                 //return iterator(--mData.end());
//                 return iterator(mData.end()-1);
//             }
//         **i = rData;
//         return i;
//     }

    iterator insert(iterator Position, const TPointerType pData)
    {
        ptr_iterator sorted_part_end;

        key_type key = KeyOf(*pData);

        if(mData.size() - mSortedPartSize >= mMaxBufferSize)
        {
            Sort();
            sorted_part_end = mData.end();
        }
        else
            sorted_part_end	= mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, key, CompareKey()));
        if (i == sorted_part_end)
        {
            mSortedPartSize++;
            return mData.insert(sorted_part_end, pData);
        }

        if (!EqualKeyTo(key)(*i))
            if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(key))) == mData.end())
            {
                mData.push_back(pData);
                return iterator(mData.end()-1);
// 	  return iterator(--mData.end());
            }

        *i = pData;
        return i;
    }

    template <class InputIterator>
    void insert(InputIterator First, InputIterator Last)
    {
        for(; First != Last; ++First)
            insert(begin(),*First);
    }


    iterator erase(iterator pos)
    {
		if (pos.base() == mData.end())
			return mData.end();
        iterator new_end = iterator( mData.erase( pos.base() ) );
        mSortedPartSize = mData.size();
        return new_end;
    }

    iterator erase( iterator first, iterator last )
    {
        iterator new_end = iterator( mData.erase( first.base(), last.base() ) );
        // TODO: Sorted part size must change
        mSortedPartSize = mData.size();
        return new_end;
    }

    iterator erase(const key_type& k)
    {
        return erase(find(k));
    }

    void clear()
    {
        mData.clear();
        mSortedPartSize = size_type();
        mMaxBufferSize = 1;
    }

    iterator find(const key_type& Key)
    {
        ptr_iterator sorted_part_end;

        if(mData.size() - mSortedPartSize >= mMaxBufferSize)
        {
            Sort();
            sorted_part_end = mData.end();
        }
        else
            sorted_part_end	= mData.begin() + mSortedPartSize;

        ptr_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end || (!EqualKeyTo(Key)(*i)))
            if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
                return mData.end();

        return i;
    }

    const_iterator find(const key_type& Key) const
    {
        ptr_const_iterator sorted_part_end(mData.begin() + mSortedPartSize);

        ptr_const_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
        if (i == sorted_part_end || (!EqualKeyTo(Key)(*i)))
            if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
                return mData.end();

        return const_iterator(i);
    }

    size_type count(const key_type& Key)
    {
        return find(Key) == mData.end() ? 0 : 1;
    }

    void reserve(int reservedsize)
    {
        mData.reserve(reservedsize);
    }
    int capacity()
    {
        return mData.capacity();
    }

    void Sort()
    {
        std::sort(mData.begin(), mData.end(), CompareKey());
        mSortedPartSize = mData.size();
    }

    void Unique()
    {
        typename TContainerType::iterator end_it = mData.end();
//#ifndef _OPENMP
        std::sort(mData.begin(), mData.end(), CompareKey());
//#else
//	if(mData.size() > 2000)
//	  omptl::sort(mData.begin(), mData.end(), CompareKey());
//	else
//	  std::sort(mData.begin(), mData.end(), CompareKey());
//#endif
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


    /** Set the maximum size of buffer used in the container.

    This container uses a buffer which keep data unsorted. After
    buffer size arrived to the MaxBufferSize it will sort all
    container and empties buffer.

    @param NewSize Is the new buffer maximum size. */
    void SetMaxBufferSize(size_type NewSize)
    {
        mMaxBufferSize = NewSize;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Pointer vector set (size = " << size() << ") : ";

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

    class CompareKey : public std::binary_function<TPointerType, TPointerType, typename TCompareType::result_type>
    {
    public:
        typename TCompareType::result_type operator()(key_type a, TPointerType b) const
        {
            return TCompareType()(a, TGetKeyType()(*b));
        }
        typename TCompareType::result_type operator()(TPointerType a, key_type b) const
        {
            return TCompareType()(TGetKeyType()(*a), b);
        }
        typename TCompareType::result_type operator()(TPointerType a, TPointerType b) const
        {
            return TCompareType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };

    class EqualKeyTo : public std::binary_function<TPointerType, TPointerType, typename TEqualType::result_type>
    {
        key_type mKey;
    public:
        EqualKeyTo() : mKey() {}
        EqualKeyTo(key_type Key) : mKey(Key) {}
        typename TEqualType::result_type operator()(TPointerType a) const
        {
            return TEqualType()(mKey, TGetKeyType()(*a));
        }
        typename TEqualType::result_type operator()(TPointerType a, TPointerType b) const
        {
            return TEqualType()(TGetKeyType()(*a), TGetKeyType()(*b));
        }
    };
//        static typename TCompareType::result_type CompareKey(TDataType const & a, TDataType const & b)
//        {return TCompareType()(KeyOf(a), KeyOf(b));}

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

#endif // KRATOS_POINTER_VECTOR_SET_H_INCLUDED  defined 
