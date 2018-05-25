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



#if !defined(KRATOS_POINTER_HASH_MAP_SET_H_INCLUDED )
#define  KRATOS_POINTER_HASH_MAP_SET_H_INCLUDED



// System includes


// External includes
#include <unordered_map>

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

/// PointerHashMapSet is a hash implemenetation of the PointerVectorSet.
/** This container is like a set but is built over a hash map in order 
	to allow the key to be a part of the value. It is important to mention
	that the value is not constant and if the key inside the value changed
	outside results in inconsistence condition.

    This Container does not free the memory by
    itself and relies on using of counted pointers or manual
    deleting.
 */
template<class TDataType,
         class THashType = std::hash<TDataType>,
         class TGetKeyType = SetIdentityFunction<TDataType>,
         class TPointerType = Kratos::shared_ptr<TDataType> >
class PointerHashMapSet
{


public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointerHashMapSet
    KRATOS_CLASS_POINTER_DEFINITION(PointerHashMapSet);

    /// Key type for searching in this container.
    typedef typename TGetKeyType::result_type key_type;

    /// data type stores in this container.
    typedef TDataType data_type;
    typedef TDataType value_type;
    typedef THashType hasher;
    typedef TPointerType pointer_type;
    typedef TDataType& reference;
    typedef const TDataType& const_reference;
	typedef std::unordered_map<key_type, TPointerType, hasher> ContainerType;

    typedef typename ContainerType::size_type size_type;
    typedef typename ContainerType::iterator ptr_iterator;
    typedef typename ContainerType::const_iterator ptr_const_iterator;
    typedef typename ContainerType::difference_type difference_type;

	///@}

private:
	///@name Nested clases
	///@{
	class iterator_adaptor : public std::iterator<std::forward_iterator_tag, data_type>
	{
		ptr_iterator map_iterator;
	public:
		iterator_adaptor(ptr_iterator it) :map_iterator(it) {}
		iterator_adaptor(const iterator_adaptor& it) : map_iterator(it.map_iterator) {}
		iterator_adaptor& operator++() { map_iterator++; return *this; }
		iterator_adaptor operator++(int) { iterator_adaptor tmp(*this); operator++(); return tmp; }
		bool operator==(const iterator_adaptor& rhs) const { return map_iterator == rhs.map_iterator; }
		bool operator!=(const iterator_adaptor& rhs) const { return map_iterator != rhs.map_iterator; }
		data_type& operator*() const { return *(map_iterator->second); }
		pointer_type operator->() const { return map_iterator->second; }
		ptr_iterator& base() { return map_iterator; }
		ptr_iterator const& base() const { return map_iterator; }
	};

	class const_iterator_adaptor : public std::iterator<std::forward_iterator_tag, data_type>
	{
		ptr_const_iterator map_iterator;
	public:
		const_iterator_adaptor(ptr_const_iterator it) :map_iterator(it) {}
		const_iterator_adaptor(const const_iterator_adaptor& it) : map_iterator(it.map_iterator) {}
		const_iterator_adaptor& operator++() { map_iterator++; return *this; }
		const_iterator_adaptor operator++(int) { const_iterator_adaptor tmp(*this); operator++(); return tmp; }
		bool operator==(const const_iterator_adaptor& rhs) const { return map_iterator == rhs.map_iterator; }
		bool operator!=(const const_iterator_adaptor& rhs) const { return map_iterator != rhs.map_iterator; }
		data_type const& operator*() const { return *(map_iterator->second); }
		pointer_type operator->() const { return map_iterator->second; }
		ptr_const_iterator& base() { return map_iterator; }
		ptr_const_iterator const& base() const { return map_iterator; }
	};

	///@}

public:
	///@name Type Definitions
	///@{


	typedef iterator_adaptor iterator;
	typedef const_iterator_adaptor const_iterator;

	///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointerHashMapSet() : mData() {}

    template <class TInputIteratorType>

    PointerHashMapSet(TInputIteratorType First, TInputIteratorType Last, size_type NewMaxBufferSize = 1)
    {
        for(; First != Last; ++First)
            insert(begin(), *First);
    }

    PointerHashMapSet(const PointerHashMapSet& rOther)
        :  mData(rOther.mData) {}

    PointerHashMapSet(const ContainerType& rContainer) :  mData(rContainer)
    {
    }

    /// Destructor.
    virtual ~PointerHashMapSet() {}


    ///@}
    ///@name Operators
    ///@{

    PointerHashMapSet& operator=(const PointerHashMapSet& rOther)
    {
        mData = rOther.mData;

		return *this;
    }

    TDataType& operator[](const key_type& Key)
    {
		return *(mData[Key].second);
    }

    pointer_type& operator()(const key_type& Key)
    {
		return mData[Key].second;
	}

    bool operator==( const PointerHashMapSet& r ) const // nothrow
    {
       return (mData == r.mData);
    }

	bool operator!=(const PointerHashMapSet& r) const // nothrow
	{
		return (mData != r.mData);
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

    reference        front()       /* nothrow */
    {
        //assert( !empty() );
        return *(mData.front().second);
    }
    const_reference  front() const /* nothrow */
    {
        //assert( !empty() );
        return *(mData.front().second);
    }
    reference        back()        /* nothrow */
    {
        //assert( !empty() );
        return *(mData.back().second);
    }
    const_reference  back() const  /* nothrow */
    {
        //assert( !empty() );
        return *(mData.back().second);
    }

    size_type size() const
    {
        return mData.size();
    }

    //size_type max_size() const
    //{
    //    return mData.max_size();
    //}

    //key_compare key_comp() const
    //{
    //    return TCompareType();
    //}

    void swap(PointerHashMapSet& rOther)
    {
        mData.swap(rOther.mData);
    }

    template<class TOtherDataType>
    iterator insert(const TOtherDataType& rData)
    {
		TDataType* p_new_data = new TDataType(rData);
		return mData.insert(ContainerType::value_type(TGetKeyOf(rData), p_new_data));
    }

    iterator insert(TPointerType pData)
    {
		std::string key=KeyOf(*pData);
		typename ContainerType::value_type item(key, pData);
		std::pair<typename ContainerType::iterator, bool> result = mData.insert(item);
	// TODO: I should enable this after adding the KRATOS_ERROR to define.h. Pooyan.
	//if(result.second != true)
	//	KRATOS_ERROR << "Error in adding the new item" << std::endl
		return result.first;
	}

    template <class InputIterator>
    void insert(InputIterator First, InputIterator Last)
    {
        for(; First != Last; ++First)
            insert(begin(),*First);
    }


	iterator erase(iterator pos)
	{
		return mData.erase(pos.base());
	}

	size_type erase(key_type const& Key)
	{
		return mData.erase(Key);
	}

	iterator erase( iterator first, iterator last )
    {
        return mData.erase( first.base(), last.base() );
    }

    void clear()
    {
        mData.clear();
    }

    iterator find(const key_type& Key)
    {
         return mData.find(Key);
    }

    const_iterator find(const key_type& Key) const
    {
		return mData.find(Key);
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

 
    ///@}
    ///@name Access
    ///@{

    /** Gives a reference to underly normal container. */
    ContainerType& GetContainer()
    {
        return mData;
    }

    /** Gives a constant reference to underly normal container. */
    const ContainerType& GetContainer() const
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
        buffer << "Pointer hash map set (size = " << size() << ") : ";

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

    ContainerType mData;

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
        size_type size = mData.size();

        rSerializer.save("size", size);

		for (ptr_const_iterator i = ptr_begin(); i != ptr_end(); i++)
            rSerializer.save("E", i->second);
    }

    virtual void load(Serializer& rSerializer)
    {
        size_type size;

        rSerializer.load("size", size);

		for (size_type i = 0; i < size; i++)
		{
 		        pointer_type p = nullptr;// new TDataType;
			rSerializer.load("E", p);
			insert(p);
		}
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

}; // Class PointerHashMapSet

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
//template<class TDataType,
//         class TGetKeyType,
//         class TCompareType,
//         class TEqualType,
//         class TPointerType,
//         class ContainerType>
//inline std::istream& operator >> (std::istream& rIStream,
//                                  PointerHashMapSet<TDataType, TGetKeyType, TCompareType, TEqualType, TPointerType, ContainerType>& rThis);

/// output stream function
template<class TDataType,
         class TGetKeyType,
         class TCompareType,
         class TEqualType,
         class TPointerType,
         class ContainerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PointerHashMapSet<TDataType, TGetKeyType, TCompareType, TPointerType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_POINTER_HASH_MAP_SET_H_INCLUDED  defined 
