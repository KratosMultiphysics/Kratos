/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_VECTOR_MAP_H_INCLUDED )
#define  KRATOS_VECTOR_MAP_H_INCLUDED



// System includes
#include <vector>
#include <functional> 
#include <string>
#include <iostream>
#include <sstream> 
#include <cstddef>

// External includes 


// Project includes
#include "includes/define.h"
#include "containers/vector_map_iterator.h"


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
  
  /// VectorMap is a sorted associative container like stl map but using a vector to store its data.
  /** VectorMap is a sorted associative container like stl map
      but using a vector to store data. 
   */
  template<class TKeyType,class TDataType,
	    
	   class TCompareType = std::less<TKeyType>, 
	   class TContainerType = std::vector<std::pair<TKeyType, TDataType> > >
  class VectorMap
  {
  public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of VectorMap
      KRATOS_CLASS_POINTER_DEFINITION(VectorMap);
    
    /// Key type for searching in this container. 
    typedef TKeyType key_type;

    /// data type stores in this container. 
    typedef TDataType data_type;
    typedef std::pair<TKeyType, TDataType> value_type;
    typedef TCompareType key_compare;
    typedef TDataType& reference;
    typedef const TDataType& const_reference;
    typedef TContainerType ContainerType;

    typedef VectorMapIterator<typename TContainerType::iterator, TDataType>                iterator;
    typedef VectorMapIterator<typename TContainerType::const_iterator, TDataType>          const_iterator;                
    typedef VectorMapIterator<typename TContainerType::reverse_iterator, TDataType>        reverse_iterator;
    typedef VectorMapIterator<typename TContainerType::const_reverse_iterator, TDataType>  const_reverse_iterator;              

    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::iterator pair_iterator;
    typedef typename TContainerType::const_iterator pair_const_iterator;
    typedef typename TContainerType::reverse_iterator pair_reverse_iterator;
    typedef typename TContainerType::const_reverse_iterator pair_const_reverse_iterator;
    
    ///@}
    ///@name Life Cycle 
    ///@{ 
    
    /// Default constructor.
    VectorMap() : mData(), mSortedPartSize(size_type()), mMaxBufferSize(100) {}
    
    VectorMap(const VectorMap& rOther) :  mData(rOther.mData), mSortedPartSize(rOther.mSortedPartSize), mMaxBufferSize(rOther.mMaxBufferSize) {}
    
    VectorMap(const TContainerType& rContainer) :  mData(rContainer), mSortedPartSize(size_type()), mMaxBufferSize(100) 
    {
      Sort();
      std::unique(mData.begin(), mData.end(), EqualKeyTo());
    }
    
    /// Destructor.
    virtual ~VectorMap(){}
      

    ///@}
    ///@name Operators 
    ///@{
    
    VectorMap& operator=(const VectorMap& rOther)
    {
      mData = rOther.mData;
      mSortedPartSize = rOther.mSortedPartSize;

      return *this;
    }
    
	TDataType& operator[](const key_type& Key)
	{
		pair_iterator sorted_part_end;

		if(mData.size() - mSortedPartSize >= mMaxBufferSize)
		{
			Sort();
			sorted_part_end = mData.end();
		}
		else
			sorted_part_end	= mData.begin() + mSortedPartSize;

		pair_iterator i(std::lower_bound(mData.begin(), sorted_part_end, value_type(Key,data_type()), CompareKey()));
		if (i == sorted_part_end)
		{
			mSortedPartSize++;
			return mData.insert(sorted_part_end, value_type(Key, TDataType()))->second;
		}

		if (Key != i->first)
			if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
			{
				mData.push_back(value_type(Key, TDataType()));
				return ((--mData.end())->second);
			}

			return (i->second);
	}
    
    TDataType& operator()(const key_type& Key)
    {
		pair_iterator sorted_part_end;

		if(mData.size() - mSortedPartSize >= mMaxBufferSize)
		{
			Sort();
			sorted_part_end = mData.end();
		}
		else
			sorted_part_end	= mData.begin() + mSortedPartSize;

		pair_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
		if (i == sorted_part_end)
		{
			mSortedPartSize++;
			return (mData.insert(sorted_part_end, value_type(Key, TDataType()))->second);
		}

		if (Key != i->first)
			if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
			{
				mData.push_back(value_type(Key, TDataType()));
				return ((--mData.end())->second);
			}

			return (i->second);
    }
    
        
    ///@}
    ///@name Operations
    ///@{
    
    iterator                   begin()            { return iterator( mData.begin() ); }
    const_iterator             begin() const      { return const_iterator( mData.begin() ); }
    iterator                   end()              { return iterator( mData.end() ); }
    const_iterator             end() const        { return const_iterator( mData.end() ); }
    reverse_iterator           rbegin()           { return reverse_iterator( mData.rbegin() ); } 
    const_reverse_iterator     rbegin() const     { return const_reverse_iterator( mData.rbegin() ); } 
    reverse_iterator           rend()             { return reverse_iterator( mData.rend() ); } 
    const_reverse_iterator     rend() const       { return const_reverse_iterator( mData.rend() ); } 
    pair_iterator               pair_begin()        { return mData.begin(); }
    pair_const_iterator         pair_begin() const  { return mData.begin(); }
    pair_iterator               pair_end()          { return mData.end(); }
    pair_const_iterator         pair_end() const    { return mData.end(); }
    pair_reverse_iterator       pair_rbegin()       { return mData.rbegin(); }
    pair_const_reverse_iterator pair_rbegin() const { return mData.rbegin(); }
    pair_reverse_iterator       pair_rend()         { return mData.rend(); }
    pair_const_reverse_iterator pair_rend() const   { return mData.rend(); }

    reference        front()       /* nothrow */ { assert( !empty() ); return (mData.front().second); }
    const_reference  front() const /* nothrow */ { assert( !empty() ); return (mData.front().second); }
    reference        back()        /* nothrow */ { assert( !empty() ); return (mData.back().second); }
    const_reference  back() const  /* nothrow */ { assert( !empty() ); return (mData.back().second); }
    
    size_type size() const {return mData.size();}
    
    size_type max_size() const {return mData.max_size();}
    
    key_compare key_comp() const {return TCompareType();}
    
    void swap(VectorMap& rOther) {mData.swap(rOther.mData);}
    
    void push_back(value_type x)                
    {
      mData.push_back(x);
    } 
        
    iterator insert(key_type const& Key, const TDataType& rData)
    {
		pair_iterator sorted_part_end;

		if(mData.size() - mSortedPartSize >= mMaxBufferSize)
		{
			Sort();
			sorted_part_end = mData.end();
		}
		else
			sorted_part_end	= mData.begin() + mSortedPartSize;

		pair_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
		if (i == sorted_part_end)
		{
			mSortedPartSize++;
			return (mData.insert(sorted_part_end, value_type(Key, TDataType(rData))));
		}

		if (Key != i->first)
			if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
			{
				mData.push_back(value_type(Key, TDataType()));
				return iterator(--mData.end());
			}
			(i->second) = rData;
			return i;
    }

    iterator erase(iterator pos) {return iterator(mData.erase(pos.base()));}
    
    iterator erase( iterator first, iterator last )
    {
      return iterator( mData.erase( first.base(), last.base() ) );
    }

    iterator erase(const key_type& k) {return erase(find(k));}
    
    void clear() {mData.clear();}
    
    iterator find(const key_type& Key)
    {
      pair_iterator sorted_part_end;

      if(mData.size() - mSortedPartSize >= mMaxBufferSize)
      {
	Sort();
	sorted_part_end = mData.end();
      }
      else
        sorted_part_end	= mData.begin() + mSortedPartSize;

      pair_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
       if (i == sorted_part_end || (Key != i->first))
			if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
	  return mData.end();
      
      return i;
    }

    const_iterator find(const key_type& Key) const
    {
      pair_const_iterator sorted_part_end(mData.begin() + mSortedPartSize);

       pair_iterator i(std::lower_bound(mData.begin(), sorted_part_end, Key, CompareKey()));
       if (i == sorted_part_end || (Key != i->first))
			if((i = std::find_if(sorted_part_end, mData.end(), EqualKeyTo(Key))) == mData.end())
	  return mData.end();
      
      return const_iterator(i);
    }

    size_type count(const key_type& Key)
    {
      return find(Key) == mData.end() ? 0 : 1;
    }
    
    
    void Sort()
    {
      std::sort(mData.begin(), mData.end(), CompareValue());
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
    
    bool empty() const {return mData.empty();}
    
    bool IsSorted() const {return (mSortedPartSize == mData.size());}
    
    ///@}      
    ///@name Input and output
    ///@{
    
      /// Turn back information as a string.
    virtual std::string Info() const
    {
      std::stringstream buffer;
      buffer << "vector map (size = " << size() << ") : ";
      
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
		for(typename TContainerType::const_iterator i = mData.begin() ; i != mData.end() ; i++)
			rOStream << "(" << i->first << " , " << (i->second) << ")" << std::endl;
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
     class CompareKey : public std::binary_function<value_type, value_type, bool>{
     public:
       bool operator()(value_type a, value_type b) const
       {return TCompareType()(a.first, b.first);}
       bool operator()(value_type a, key_type b) const
       {return TCompareType()(a.first, b);}
       bool operator()(key_type a, value_type b) const
       {return TCompareType()(a, b.first);}
     };
     class CompareValue : public std::binary_function<value_type, value_type, bool>{
     public:
       bool operator()(value_type a, value_type b) const
       {return TCompareType()(a.first, b.first);}
     };
     class EqualKeyTo : public std::binary_function<value_type, value_type, bool>{
		 key_type mKey;
     public:
		 EqualKeyTo(key_type k) : mKey(k){}
       bool operator()(value_type a, value_type b) const
       {return a.first == b.first;}
       bool operator()(value_type a) const
       {return a.first == mKey;}
     };

//        static typename TCompareType::result_type TCompareType(TDataType const & a, TDataType const & b)
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
        
    }; // Class VectorMap 

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
	   class TContainerType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    VectorMap<TDataType, TGetKeyType, TCompareType, TContainerType>& rThis);

  /// output stream function
  template<class TDataType,
	   class TGetKeyType, 
	   class TCompareType, 
	   class TContainerType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const VectorMap<TDataType, TGetKeyType, TCompareType, TContainerType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VECTOR_MAP_H_INCLUDED  defined 
