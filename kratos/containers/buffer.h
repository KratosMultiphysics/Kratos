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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_BUFFER_H_INCLUDED )
#define  KRATOS_BUFFER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


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
  
  /// Short class definition.
  /** Detail class definition.
  */
  template<class TContainerType>
  class Buffer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Buffer
      KRATOS_CLASS_POINTER_DEFINITION(Buffer);
  
      typedef typename TContainerType::value_type value_type;
      typedef typename TContainerType::size_type size_type;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef TContainerType ContainerType;

      typedef typename TContainerType::iterator                iterator;
      typedef typename TContainerType::const_iterator          const_iterator;                
      typedef typename TContainerType::reverse_iterator        reverse_iterator;
      typedef typename TContainerType::const_reverse_iterator  const_reverse_iterator;    
          
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Buffer() : mData(), mCurrentIndex(0) {}
    
      template <class TInputIteratorType>
      Buffer(TInputIteratorType First, TInputIteratorType Last) 
	: mData(First, Last), mCurrentIndex(0)
      {
      }
      
      Buffer(const Buffer& rOther) :  mData(rOther.mData), mCurrentIndex(rOther.mCurrentIndex) {}
    
      Buffer(const TContainerType& rContainer) :  mData(rContainer) , mCurrentIndex(0)
      {
      }

      Buffer(std::size_t NewSize) :  mData(NewSize) , mCurrentIndex(0)
      {
      }

      template<class TOtherDataType>
      Buffer(std::size_t NewSize, TOtherDataType const& Value) :  mData(NewSize), mCurrentIndex(0)
      {
	for(size_type i = 0 ; i < NewSize ; i++)
	  mData[i] = pointer(new TOtherDataType(Value));
      }
    

      /// Destructor.
      virtual ~Buffer(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
    
      Buffer& operator=(const Buffer& rOther)
      {
	mData = rOther.mData;
	mCurrentIndex = rOther.mCurrentIndex;
	return *this;
      }
    
      value_type& operator[](const size_type& i)
      {
	return mData[TransformIndex(i)];
      }
    
      value_type const& operator[](const size_type& i) const
      {
	return mData[TransformIndex(i)];
      }
    
      ///@}
      ///@name Operations
      ///@{
      
      iterator                   begin()            { return mData.begin(); }
      const_iterator             begin() const      { return mData.begin(); }
      iterator                   end()              { return mData.end(); }
      const_iterator             end() const        { return mData.end(); }
      
      reference        front()       /* nothrow */ { assert( !empty() ); return mData[TransformIndex(0)]; }
      const_reference  front() const /* nothrow */ { assert( !empty() ); return mData[TransformIndex(0)]; }
      reference        back()        /* nothrow */ { assert( !empty() ); return mData[TransformIndex(size() - 1)]; }
      const_reference  back() const  /* nothrow */ { assert( !empty() ); return mData[TransformIndex(size() - 1)]; }

      size_type size() const {return mData.size();}
    
      size_type max_size() const {return mData.max_size();}
    
      void swap(Buffer& rOther) {mData.swap(rOther.mData);}
    
      void push_back(value_type& x)                
      {
	std::size_t index = TransformIndex(size());
	mData[index] = x;
	mCurrentIndex = index;
      } 
        
      void push_front(value_type const& x)                
      {
	std::size_t index = (mCurrentIndex == 0) ? size() - 1 :  mCurrentIndex - 1;
// 	std::cout << "push_front..." << std::endl;
// 	KRATOS_WATCH(mCurrentIndex)
// 	KRATOS_WATCH(size())
// 	KRATOS_WATCH(index)
// /* 	KRATOS_WATCH(mData) */
// 	KRATOS_WATCH(x)
// 	KRATOS_WATCH(mData[index])

	mData[index] = x;
	mCurrentIndex = index;
      } 
        
      iterator insert(iterator Position, const value_type& Data)
      {
	return mData.insert(Position, Data);
      }

      template <class InputIterator>
      void insert(InputIterator First, InputIterator Last)
      {
	for(;First != Last; ++First)
	  insert(*First);
      }
    
                       
      iterator erase(iterator pos) {return mData.erase(pos.base());}
    
      iterator erase( iterator first, iterator last )
      {
	return mData.erase( first.base(), last.base());
      }

      void clear() {mData.clear();}
    
      void resize(std::size_t NewSize)
      {
	mData.resize(NewSize);
      }

      void resize(std::size_t NewSize, value_type& Data)
      {
	mData.resize(NewSize, Data);
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
      
      bool empty() const {return mData.empty();}
    
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	std::stringstream buffer;
	buffer << "buffer (size = " << size() << ") : ";
	
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
	std::copy(begin(), end(), std::ostream_iterator<value_type>(rOStream, "\n "));
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
       
      std::size_t mCurrentIndex;

      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      inline std::size_t TransformIndex(std::size_t ThisIndex)
      {
	std::size_t index = ThisIndex + mCurrentIndex;
	return (index < mData.size()) ? index : index - mData.size();
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
        
    }; // Class Buffer 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TContainerType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Buffer<TContainerType>& rThis);

  /// output stream function
  template<class TContainerType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Buffer<TContainerType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_BUFFER_H_INCLUDED  defined 


