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
#include <iterator>
#include <vector>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"


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
  
  /// Buffer holds the binary value of the data.
  /** Detail class definition.
  */
  class Buffer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Buffer
      KRATOS_CLASS_POINTER_DEFINITION(Buffer);
      

      typedef std::size_t IndexType;

      typedef std::size_t SizeType;

      typedef double BlockType;
  
      /// Type of the container used for storing values 
      typedef BlockType* ContainerType;
      
      typedef BlockType* iterator;
      typedef BlockType const* const_iterator;
  
          
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Buffer(SizeType NewSize = 1) :  mpData(0) , mpBegin(mpData), mpEnd(mpData), mSize(0)
      {
	resize(BlockCompatibleSize(NewSize));

	for(SizeType i = 0 ; i < mSize ; i++)
	  mpData[i] = BlockType();
      }
    
      Buffer(const Buffer& rOther) :  mpData(new BlockType[rOther.mSize]), mpBegin(0), mpEnd(0), mSize(rOther.mSize) 
      {
	      // Setting the current position with relative source container offset
	      mpBegin = mpData + (rOther.mpBegin - rOther.mpData);
	      mpEnd = mpData + (rOther.mpEnd - rOther.mpData);
	      std::copy(rOther.mpBegin, rOther.mpEnd, mpBegin); 
      }
    
      Buffer(const ContainerType& rContainer, SizeType NewSize) : mpData(new BlockType[BlockCompatibleSize(NewSize)]) , mpBegin(mpData), mpEnd(mpData+NewSize), mSize(BlockCompatibleSize(NewSize))
      {
 	   std::copy(rContainer, rContainer + mSize, mpData); 
      }

      /// Destructor.
      virtual ~Buffer()
      {
	free(mpData);
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
    
      Buffer& operator=(const Buffer& rOther)
      {
	// here I'm just copying the active part of the buffer not the whole!
	const SizeType other_size=rOther.mpEnd - rOther.mpBegin;
	if(mSize < other_size)
	{
	  mSize=other_size;
	  mpData = (BlockType*)realloc(mpData, mSize * sizeof(BlockType));
	}
	std::copy(rOther.mpBegin, rOther.mpEnd, mpData); 
	mpBegin = mpData;
	mpEnd = mpData+mSize;

	return *this;
      }
    
    
      ///@}
      ///@name Operations
      ///@{
      
      SizeType size() const {return mSize * sizeof(BlockType);}
    
      void swap(Buffer& rOther) 
      {
	 std::swap(mSize, rOther.mSize);
	 std::swap(mpData, rOther.mpData);
	 std::swap(mpBegin, rOther.mpBegin);
	 std::swap(mpEnd, rOther.mpEnd);
      }
    
      void push_back(std::string const& rValue)                
      {
	std::size_t string_size = rValue.size() + 1; // the one is for end string null character
	SizeType data_size = BlockCompatibleSize(string_size*sizeof(char)); 
	iterator new_end = mpEnd + data_size;
	iterator max_end = mpData + mSize;
	if(new_end > max_end)
	{
	  resize(mSize+data_size);
	}
	std::copy(rValue.c_str(), rValue.c_str() + string_size, (char*)mpEnd);
	mpEnd+=data_size;
      } 
    
      void push_back(bool rValue)                
      {
	int temp(rValue);
	push_back(temp);
      } 
    
      template<class TDataType>
      void push_back(std::vector<TDataType> const& rValue)                
      {
	const SizeType size=rValue.size();
	push_back(size);
	push_back(rValue.begin(), rValue.end());
      } 
    
      template<class TDataType>
      void push_back(boost::numeric::ublas::vector<TDataType> const& rValue)                
      {
	push_back(rValue.size());
	push_back(rValue.begin(), rValue.end());
      } 
    
      template<class TDataType, std::size_t TDimenasion>
      void push_back(array_1d<TDataType,TDimenasion> const& rValue)                
      {
	push_back(rValue.size());
	push_back(rValue.begin(), rValue.end());
      } 
    
      void push_back(Matrix const& rValue)                
      {
	push_back(rValue.size1());
	push_back(rValue.size2());
	push_back(rValue.data().begin(), rValue.data().end());
      } 
    
      template<class TDataType>
      void push_back(boost::shared_ptr<TDataType> const& rValue)                
      {
	KRATOS_ERROR(std::logic_error, "You cannot store a pointer in the buffer try the Serializer instead", "" );
      } 
    
      template<class TDataType>
      void push_back(TDataType const& rValue)                
      {
	SizeType data_size = BlockCompatibleSize(sizeof(TDataType));
	iterator new_end = mpEnd + data_size;
	iterator max_end = mpData + mSize;
	if(new_end > max_end)
	{
	  resize(mSize+data_size);
	}
	
	*((TDataType *)mpEnd) = rValue;
	mpEnd+=data_size;
      } 
    
      template<class TIteratorType>
      void push_back(TIteratorType First, TIteratorType Last)                
      {
	for(;First != Last ; First++)
	  push_back(*First);
      } 
    
      std::string& pop_front(std::string& rValue)                
      {
	rValue= static_cast<char*>((void*)mpBegin);
	const std::size_t string_size=rValue.size() + 1;;
	SizeType data_size = BlockCompatibleSize(string_size);
	
	mpBegin += data_size;
	return rValue;
      } 
        
      void pop_front(bool& rValue)                
      {
	int temp;
	
	pop_front(temp);
	
 	rValue = temp << 1;
      } 
        
      template<class TDataType>
      void pop_front(std::vector<TDataType>& rValue)                
      {
	SizeType size;
	
	pop_front(size);
	
 	rValue.resize(size);
	
	pop_front(rValue.begin(), rValue.end());
      } 
        
      template<class TDataType>
      void pop_front(boost::numeric::ublas::vector<TDataType>& rValue)                
      {
	SizeType size;
	
	pop_front(size);
	
 	rValue.resize(size);
	
	pop_front(rValue.begin(), rValue.end());
      } 
        
      template<class TDataType, std::size_t TDimenasion>
      void pop_front(array_1d<TDataType, TDimenasion>& rValue)                
      {
	SizeType size;
	
	pop_front(size);
	
 	rValue.resize(size);
	
	pop_front(rValue.begin(), rValue.end());
      } 
       
      void pop_front(Matrix& rValue)                
      {
	SizeType size1;
	SizeType size2;
	
	pop_front(size1);
	pop_front(size2);
	
 	rValue.resize(size1,size2);
	
	pop_front(rValue.data().begin(), rValue.data().end());
      } 
    
      template<class TDataType>
      TDataType& pop_front(TDataType& rValue)                
      {
	SizeType data_size = BlockCompatibleSize(sizeof(TDataType));
	
	rValue = *static_cast<TDataType*>((void*)mpBegin);
	mpBegin += data_size;
	return rValue;
      } 
        
      template<class TIteratorType>
      void pop_front(TIteratorType First, TIteratorType Last)                
      {
	for(;First != Last ; First++)
	  pop_front(*First);
      } 
    
       

      void clear() 
      {
	mpBegin = mpData;
	mpEnd = mpData;
      }
    
      void resize(SizeType NewSize)
      {
	if(mSize < NewSize)
	{
	  mSize=NewSize;
	  
	  SizeType begin_offset = mpBegin - mpData;
	  SizeType end_offset = mpEnd - mpData;
	  if(mpData)
	    mpData = (BlockType*)realloc(mpData, mSize * sizeof(BlockType));
	  else
	    mpData = (BlockType*)malloc(mSize * sizeof(BlockType));
	    
	  mpBegin = mpData + begin_offset;
	  mpEnd = mpData + end_offset;
	}
      }

    
      ///@}
      ///@name Access
      ///@{ 
      
       iterator                   begin()            { return mpBegin; }
       const_iterator             begin() const      { return mpBegin; }
	
       iterator                   end()              { return mpEnd; }
       const_iterator             end()   const      { return mpEnd; }
	
     
      
      ///@}
      ///@name Inquiry
      ///@{
      
      bool empty() const {return mpBegin==mpEnd;}
    
      
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
	std::copy(begin(), end(), std::ostream_iterator<BlockType>(rOStream));
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
        
      ContainerType mpData;
       
      iterator mpBegin;
      
      iterator mpEnd;
      
      SizeType mSize;

      ///@} 
      ///@name Private Operators
      ///@{
      
	inline SizeType BlockCompatibleSize(SizeType rSize)
	{
	  const SizeType block_size = sizeof(BlockType);
	  return static_cast<SizeType>(((block_size - 1) + rSize) / block_size);	  
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
        
    }; // Class Buffer 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Buffer& rThis)
  {
      return rIStream;
  }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Buffer& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
   
  ///@} 
  
  template<class TDataType>
  inline Buffer& operator >> (Buffer& rThis, TDataType& rValue)
  {
    rThis.pop_front(rValue);
    
    return rThis;
  }

  
  template<class TDataType>
  inline Buffer& operator << (Buffer& rThis, const TDataType& rValue)
  {
    rThis.push_back(rValue);
    
    return rThis;
  }
  
  inline Buffer& operator >> (Buffer& rThis, void*& rValue)
  {
    rThis.pop_front<void*>(rValue);
    
    return rThis;
  }

  
  inline Buffer& operator << (Buffer& rThis, void* rValue)
  {
    rThis.push_back<void*>(rValue);
    
    return rThis;
  }
  
  inline Buffer& operator >> (Buffer& rThis, bool& rValue)
  {
    rThis.pop_front(rValue);
    
    return rThis;
  }

  
  inline Buffer& operator << (Buffer& rThis, bool rValue)
  {
    rThis.push_back(rValue);
    
    return rThis;
  }
  
  
/*  inline Buffer& operator << (Buffer& rThis, int const& rValue)
  {
    rThis.push_back(static_cast<Buffer::BlockType>(rValue));
    
    return rThis;
  }
  
  inline Buffer& operator << (Buffer& rThis, double const& rValue)
  {
    rThis.push_back(rValue);
    
    return rThis;
  }*/
  

}  // namespace Kratos.

#endif // KRATOS_BUFFER_H_INCLUDED  defined 


