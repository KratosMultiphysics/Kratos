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

  ///@name Kratos Classes
  ///@{
  
  /// Buffer holds the binary value of the data.
  /** This buffer is desinged to be used via push_back and pop_front methods.
      By using push_back the binary representation of the class will be added
	  to the buffer. If the buffer size exceed an automatic realocation will
	  be performed. The stored values can be extracted using pop_front method.
	  The part of the buffer which is exctracted cannot be used again until the
	  clear method is called.
  */
  class Buffer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Buffer
      KRATOS_CLASS_POINTER_DEFINITION(Buffer);
      

	  /// Type used for indexing in the buffer
      typedef std::size_t IndexType;

	  /// Type used for returning the size of the buffer.
      typedef std::size_t SizeType;

	  /** The building block of the buffer which is a double.
	      This means that even a char will be stored in place of a double.
		  Note that the extracted values are in their original type independent
		  of this type. The reason of using double is compatibility with memory arrangment.
	  */
      typedef double BlockType;
  
      /// Type of the container used for storing values which is a raw pointer to BlockType
      typedef BlockType* ContainerType;
      
	  /// The iterator is defined as a raw pointer to BlockType
      typedef BlockType* iterator;

	  /// The constant iterator is defined as a raw pointer to constant BlockType
      typedef BlockType const* const_iterator;
  
          
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /** Constructor with buffer size. 
	  
		  @param NewSize This size is for reserving the memory and by default
	      is equal to 1.
	  */
	  Buffer(SizeType NewSize = 1) :  mpData(0) , mpBegin(mpData), mpEnd(mpData), mCapacity(0)
	  {
		  resize(BlockCompatibleSize(NewSize));

		  for(SizeType i = 0 ; i < mCapacity ; i++)
			  mpData[i] = BlockType();
	  }

	  /** Copy constructor which does a deep copy of other buffer.
	      
		  @param rOther The other buffer to be copied
	  */
      Buffer(const Buffer& rOther) :  mpData(new BlockType[rOther.mCapacity]), mpBegin(0), mpEnd(0), mCapacity(rOther.mCapacity)
      {
	      // Setting the current position with relative source container offset
	      mpBegin = mpData + (rOther.mpBegin - rOther.mpData);
	      mpEnd = mpData + (rOther.mpEnd - rOther.mpData);
	      std::copy(rOther.mpBegin, rOther.mpEnd, mpBegin); 
      }
    
	  /** Constructor which does a deep copy of data given as a memory block.
	      
		  @param rContainer A pointer to the constant memory block containing the data.
		  @param NewSize The size of the memory block
	  */
      Buffer(const ContainerType& rContainer, SizeType NewSize) : mpData(new BlockType[BlockCompatibleSize(NewSize)]) , mpBegin(mpData), mpEnd(mpData+NewSize), mCapacity(BlockCompatibleSize(NewSize))
      {
 	   std::copy(rContainer, rContainer + mCapacity, mpData);
      }

	  /// Destructor releases the memory allocated for buffer.
	  virtual ~Buffer()
	  {
		  free(mpData);
	  }
      

      ///@}
      ///@name Operators 
      ///@{
      
      /** The assignment operator makes a deep copy of other Buffer.
	      @param rOther The other buffer to be copied.
	  */
	  Buffer& operator=(const Buffer& rOther)
	  {
		  // here I'm just copying the active part of the buffer not the whole!
		  const SizeType other_size=rOther.mpEnd - rOther.mpBegin;
		  if(mCapacity < other_size)
		  {
			  mCapacity=other_size;
			  mpData = (BlockType*)realloc(mpData, mCapacity * sizeof(BlockType));
		  }
		  std::copy(rOther.mpBegin, rOther.mpEnd, mpData); 
		  mpBegin = mpData;
		  mpEnd = mpData+mCapacity;

		  return *this;
	  }

    
      ///@}
      ///@name Operations
      ///@{
      
	  /** Returns the size of the buffer.
	  */
          SizeType size() const {return ((mpEnd-mpBegin) * sizeof(BlockType));}

          SizeType capacity() const {return mCapacity * sizeof(BlockType);}

 
	  /** Swaps this buffer by another one by swapping the pointer 
	      to data which make it efficient.
	  */
	  void swap(Buffer& rOther) 
	  {
		  std::swap(mCapacity, rOther.mCapacity);
		  std::swap(mpData, rOther.mpData);
		  std::swap(mpBegin, rOther.mpBegin);
		  std::swap(mpEnd, rOther.mpEnd);
	  }

	  /// Adds a string to the end of the buffer.
	  void push_back(std::string const& rValue)                
	  {
		  std::size_t string_size = rValue.size() + 1; // the one is for end string null character
		  SizeType data_size = BlockCompatibleSize(string_size*sizeof(char)); 
		  iterator new_end = mpEnd + data_size;
		  iterator max_end = mpData + mCapacity;
		  if(new_end > max_end)
		  {
			  resize(mCapacity+data_size);
		  }
		  std::copy(rValue.c_str(), rValue.c_str() + string_size, (char*)mpEnd);
		  mpEnd+=data_size;
	  } 

	  /// Adds a bool to the end of the buffer.
	  void push_back(bool rValue)                
	  {
		  int temp(rValue);
		  push_back(temp);
	  } 
    
	  /** Adds a std::vector to the end of the buffer. This method stores the size and all elements of the vector.
	  */
	  template<class TDataType>
	  void push_back(std::vector<TDataType> const& rValue)                
	  {
		  const SizeType size=rValue.size();
		  push_back(size);
		  push_back(rValue.begin(), rValue.end());
	  } 

	  /** Adds a ublas::vector to the end of the buffer. This method stores the size and all elements of the vector.
	  */
	  template<class TDataType>
	  void push_back(boost::numeric::ublas::vector<TDataType> const& rValue)                
	  {
		  push_back(rValue.size());
		  push_back(rValue.begin(), rValue.end());
	  } 

	  /** Adds an array_1d to the end of the buffer. This method stores the size and all elements of the array.
	      The size is stored for compatibility with other vectors but can be omitted.
	  */
	  template<class TDataType, std::size_t TDimenasion>
	  void push_back(array_1d<TDataType,TDimenasion> const& rValue)                
	  {
		  push_back(rValue.size());
		  push_back(rValue.begin(), rValue.end());
	  } 

	  /** Adds a Matrix to the end of the buffer. This method stores the number of 
	      rows and columns and then all elements of the Matrix.
	  */
	  void push_back(Matrix const& rValue)                
	  {
		  push_back(rValue.size1());
		  push_back(rValue.size2());
		  push_back(rValue.data().begin(), rValue.data().end());
	  } 

	  /** This buffer is not prepared for storing pointers. So a logical error 
	      will be sent to notify the developers to use the serialization which
		  is in charge of storing pointers correctly.
	  */
	  template<class TDataType>
	  void push_back(boost::shared_ptr<TDataType> const& rValue)                
	  {
		  KRATOS_ERROR(std::logic_error, "You cannot store a pointer in the buffer try the Serializer instead", "" );
	  } 

	  /** A generic push back to cover all other type of data wich are not specified before.
	  */
	  template<class TDataType>
	  void push_back(TDataType const& rValue)                
	  {
		  SizeType data_size = BlockCompatibleSize(sizeof(TDataType));
		  iterator new_end = mpEnd + data_size;
		  iterator max_end = mpData + mCapacity;
		  if(new_end > max_end)
		  {
			  resize(mCapacity+data_size);
		  }

		  *((TDataType *)mpEnd) = rValue;
		  mpEnd+=data_size;
	  } 

	  /** Adds elements of a container specified by its iterators.
	  */
	  template<class TIteratorType>
	  void push_back(TIteratorType First, TIteratorType Last)                
	  {
		  for(;First != Last ; First++)
			  push_back(*First);
	  } 

	  /// Extracts a string from the beginning of the buffer.
	  std::string& pop_front(std::string& rValue)                
	  {
		  rValue= static_cast<char*>((void*)mpBegin);
		  const std::size_t string_size=rValue.size() + 1;;
		  SizeType data_size = BlockCompatibleSize(string_size);

		  mpBegin += data_size;
		  return rValue;
	  } 

	  /// Extracts a boolian from the beginning of the buffer.
	  void pop_front(bool& rValue)                
	  {
		  int temp;

		  pop_front(temp);

		  rValue = temp << 1;
	  } 


	  /** Extracts a std::vector from the beginning of the buffer.
	      It reads the size and resizes the vector. Then it fills the vector element by element.
	  */
	  template<class TDataType>
	  void pop_front(std::vector<TDataType>& rValue)                
	  {
		  SizeType size;

		  pop_front(size);

		  rValue.resize(size);

		  pop_front(rValue.begin(), rValue.end());
	  } 

	  /** Extracts a ublas::vector from the beginning of the buffer.
	      It reads the size and resizes the vector. Then it fills the vector element by element.
	  */
	  template<class TDataType>
	  void pop_front(boost::numeric::ublas::vector<TDataType>& rValue)                
	  {
		  SizeType size;

		  pop_front(size);

		  rValue.resize(size,false);

		  pop_front(rValue.begin(), rValue.end());
	  } 

	  /** Extracts a array_1d from the beginning of the buffer.
	  */
	  template<class TDataType, std::size_t TDimenasion>
	  void pop_front(array_1d<TDataType, TDimenasion>& rValue)                
	  {
		  /// TODO: I have to take out the resize. Pooyan.
		  SizeType size;

		  pop_front(size);

		  rValue.resize(size);

		  pop_front(rValue.begin(), rValue.end());
	  } 

	  /** Extracts a Matrix from the beginning of the buffer.
	      It reads the sizes and resizes the Matrix. Then it fills the Matrix element by element.
	  */
	  void pop_front(Matrix& rValue)                
	  {
		  SizeType size1;
		  SizeType size2;

		  pop_front(size1);
		  pop_front(size2);

		  rValue.resize(size1,size2);

		  pop_front(rValue.data().begin(), rValue.data().end());
	  } 

	  /** A generic pop front to cover all other type of data wich are not specified before.
	  */
	  template<class TDataType>
	  TDataType& pop_front(TDataType& rValue)                
	  {
		  SizeType data_size = BlockCompatibleSize(sizeof(TDataType));

		  rValue = *static_cast<TDataType*>((void*)mpBegin);
		  mpBegin += data_size;
		  return rValue;
	  } 

	  /** Extracts elements of a container specified by its iterators.
	  */
	  template<class TIteratorType>
	  void pop_front(TIteratorType First, TIteratorType Last)                
	  {
		  for(;First != Last ; First++)
			  pop_front(*First);
	  } 



	  /** Clears the buffer by reseting the begin and end iterators.
	      The memory won't be freed by this method.
	  */
	  void clear() 
	  {
		  mpBegin = mpData;
		  mpEnd = mpData;
	  }

	  /** Resizes the buffer to a larger size keeping the current data if is not empty.
	      This method does nothing in case of smaller new size than the current size.
		  @param NewSize The new size of the buffer.
	  */
	  void resize(SizeType NewSize)
	  {
		  if(mCapacity < NewSize)
		  {
			  mCapacity=NewSize;

			  SizeType begin_offset = mpBegin - mpData;
			  SizeType end_offset = mpEnd - mpData;
			  if(mpData)
				  mpData = (BlockType*)realloc(mpData, mCapacity * sizeof(BlockType));
			  else
				  mpData = (BlockType*)malloc(mCapacity * sizeof(BlockType));

			  mpBegin = mpData + begin_offset;
			  mpEnd = mpData + end_offset;
		  }
	  }

    
      ///@}
      ///@name Access
      ///@{ 
      
	  /// Returns an iterator to the beginning of the buffer.
       iterator                   begin()            { return mpBegin; }
	  /// Returns a constant iterator to the beginning of the buffer.
       const_iterator             begin() const      { return mpBegin; }
	
	  /// Returns an iterator to the end of the buffer.
       iterator                   end()              { return mpEnd; }
	  /// Returns a constant iterator to the end of the buffer.
       const_iterator             end()   const      { return mpEnd; }

       ContainerType data()
       {
           return mpData;
       }
	
     
      
      ///@}
      ///@name Inquiry
      ///@{
      
	   /// Returns true if the buffer is empty and false if not.
      bool empty() const {return mpBegin==mpEnd;}
    
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	std::stringstream buffer;
	buffer << "buffer (capacity = " << capacity() << ") : ";
	
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
      
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      ContainerType mpData;
       
      iterator mpBegin;
      
      iterator mpEnd;
      
      SizeType mCapacity;

      ///@} 
      ///@name Private Operators
      ///@{
      
	inline SizeType BlockCompatibleSize(SizeType rSize)
	{
	  const SizeType block_size = sizeof(BlockType);
	  return static_cast<SizeType>(((block_size - 1) + rSize) / block_size);	  
	}
        
        
      
        
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


