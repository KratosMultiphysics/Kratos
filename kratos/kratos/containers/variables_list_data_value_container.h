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
//   Last Modified by:    $Author: nelson $
//   Date:                $Date: 2008-12-09 15:23:36 $
//   Revision:            $Revision: 1.10 $
//
//


#if !defined(KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED 



// System includes
#include <string>
#include <iostream> 
#include <cstddef>
#include <cstring>
 



// External includes 


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/variables_list.h"
#include "includes/global_variables.h"


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
  class VariablesListDataValueContainer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VariablesListDataValueContainer
      KRATOS_CLASS_POINTER_DEFINITION(VariablesListDataValueContainer);

      typedef VariablesList::BlockType BlockType;
  
      /// Type of the container used for variables 
      typedef BlockType* ContainerType;

      typedef std::size_t IndexType;

      typedef std::size_t SizeType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      VariablesListDataValueContainer(SizeType NewQueueSize = 1) 
	  : mQueueSize(NewQueueSize), mpCurrentPosition(0), 
	    mpData(0), mpVariablesList(&Globals::DefaultVariablesList)
	  {
	      // Allcating memory
	      Allocate();

	      // Setting the current position at the begining of data
	      mpCurrentPosition = mpData;
	      
	      SizeType size = mpVariablesList->DataSize();
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  BlockType* position = Position(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		      i_variable->AssignZero(position + i * size);
	      }
	  }

      /// Copy constructor.
      VariablesListDataValueContainer(VariablesListDataValueContainer const& rOther) 
	  : mQueueSize(rOther.mQueueSize), mpCurrentPosition(0),
	    mpData(0), mpVariablesList(rOther.mpVariablesList)
	  {
	      // Allcating memory
	      Allocate();

	      // Setting the current position with relative source container offset
	      mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);
  
	      SizeType size = mpVariablesList->DataSize();
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  SizeType offset = LocalOffset(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		  {
		      SizeType total_offset =  offset + i * size;
		      i_variable->Copy(rOther.mpData + total_offset, mpData + total_offset);
		  }
	      }
	  }

	  /// Variables list constructor.
	  VariablesListDataValueContainer(VariablesList*  pVariablesList, SizeType NewQueueSize = 1) 
	      : mQueueSize(NewQueueSize), mpCurrentPosition(0),
	        mpData(0), mpVariablesList(pVariablesList)
	  {
	      // Allcating memory
	      Allocate();

	      // Setting the current position at the begining of data
	      mpCurrentPosition = mpData;
	      
	      SizeType size = mpVariablesList->DataSize();
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  BlockType*  position = Position(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		      i_variable->AssignZero(position + i * size);
	      }
	  }

	  /// Variables list and data constructor
	  VariablesListDataValueContainer(VariablesList*  pVariablesList, BlockType const * ThisData, SizeType NewQueueSize = 1) 
	      : mQueueSize(NewQueueSize), mpCurrentPosition(0),
	        mpData(0), mpVariablesList(pVariablesList)
	  {
	      // Allcating memory
	      Allocate();

	      // Setting the current position at the begining of data
	      mpCurrentPosition = mpData;
	      
	      SizeType size = mpVariablesList->DataSize();
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  SizeType offset = LocalOffset(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		  {
		      SizeType total_offset =  offset + i * size;
		      i_variable->Copy(ThisData + total_offset, mpData + total_offset);
		  }
	      }
	  }

	  /// Destructor.
	  virtual ~VariablesListDataValueContainer()
	  {
		  Clear();
	  }
      

      ///@}
      ///@name Operators 
      ///@{
      
//       template<class TDataType> const TDataType& operator()(const VariableData& rThisVariable, SizeType QueueIndex = 0) const
// 	{
// 	  return GetValue(rThisVariable, QueueIndex);
// 	}
      
      template<class TDataType> 
      TDataType& operator()(const Variable<TDataType>& rThisVariable)
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TDataType>
      TDataType& operator()(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
	{
	  return GetValue(rThisVariable, QueueIndex);
	}
      
      template<class TDataType> 
      const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TDataType> 
      const TDataType& operator()(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
	{
	  return GetValue(rThisVariable, QueueIndex);
	}
      
      template<class TAdaptorType> 
      typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> 
      typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	}
      
      template<class TAdaptorType> 
      const typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> 
      const typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	}
      
//       template<class TDataType> TDataType& operator[](const VariableData& rThisVariable)
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}
      
//       template<class TDataType> const TDataType& operator[](const VariableData& rThisVariable) const
// 	{
// 	  return GetValue(rThisVariable, 0);
// 	}
      
      template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
	{
	  return GetValue(rThisVariable);
	}
      
      template<class TAdaptorType> typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), 0));
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), 0));
	}
      
      /// Assignment operator.
      VariablesListDataValueContainer& operator=(const VariablesListDataValueContainer& rOther)
	{
	    if(rOther.mpVariablesList == 0)
		Clear();
	    else if((mpVariablesList == rOther.mpVariablesList) && 
		    (mQueueSize == rOther.mQueueSize))
	    {
		mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);  

		SizeType size = mpVariablesList->DataSize();
		for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		    i_variable != mpVariablesList->end() ; i_variable++)
		{
		    SizeType offset = LocalOffset(*i_variable);
		    for(SizeType i = 0 ; i < mQueueSize ; i++)
		    {
			SizeType total_offset =  offset + i * size;
			i_variable->Assign(rOther.mpData + total_offset, mpData + total_offset);
		    }
		}
	    }
	  else
	  {
	      DestructAllElements();

	      mQueueSize = rOther.mQueueSize;
	      mpVariablesList = rOther.mpVariablesList;
	      
	      Reallocate();
	    
	      mpCurrentPosition = mpData + (rOther.mpCurrentPosition - rOther.mpData);  

	      SizeType size = mpVariablesList->DataSize();
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  SizeType offset = LocalOffset(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		    {
			SizeType total_offset =  offset + i * size;
			i_variable->Copy(rOther.mpData + total_offset, mpData + total_offset);
		    }
	      }
	  }
	  
	  return *this;
	}
      
      ///@}
      ///@name Operations
      ///@{
      
	template<class TDataType> 
	TDataType& GetValue(const Variable<TDataType>& rThisVariable)
	    {
		if(!mpVariablesList->Has(rThisVariable))
		    KRATOS_ERROR(std::invalid_argument, "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:",rThisVariable);
		return *(TDataType*)Position(rThisVariable);
	    }
      
	template<class TDataType> 
	TDataType& GetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
	    {
		if(!mpVariablesList->Has(rThisVariable))
		    KRATOS_ERROR(std::invalid_argument, "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:",rThisVariable);
		return *(TDataType*)Position(rThisVariable, QueueIndex);
	    }

	template<class TDataType> 
	const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
	    {
		if(!mpVariablesList->Has(rThisVariable))
		    KRATOS_ERROR(std::invalid_argument, "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:",rThisVariable);
		return *(const TDataType*)Position(rThisVariable);
	    }

	template<class TDataType> 
	const TDataType& GetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
	    {
		if(!mpVariablesList->Has(rThisVariable))
		    KRATOS_ERROR(std::invalid_argument, "This container only can store the variables specified in its variables list. The variables list doesn't have this variable:",rThisVariable);
		return *(const TDataType*)Position(rThisVariable, QueueIndex);
	    }
	

	template<class TAdaptorType> 
	typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable)
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	    }

	template<class TAdaptorType> 
	typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex)
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	    }
	
	template<class TAdaptorType> 
	const typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable) const
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	    }
	
	template<class TAdaptorType> 
	const typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex) const
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	    }
	
	template<class TDataType> 
	TDataType& FastGetValue(const Variable<TDataType>& rThisVariable)
	    {
		return *(TDataType*)Position(rThisVariable);
	    }
	
	template<class TDataType> 
	TDataType* pFastGetValue(const Variable<TDataType>& rThisVariable)
	    {
		return (TDataType*)Position(rThisVariable);
	    }
	
	template<class TDataType> 
	TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex)
	    {
		return *(TDataType*)Position(rThisVariable, QueueIndex);
	    }
	
	template<class TDataType> 
	TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex, SizeType ThisPosition)
	    {
		return *(TDataType*)(Position(QueueIndex) + ThisPosition);
	    }
	
	template<class TDataType> 
	TDataType& FastGetCurrentValue(const Variable<TDataType>& rThisVariable, SizeType ThisPosition)
	    {
		return *(TDataType*)(mpCurrentPosition + ThisPosition);
	    }
	
	template<class TDataType>
	const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable) const
	    {
		return *(const TDataType*)Position(rThisVariable);
	    }
	
	template<class TDataType>
	const TDataType* pFastGetValue(const Variable<TDataType>& rThisVariable) const
	    {
		return (const TDataType*)Position(rThisVariable);
	    }

	template<class TDataType>
	const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex) const
	    {
		return *(const TDataType*)Position(rThisVariable, QueueIndex);
	    }

	template<class TDataType> 
	const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable, SizeType QueueIndex, SizeType ThisPosition) const
	    {
		return *(TDataType*)(Position(QueueIndex) + ThisPosition);
	    }
	
	
	template<class TDataType> 
	const TDataType& FastGetCurrentValue(const Variable<TDataType>& rThisVariable, SizeType ThisPosition) const
	    {
		return *(TDataType*)(mpCurrentPosition + ThisPosition);
	    }

	template<class TAdaptorType>
	typename TAdaptorType::Type& FastGetValue(const VariableComponent<TAdaptorType>& rThisVariable)
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	    }
      
	template<class TAdaptorType>
	typename TAdaptorType::Type& FastGetValue(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex)
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	    }
      
      
      template<class TAdaptorType>
      const typename TAdaptorType::Type& FastGetValue(const VariableComponent<TAdaptorType>& rThisVariable) const
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	    }
	
      template<class TAdaptorType>
      const typename TAdaptorType::Type& FastGetValue(const VariableComponent<TAdaptorType>& rThisVariable, SizeType QueueIndex) const
	    {
		return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex));
	    }
	

	SizeType Size() const
	{
	    return mpVariablesList->DataSize();
	}
	
      SizeType QueueSize() const
	{
	  return mQueueSize;
	}

      SizeType TotalSize() const
	{
	  return mQueueSize * mpVariablesList->DataSize();
	}

      template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
	{
	  GetValue(rThisVariable) = rValue;
	}

      template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue, SizeType QueueIndex)
	{
	  GetValue(rThisVariable, QueueIndex) = rValue;
	}

      template<class TAdaptorType> void SetValue(const VariableComponent<TAdaptorType>& rThisVariable, typename TAdaptorType::Type const& rValue)
	{
	  rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable())) = rValue;
	}

      template<class TAdaptorType> void SetValue(const VariableComponent<TAdaptorType>& rThisVariable, typename TAdaptorType::Type const& rValue, SizeType QueueIndex)
	{
	  rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable(), QueueIndex)) = rValue;
	}
	
//       template<class TDataType> void Erase(const Variable<TDataType>& rThisVariable)
// 	{
// 	  typename ContainerType::iterator i;
      
// 	  if ((i = std::find_if(mpData.begin(), mpData.end(), IndexCheck(rThisVariable.Key())))  != mpData.end())
// 	    {
// 	      i->first->Delete(i->second);
// 	      mpData.erase(i);
// 	    }
// 	}

      void Clear()
	  {
	      DestructAllElements();
	      if(mpData)
		  free(mpData);

	      mpData = 0;
	  }


	///@}
      ///@name Access
      ///@{ 

      VariablesList* pGetVariablesList()
	  {
		  return mpVariablesList;
	  }

      const VariablesList * pGetVariablesList() const
	  {
		  return mpVariablesList;
	  }



	void SetVariablesList(VariablesList* pVariablesList)
	{
	    DestructAllElements();
	    
	    mpVariablesList = pVariablesList;
	    
	    Reallocate();

	    mpCurrentPosition = mpData;

	    SizeType size = mpVariablesList->DataSize();
	    for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		i_variable != mpVariablesList->end() ; i_variable++)
	    {
		BlockType*  position = Position(*i_variable);
		for(SizeType i = 0 ; i < mQueueSize ; i++)
		    i_variable->AssignZero(position + i * size);
	    }
	}

	void SetVariablesList(VariablesList* pVariablesList, SizeType ThisQueueSize)
	{
	    DestructAllElements();
	    
	    mpVariablesList = pVariablesList;

	    mQueueSize = ThisQueueSize;

	    Reallocate();

	    mpCurrentPosition = mpData;
	    
	    SizeType size = mpVariablesList->DataSize();
	    for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		i_variable != mpVariablesList->end() ; i_variable++)
	    {
		BlockType*  position = Position(*i_variable);
		for(SizeType i = 0 ; i < mQueueSize ; i++)
		    i_variable->AssignZero(position + i * size);
	    }
	}

	VariablesList& GetDefaultVariablesList()
	  {
		  return Globals::DefaultVariablesList;
	  }

	void Resize(SizeType NewSize)
	{
	    if(mQueueSize == NewSize) // if new size is equal to the provious one
		return; // Do nothing!!
	    
	    if(mQueueSize > NewSize) // if new size is smaller
	    {
		// Destructing elements out of the new size
		for(SizeType i = NewSize ; i < mQueueSize ; i++)
		    DestructElements(i);

		SizeType size = mpVariablesList->DataSize();
			
		//allocating memory
		BlockType* temp = (BlockType*)malloc(size * sizeof(BlockType) * NewSize);

		//Copying data to allocated memory
		for(SizeType i = 0 ; i < NewSize ; i++)
		    memcpy(temp + i * size, Position(i), size * sizeof(BlockType));

		// Updating the queue size
		mQueueSize = NewSize;

		// freeing the old memory
		free(mpData);

		// Setting data pointer to the allocated memory
		mpData = temp;

		// updating the current position
		mpCurrentPosition = mpData;
		
	    }
	    else
	    {
		// Calculating the difference of new and old sizes
		SizeType difference = NewSize - mQueueSize;

		//keeping the old queue size
		SizeType old_size = mQueueSize;

		// Getting the relative offset of current position
		SizeType current_offset = mpCurrentPosition - mpData;

		//Updating the queue size
		mQueueSize = NewSize;

		// Reallocating for new size
		Reallocate();

		SizeType size = mpVariablesList->DataSize();

		// Updating the mpCurrentPosition
		mpCurrentPosition = mpData + current_offset; 

		// moving the region after current position to the end
		SizeType region_size = old_size * size - current_offset;
		memmove(mpCurrentPosition + difference * size, mpCurrentPosition, region_size * sizeof(BlockType));

		// Assigning zero to the added elements
		for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		    i_variable != mpVariablesList->end() ; i_variable++)
		{
		    BlockType*  position = mpCurrentPosition + LocalOffset(*i_variable);
		    for(SizeType i = 0 ; i < difference ; i++)
			i_variable->AssignZero(position + i * size);
		}

		//moving the current position to the moved place
		mpCurrentPosition +=  difference * size;
	    }
	}


// 	void Resize(SizeType NewSize, BlockType* SourceData)
// 	{
// 	    if(mQueueSize == NewSize)
// 		return;
	    
// 	    if(mQueueSize > NewSize)
// 	    {
// 		for(SizeType i = NewSize ; i < mQueueSize ; i++)
// 		    DestructElements(i);
			
// 		mQueueSize = NewSize;

// 		Reallocate();
// 	    }
// 	    else
// 	    {
// 		SizeType old_size = mQueueSize;
// 		mQueueSize = NewSize;
// 		Reallocate();
// 		SizeType size = mpVariablesList->DataSize();
// 		for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
// 		    i_variable != mpVariablesList->end() ; i_variable++)
// 		{
// 		    SizeType  offset = LocalOffset(*i_variable);
// 		    for(SizeType i = old_size ; i < mQueueSize ; i++)
// 			i_variable->Copy(SourceData + offset, mpData + offset + i * size);
// 		}
// 	    }
// 	}


	BlockType* Data()
	    {
		return mpData;
	    }

	BlockType* Data(SizeType QueueIndex)
	    {
		return Position(QueueIndex);
	    }

	SizeType DataSize()
	    {
		return mpVariablesList->DataSize() * sizeof(BlockType);
	    }

	SizeType TotalDataSize()
	    {
		return mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize;
	    }

	void AssignData(BlockType* Source, SizeType QueueIndex)
	    {
		AssignData(Source, Position(QueueIndex));
	    }

//	void SetData(BlockType* pThisData)
//	    {
//		if(mpData)
//			free(mpData);
//		mpData = pThisData;
//	    }

	void CloneFront()
	{
	    if(mQueueSize == 0)
	    {
		Resize(1);
		return;
	    }
	    
	    if(mQueueSize == 1)
		return;
	    
	    SizeType size = mpVariablesList->DataSize();
	    BlockType* position = (mpCurrentPosition == mpData) ? mpData + TotalSize() - size :  mpCurrentPosition - size;
	    AssignData(mpCurrentPosition, position);
	    mpCurrentPosition = position;
	    
	}
	
	void PushFront()
	{
	    if(mQueueSize == 0)
	    {
		Resize(1);
		return;
	    }
	    
	    if(mQueueSize == 1)
		return;
	    
	    SizeType size = mpVariablesList->DataSize();
	    mpCurrentPosition = (mpCurrentPosition == mpData) ? mpData + TotalSize() - size :  mpCurrentPosition - size;
	    AssignZero();

	}
      

	void AssignZero()
	{
	    for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		i_variable != mpVariablesList->end() ; i_variable++)
		i_variable->AssignZero(mpCurrentPosition + LocalOffset(*i_variable));
	}

	void AssignZero(SizeType QueueIndex)
	{
	    BlockType* position = Position(QueueIndex);
	    for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		i_variable != mpVariablesList->end() ; i_variable++)
		i_variable->AssignZero(position + LocalOffset(*i_variable));
	}
	
      ///@}
      ///@name Inquiry
      ///@{
      
      template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
	{
	  return mpVariablesList->Has(rThisVariable);
	}

      template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return mpVariablesList->Has(rThisVariable.GetSourceVariable());
	}
      
      bool IsEmpty()
	{
	  return mpVariablesList->IsEmpty();
	}
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return std::string("variables list data value container");
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "variables list data value container";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	  for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
	      i_variable != mpVariablesList->end() ; i_variable++)
	  {
	      rOStream <<"    ";
	      for(SizeType i = 0 ; i < mQueueSize ; i++)
	      {
		  rOStream << i << ": ";
		  i_variable->Print(Position(*i_variable, i), rOStream);
		  rOStream << "  ";
	      }
	      rOStream << std::endl;
	  }
	  
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
	  
      SizeType mQueueSize;
      
      BlockType* mpCurrentPosition;

      ContainerType mpData;
      
      VariablesList* mpVariablesList;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
      inline void Allocate()
	  {
	      mpData = (BlockType*)malloc(mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
	  }

	inline void Reallocate()
	    {
		mpData = (BlockType*)realloc(mpData, mpVariablesList->DataSize() * sizeof(BlockType) * mQueueSize);
	    }
        
      void DestructElements(SizeType ThisIndex)
	  {
	      if(mpData == 0)
		  return;
	      BlockType* position = Position(ThisIndex);
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
		  i_variable->Destruct(position + LocalOffset(*i_variable));
	  }
      
      
      void DestructAllElements()
	  {
	      if(mpData == 0)
		  return;
	      for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		  i_variable != mpVariablesList->end() ; i_variable++)
	      {
		  SizeType size = mpVariablesList->DataSize();
		  BlockType*  position = mpData + LocalOffset(*i_variable);
		  for(SizeType i = 0 ; i < mQueueSize ; i++)
		      i_variable->Destruct(position + i * size);
	      }
	  }
      

	void AssignData(BlockType* Source, BlockType* Destination)
	    {
		for(VariablesList::const_iterator i_variable = mpVariablesList->begin() ;
		    i_variable != mpVariablesList->end() ; i_variable++)
		{
		    SizeType offset = LocalOffset(*i_variable);
		    i_variable->Assign(Source + offset, Destination + offset);
		}
	    }
        
      inline SizeType LocalOffset(VariableData const & rThisVariable) const
	  {
	      return mpVariablesList->Index(rThisVariable);
	  }

      inline BlockType* Position(VariableData const & rThisVariable) const
	  {
	      return mpCurrentPosition + mpVariablesList->Index(rThisVariable);
	  }

      inline BlockType* Position(VariableData const & rThisVariable, SizeType ThisIndex) const
	  {
	      return Position(ThisIndex) + mpVariablesList->Index(rThisVariable);
	  }
      
      inline BlockType* Position() const
	  {
	      return mpCurrentPosition;
	  }

      inline BlockType* Position(SizeType ThisIndex) const
	  {
	      SizeType total_size = TotalSize();
	      BlockType* position = mpCurrentPosition + ThisIndex * mpVariablesList->DataSize();
	      return (position < mpData + total_size) ? position : position - total_size;
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
        
    }; // Class VariablesListDataValueContainer 

  ///@} 
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    VariablesListDataValueContainer& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const VariablesListDataValueContainer& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VARIABLES_LIST_DATA_VALUE_CONTAINER_H_INCLUDED  defined 


