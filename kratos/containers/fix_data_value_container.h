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
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_FIX_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_FIX_DATA_VALUE_CONTAINER_H_INCLUDED 



// System includes
#include <string>
#include <iostream> 
#include <cstddef> 


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
  class FixDataValueContainer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of FixDataValueContainer
      KRATOS_CLASS_POINTER_DEFINITION(FixDataValueContainer);
  
      /// Type of the container used for variables 
      typedef double* ContainerType;

      typedef std::size_t IndexType;

      typedef std::size_t SizeType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      FixDataValueContainer() : mSize(0), mpData(0), mpVariablesList(&Globals::DefaultVariablesList){}

      /// Copy constructor.
      FixDataValueContainer(FixDataValueContainer const& rOther) : mSize(rOther.Size()), mpData(0), mpVariablesList(rOther.mpVariablesList)
	{
		if(mSize == 0) 
			return;
	  mpData = (double*)malloc(mSize * sizeof(double));

	  VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	  SizeType size = variables.size();

	  for(IndexType i = 0 ; i < size ; ++i)
	    {
	      const VariableData* p_variable = variables[i];
	      std::size_t offset = mpVariablesList->Index(p_variable);
	      if(offset < rOther.mSize)
		p_variable->Copy(rOther.mpData + offset, mpData + offset);
	      else
		p_variable->AssignZero(mpData + offset);
	    }
	}

      /// Variables list constructor.
      FixDataValueContainer(VariablesList*  pVariablesList) : mSize(0), mpData(0), mpVariablesList(pVariablesList)
	{
	  mpVariablesList = pVariablesList;

	  mSize = mpVariablesList->DataSize();

	  mpData = (double*)malloc(mSize * sizeof(double));

	  VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	  SizeType size = variables.size();

	  for(IndexType i = 0 ; i < size ; ++i)
	    {
	      const VariableData* p_variable = variables[i];
	      std::size_t offset = mpVariablesList->Index(p_variable);
	      p_variable->AssignZero(mpData + offset);
	    }
		mSize = variables.size();
	}

      /// Destructor.
      virtual ~FixDataValueContainer()
      {

	  VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	  SizeType size = variables.size();

	  for(IndexType i = 0 ; i < size ; ++i)
	    {
	      const VariableData* p_variable = variables[i];
	      std::size_t offset = mpVariablesList->Index(p_variable);
		  if(offset < mSize)
	      p_variable->Destruct(mpData + offset);
	    }

	  if(mpData)
	    free(mpData);
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
      template<class TDataType> const TDataType& operator()(const VariableData& rThisVariable) const
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TDataType> TDataType& operator()(const Variable<TDataType>& rThisVariable)
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TDataType> const TDataType& operator()(const Variable<TDataType>& rThisVariable) const
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TAdaptorType> typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& operator()(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TDataType> TDataType& operator[](const VariableData& rThisVariable)
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TDataType> const TDataType& operator[](const VariableData& rThisVariable) const
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TDataType> TDataType& operator[](const Variable<TDataType>& rThisVariable)
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TDataType> const TDataType& operator[](const Variable<TDataType>& rThisVariable) const
	{
	  return GetValue<TDataType>(rThisVariable);
	}
      
      template<class TAdaptorType> typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& operator[](const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      /// Assignment operator.
      FixDataValueContainer& operator=(const FixDataValueContainer& rOther)
	{
	  if(rOther.mpVariablesList == 0)
	    {
	      //TODO: empty this container
	    }
	    
	  if(mpVariablesList == rOther.mpVariablesList)
	    {
// 	      std::cout << "using same variables list" << std::endl;
	      VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	      SizeType size = variables.size();
	      
// 	      std::cout << size << " ?= " << mSize << std::endl;
	      if(mpVariablesList->DataSize() != mSize)
		{
		  mSize = mpVariablesList->DataSize();
		  mpData = (double*)realloc(mpData, mSize * sizeof(double));
		}

	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		  
		  if(offset < rOther.mSize)
		    p_variable->Assign(rOther.mpData + offset, mpData + offset);
		  else
		    p_variable->AssignZero(mpData + offset);
		}
	    }
	  else
	    {
	      VariablesList::VariablesContainerType const& old_variables = mpVariablesList->Variables();

	      SizeType size = old_variables.size();

	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = old_variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		if(offset < mSize)
		  p_variable->Destruct(mpData + offset);
		}

	      mpVariablesList = rOther.mpVariablesList;

	      mSize = mpVariablesList->DataSize();

	      mpData = (double*)realloc(mpData, mSize * sizeof(double));
	      
	      VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	      size = variables.size();
	      
	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		  if(offset < rOther.mSize)
		    p_variable->Copy(rOther.mpData + offset, mpData + offset);
		  else
		    p_variable->AssignZero(mpData + offset);
		}
	      
	    }
	  
	  return *this;
	}
      
      ///@}
      ///@name Operations
      ///@{
      
      template<class TDataType> TDataType& GetValue(const Variable<TDataType>& rThisVariable)
	{
	  if(!mpVariablesList->Has(rThisVariable))
	    mpVariablesList->Add(rThisVariable);

//  	  KRATOS_WATCH(rThisVariable)

	  IndexType index = mpVariablesList->Index(rThisVariable);


	  if((index >= mSize)||(mSize == 0))
	    Update();

//  	  KRATOS_WATCH(*(TDataType*)(mpData + index))
	  return *(TDataType*)(mpData + index);
	}

      template<class TDataType> const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
	{
	  if(!mpVariablesList->Has(rThisVariable))
	    return rThisVariable.Zero();
	  //std:: cout << "oh oh" << std::endl;
 	 // KRATOS_WATCH(rThisVariable)
	  IndexType index = mpVariablesList->Index(rThisVariable);

//  	  KRATOS_WATCH(index)
//  	  KRATOS_WATCH(mSize)
	  if(index >= mSize)
	    return rThisVariable.Zero();

	  return *(const TDataType*)(mpData + index);
	}

	//*******************************************************************************************
	//by Riccardo: variables are not added
     template<class TDataType> TDataType& FastGetValue(const Variable<TDataType>& rThisVariable)
	{
	  IndexType index = mpVariablesList->Index(rThisVariable);
#ifdef _DEBUG
	  //KRATOS_WATCH("attention printing FastGetValue");
	  if(!mpVariablesList->Has(rThisVariable))
		  KRATOS_ERROR(std::logic_error,"","");
	  if(index >= mSize)
	    KRATOS_ERROR(std::logic_error,"","");
#endif
	  return *(TDataType*)(mpData + index);
	}
      template<class TDataType> const TDataType& FastGetValue(const Variable<TDataType>& rThisVariable) const
	{
	  IndexType index = mpVariablesList->Index(rThisVariable);
#ifdef _DEBUG
	  //KRATOS_WATCH("attention printing FastGetValue");
	  if(!mpVariablesList->Has(rThisVariable))
		  KRATOS_ERROR(std::logic_error,"","");
	  if(index >= mSize)
	    KRATOS_ERROR(std::logic_error,"","");
#endif
	  return *(const TDataType*)(mpData + index);
	}
	//*******************************************************************************************

      template<class TAdaptorType> typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}

      SizeType Size() const
	{
	  return mpVariablesList->DataSize();
	}

      template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
	{
	  GetValue(rThisVariable) = rValue;
	}

      template<class TAdaptorType> void SetValue(const VariableComponent<TAdaptorType>& rThisVariable, typename TAdaptorType::Type const& rValue)
	{
	  rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable())) = rValue;
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
	  VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	  SizeType size = variables.size();

	  for(IndexType i = 0 ; i < size ; ++i)
	    {
	      const VariableData* p_variable = variables[i];
	      std::size_t offset = mpVariablesList->Index(p_variable);
		  if(offset < mSize)
	      p_variable->Destruct(mpData + offset);
	    }
	  if(mpData)
	    free(mpData);
	}

      void Update()
	{
	  SizeType old_size = mSize;

	      mSize = mpVariablesList->DataSize();

	      if(mpData)
		mpData = (double*)realloc(mpData, mSize * sizeof(double));
	      else
		mpData = (double*)malloc(mSize * sizeof(double));

	      
	      VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	      SizeType size = variables.size();
	      
	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		  if(offset >= old_size)
		    p_variable->AssignZero(mpData + offset);
		}
	}
      
      void UpdateVariablesList(VariablesList* pVariablesList)
	{
	      VariablesList::VariablesContainerType const& old_variables = mpVariablesList->Variables();

	      SizeType size = old_variables.size();

	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = old_variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		  if(offset < mSize)
		  p_variable->Destruct(mpData + offset);
		}

	      mpVariablesList = pVariablesList;
// To see if we need to just free the memory when mSize == 0
//	  if(mpData)
//	    free(mpData); 

	  mpData = (double*)realloc(mpData, mpVariablesList->DataSize() * sizeof(double));
	      
	      VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	      size = variables.size();
	      
	      for(IndexType i = 0 ; i < size ; ++i)
		{
		  const VariableData* p_variable = variables[i];
		  std::size_t offset = mpVariablesList->Index(p_variable);
		  p_variable->AssignZero(mpData + offset);
		}

		//mSize = variables.size();
	      mSize = mpVariablesList->DataSize();
	}
      
      ///@}
      ///@name Access
      ///@{ 

      VariablesList* pGetVariablesList()
	  {
		  return mpVariablesList;
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
	  return std::string("data value container");
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "data value container";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	VariablesList::VariablesContainerType const& variables = mpVariablesList->Variables();

	SizeType size = variables.size();
	      
	for(IndexType i = 0 ; i < size ; ++i)
	  {
	    rOStream <<"    ";
	    const VariableData* p_variable = variables[i];
	    std::size_t offset = mpVariablesList->Index(p_variable);
	    p_variable->Print(mpData + offset, rOStream);
	    rOStream << std::endl;
	  }
// 	for(ConstantIteratorType i = mpData.begin() ; i != mpData.end() ; ++i)
// 	  {
// 	    rOStream <<"    ";
// 	    i->first->Print(i->second, rOStream);
// 	    rOStream << std::endl;
// 	  }
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

      SizeType mSize;
        
      ContainerType mpData;
      
      VariablesList* mpVariablesList;
        
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
        
    }; // Class FixDataValueContainer 

  ///@} 
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    FixDataValueContainer& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const FixDataValueContainer& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FIX_DATA_VALUE_CONTAINER_H_INCLUDED  defined 


