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


#if !defined(KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED )
#define KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED 



// System includes
#include <string>
#include <iostream> 
#include <cstddef> 
#include <vector> 


// External includes 


// Project includes
#include "includes/define.h"
#include "containers/variable.h"
#include "containers/variable_component.h"


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
  class DataValueContainer
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of DataValueContainer
      KRATOS_CLASS_POINTER_DEFINITION(DataValueContainer);
  
      /// Type of the container used for variables 
      typedef std::pair<const VariableData*, void*> ValueType;

      /// Type of the container used for variables 
      typedef std::vector<ValueType> ContainerType;

      /// Type of the container used for variables 
      typedef std::vector<ValueType>::iterator IteratorType;

      /// Type of the container used for variables 
      typedef std::vector<ValueType>::const_iterator ConstantIteratorType;

      /// Type of the container used for variables 
      typedef std::vector<ValueType>::size_type SizeType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      DataValueContainer() {}

      /// Copy constructor.
      DataValueContainer(DataValueContainer const& rOther)
	{
	  for(ConstantIteratorType i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
 	    mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
	}

      /// Destructor.
      virtual ~DataValueContainer()
      {
	for(IteratorType i = mData.begin() ; i != mData.end() ; ++i)
	    i->first->Delete(i->second);
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
      DataValueContainer& operator=(const DataValueContainer& rOther)
	{
	  Clear();

	  for(ConstantIteratorType i = rOther.mData.begin() ; i != rOther.mData.end() ; ++i)
	    mData.push_back(ValueType(i->first, i->first->Clone(i->second)));
	  
	  return *this;
	}
      
      ///@}
      ///@name Operations
      ///@{
      
      template<class TDataType> TDataType& GetValue(const Variable<TDataType>& rThisVariable)
	{
	  typename ContainerType::iterator i;
      
	  if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key())))  != mData.end())
	    return *static_cast<TDataType*>(i->second);
      
	  mData.push_back(ValueType(&rThisVariable,new TDataType(rThisVariable.Zero())));
      
	  return *static_cast<TDataType*>(mData.back().second);
	}

      template<class TDataType> const TDataType& GetValue(const Variable<TDataType>& rThisVariable) const
	{
	  typename ContainerType::const_iterator i;
      
	  if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key())))  != mData.end())
	    return *static_cast<const TDataType*>(i->second);
      
	  return rThisVariable.Zero();
	}

      template<class TAdaptorType> typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable)
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}
      
      template<class TAdaptorType> const typename TAdaptorType::Type& GetValue(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable()));
	}

      SizeType Size()
	{
	  return mData.size();
	}

      template<class TDataType> void SetValue(const Variable<TDataType>& rThisVariable, TDataType const& rValue)
	{
	  typename ContainerType::iterator i;
      
	  if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key())))  != mData.end())
	    *static_cast<TDataType*>(i->second) = rValue;
	  else
	    mData.push_back(ValueType(&rThisVariable,new TDataType(rValue)));
	}

      template<class TAdaptorType> void SetValue(const VariableComponent<TAdaptorType>& rThisVariable, typename TAdaptorType::Type const& rValue)
	{
	  rThisVariable.GetValue(GetValue(rThisVariable.GetSourceVariable())) = rValue;
	}
      
      template<class TDataType> void Erase(const Variable<TDataType>& rThisVariable)
	{
	  typename ContainerType::iterator i;
      
	  if ((i = std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key())))  != mData.end())
	    {
	      i->first->Delete(i->second);
	      mData.erase(i);
	    }
	}

      void Clear()
	{
	  for(ContainerType::iterator i = mData.begin() ; i != mData.end() ; i++)
	      i->first->Delete(i->second);

	  mData.clear();
	}
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
	{
	  return (std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.Key())) != mData.end());
	}

      template<class TAdaptorType> bool Has(const VariableComponent<TAdaptorType>& rThisVariable) const
	{
	  return (std::find_if(mData.begin(), mData.end(), IndexCheck(rThisVariable.GetSourceVariable().Key())) != mData.end());
	}
      
      bool IsEmpty()
	{
	  return mData.empty();
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
	for(ConstantIteratorType i = mData.begin() ; i != mData.end() ; ++i)
	  {
	    rOStream <<"    ";
	    i->first->Print(i->second, rOStream);
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

      class IndexCheck
	{
		std::size_t mI;
	public:
	  IndexCheck(int I) : mI(I){}
	  bool operator()(const ValueType& I){return I.first->Key() == mI;}
	};

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
        
    }; // Class DataValueContainer 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    DataValueContainer& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const DataValueContainer& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_DATA_VALUE_CONTAINER_H_INCLUDED  defined 


