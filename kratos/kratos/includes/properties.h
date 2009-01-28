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
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PROPERTIES_H_INCLUDED )
#define  KRATOS_PROPERTIES_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/data_value_container.h"
#include "includes/process_info.h"


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
  class Properties : public IndexedObject
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Properties
      KRATOS_CLASS_POINTER_DEFINITION(Properties);
  
      typedef IndexedObject BaseType;

      typedef DataValueContainer ContainerType;
      
      typedef Node<3> NodeType;

      typedef NodeType::IndexType IndexType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Properties(IndexType NewId = 0) : BaseType(NewId), mData(){}

      /// Copy constructor.
      Properties(const Properties& rOther) : BaseType(rOther), mData(rOther.mData){}

      /// Destructor.
      virtual ~Properties(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Properties& operator=(const Properties& rOther)
	{
		BaseType::operator=(rOther);
	  mData = rOther.mData;
	  return *this;
	}
      
      template<class TVariableType> 
	typename TVariableType::Type& operator()(const TVariableType& rV)
	{return GetValue(rV);}

      template<class TVariableType> 
	typename TVariableType::Type const& operator()(const TVariableType& rV) const
	{return GetValue(rV);}

      template<class TVariableType> 
	typename TVariableType::Type& operator[](const TVariableType& rV)
	{return GetValue(rV);}

      template<class TVariableType> 
	typename TVariableType::Type const& operator[](const TVariableType& rV) const
	{return GetValue(rV);}

      template<class TVariableType> 
	typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode)
	{return GetValue(rV, rThisNode);}

      template<class TVariableType> 
	typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode) const
	{return GetValue(rV, rThisNode);}

      template<class TVariableType> 
	typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode, IndexType SolutionStepIndex)
	{return GetValue(rV, rThisNode, SolutionStepIndex);}

      template<class TVariableType> 
	typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode, IndexType SolutionStepIndex) const
	{return GetValue(rV, rThisNode, SolutionStepIndex);}

      template<class TVariableType> 
	typename TVariableType::Type& operator()(const TVariableType& rV, NodeType& rThisNode, ProcessInfo const& rCurrentProcessInfo)
	{return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());}

      template<class TVariableType> 
	typename TVariableType::Type const& operator()(const TVariableType& rV, NodeType const& rThisNode, ProcessInfo const& rCurrentProcessInfo) const
	{return GetValue(rV, rThisNode, rCurrentProcessInfo.GetSolutionStepIndex());}

      ///@}
      ///@name Operations
      ///@{
      
      template<class TVariableType> 
	typename TVariableType::Type& GetValue(const TVariableType& rV)
	{
	    return mData.GetValue(rV);
	}
      
      template<class TVariableType> 
	typename TVariableType::Type const& GetValue(const TVariableType& rV) const
	{
	    return mData.GetValue(rV);
	}
      
      template<class TVariableType> 
	typename TVariableType::Type& GetValue(const TVariableType& rV, NodeType& rThisNode)
	{
	  if(mData.Has(rV))
	    return mData.GetValue(rV);

	  return rThisNode.GetValue(rV);
	}
      
      template<class TVariableType> 
	typename TVariableType::Type const& GetValue(const TVariableType& rV, NodeType const& rThisNode) const
	{
	  if(mData.Has(rV))
	    return mData.GetValue(rV);

	  return rThisNode.GetValue(rV);
	}
      
      template<class TVariableType> 
	typename TVariableType::Type& GetValue(const TVariableType& rV, NodeType& rThisNode, IndexType SolutionStepIndex)
	{
	  if(mData.Has(rV))
	    return mData.GetValue(rV);

	  return rThisNode.GetValue(rV, SolutionStepIndex);
	}
      
      template<class TVariableType> 
	typename TVariableType::Type const& GetValue(const TVariableType& rV, NodeType const& rThisNode, IndexType SolutionStepIndex) const
	{
	  if(mData.Has(rV))
	    return mData.GetValue(rV);

	  return rThisNode.GetValue(rV, SolutionStepIndex);
	}
      
      template<class TVariableType> 
	void SetValue(TVariableType const& rV, typename TVariableType::Type const& rValue)
	{
	  mData.GetValue(rV) = rValue;
	}

      ///@}
      ///@name Access
      ///@{
      
      DataValueContainer& Data()
	{
	  return mData;
	}
      
      DataValueContainer const& Data() const
	{
	  return mData;
	}
      
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      template<class TVariableType> 
      bool Has(TVariableType const& rThisVariable) const
	{
	  return mData.Has(rThisVariable);
	}

      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return "Properties";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream <<  "Properties";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  mData.PrintData(rOStream);
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
        
    }; // Class Properties 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Properties& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Properties& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PROPERTIES_H_INCLUDED  defined 


