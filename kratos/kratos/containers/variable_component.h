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


#if !defined(KRATOS_VARIABLE_COMPONENT_H_INCLUDED )
#define  KRATOS_VARIABLE_COMPONENT_H_INCLUDED 



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "variable_data.h"


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
  
  /// Provide information for store or retrive a component of a variable in data container.
  /** Provide information for store or retrive a component of a
      variable in data container. This class also provide a method to
      extract its component value from the source variable value in
      container.
  */
  template<class TAdaptorType>
  class VariableComponent : public VariableData
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VariableComponent
      KRATOS_CLASS_POINTER_DEFINITION(VariableComponent);
  
      /// Type of this component.
      typedef typename TAdaptorType::Type Type;

      /// Data value type of this component.
      typedef typename TAdaptorType::Type DataType;

      /// Source value type of this component.
      typedef typename TAdaptorType::SourceType SourceType;

      /// Source variable type of this component.
      typedef Variable<typename TAdaptorType::SourceType> SourceVariableType;

      /// Base class type definition.
      typedef VariableData BaseType;

      /// Adaptor type.
      typedef TAdaptorType AdaptorType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      VariableComponent(const std::string& NewName, const AdaptorType& NewAdaptor) 
	: BaseType(NewName, sizeof(DataType)), mAdaptor(NewAdaptor)
	{
	}

      /// Copy constructor.
      VariableComponent(const VariableComponent& rOther)
	: BaseType(rOther), mAdaptor(rOther.mAdaptor){}
      
      /// Destructor.
      virtual ~VariableComponent(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      const SourceVariableType& GetSourceVariable() const
	{
	  return mAdaptor.GetSourceVariable();
	}
      
      const AdaptorType& GetAdaptor() const
	{
	  return mAdaptor;
	}

      DataType& GetValue(SourceType& SourceValue) const
	{
	  return mAdaptor.GetValue(SourceValue);
	}
      
      const DataType& GetValue(const SourceType& SourceValue) const
	{
	  return mAdaptor.GetValue(SourceValue);
	}

      static VariableComponent const& StaticObject()
	{
	  return msStaticObject;
	}

      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << Name() << " component of " <<  mAdaptor.GetSourceVariable().Name() << " variable";
	  return buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Name() << " component of " <<  mAdaptor.GetSourceVariable().Name() << " variable";
	}

      /// Print object's data.
//       virtual void PrintData(std::ostream& rOStream) const;
      
            
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
        
      /// Assignment operator.
      VariableComponent& operator=(const VariableComponent& rOther)
	{
	  BaseType::operator=(rOther);
	  mAdaptor = rOther.mAdaptor;
	}
        
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

      static const VariableComponent  msStaticObject;
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      TAdaptorType mAdaptor;
	  
      ///@} 
      ///@name Serialization
      ///@{ 
        
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
      
      /// Default constructor.
      VariableComponent(){}
        
      ///@}    
        
    }; // Class VariableComponent 

  ///@} 

  template<class TAdaptorType>
  const VariableComponent<TAdaptorType> VariableComponent<TAdaptorType>::msStaticObject("NONE", TAdaptorType::StaticObject());
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TAdaptorType>
  inline std::istream& operator >> (std::istream& IStream, 
				    VariableComponent<TAdaptorType>& rThis);

  /// output stream function
  template<class TAdaptorType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const VariableComponent<TAdaptorType>& rThis)
    {
      rThis.PrintInfo(OStream);
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined 


