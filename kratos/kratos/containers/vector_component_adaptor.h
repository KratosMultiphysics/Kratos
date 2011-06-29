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


#if !defined(KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED )
#define  KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "variable.h"


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
  template<class TVectorType>
    class VectorComponentAdaptor
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VectorComponentAdaptor
      KRATOS_CLASS_POINTER_DEFINITION(VectorComponentAdaptor);
  
      typedef typename TVectorType::value_type Type;

      typedef TVectorType SourceType;

      typedef Variable<TVectorType>  SourceVariableType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
      VectorComponentAdaptor(const SourceVariableType& rSourceVariable, int ComponentIndex)
	: mpSourceVariable(&rSourceVariable), mComponentIndex(ComponentIndex)
	{}

      /// Copy constructor.
      VectorComponentAdaptor(const VectorComponentAdaptor& rOther)
	: mpSourceVariable(rOther.mpSourceVariable), mComponentIndex(rOther.mComponentIndex)
	{}

      /// Destructor.
      virtual ~VectorComponentAdaptor(){}
     

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      Type& GetValue(SourceType& rValue) const
      {
	return rValue[mComponentIndex];
      }
      
      const Type& GetValue(const SourceType& rValue) const
      {
	return rValue[mComponentIndex];
      }
      

      static VectorComponentAdaptor const& StaticObject()
	{
	  return msStaticObject;
	}
      
      ///@}
      ///@name Access
      ///@{ 
      
      const SourceVariableType& GetSourceVariable() const
      {
	return *mpSourceVariable;
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
	  buffer << mpSourceVariable->Name() << " vector component " << mComponentIndex << " adaptor";
	  return buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << mpSourceVariable->Name() << " vector component " << mComponentIndex << " adaptor";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
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
      
      static VectorComponentAdaptor const msStaticObject;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      const SourceVariableType* mpSourceVariable;
      
      int mComponentIndex;
	  
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
      VectorComponentAdaptor();
        
      ///@}    
        
    }; // Class VectorComponentAdaptor 

  ///@} 
  
  template<class TVectorType>
  const VectorComponentAdaptor<TVectorType> VectorComponentAdaptor<TVectorType>::msStaticObject = VectorComponentAdaptor<TVectorType>(VectorComponentAdaptor<TVectorType>::SourceVariableType::StaticObject(), 0);

  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TDataType>
  inline std::istream& operator >> (std::istream& IStream, 
				    VectorComponentAdaptor<TDataType>& rThis);

  /// output stream function
  template<class TDataType>
  inline std::ostream& operator << (std::ostream& OStream, 
				    const VectorComponentAdaptor<TDataType>& rThis)
    {
      rThis.PrintInfo(OStream);
      rThis.PrintData(OStream);

      return OStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VECTOR_COMPONENT_ADAPTOR_H_INCLUDED  defined 


