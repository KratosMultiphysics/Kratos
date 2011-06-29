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


#if !defined(KRATOS_VARIABLE_H_INCLUDED )
#define  KRATOS_VARIABLE_H_INCLUDED



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
  
  /// Variable class contains all information needed to store and retrive data from a data container.
  /** Variable class contains all information needed to store and
      retrive data from a data container.  It contains key value which
      is needed for searching in data container. Also a zero value to
      use as a default vlaue in container. Finally it has the type of
      the Variable as its template parameter so container can find it
      int its relative part.
  */
  template<class TDataType>
    class Variable : public VariableData
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Variable
      KRATOS_CLASS_POINTER_DEFINITION(Variable);

      /// type of this variable 
      typedef TDataType Type;

      // Type used for key values which defined in VariableData
      typedef VariableData::KeyType KeyType;

      // Type of this varible with given TDataType
      typedef Variable<TDataType> VariableType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /** Constructor with specific name and zero value */
      Variable(const std::string& NewName, const TDataType Zero = TDataType()) 
	: VariableData(NewName, sizeof(TDataType)), mZero(Zero)
	{
	}
        
      /// Copy constructor.
      Variable(const VariableType& rOtherVariable) : VariableData(rOtherVariable), mZero(rOtherVariable.mZero){}

      /// Destructor.
      virtual ~Variable(){}

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      void* Clone(const void* pSource) const
	{
	  return new TDataType(*static_cast<const TDataType* >(pSource) );
	}
      
      void* Copy(const void* pSource, void* pDestination) const 
	{
	  return new(pDestination) TDataType(*static_cast<const TDataType* >(pSource) );
	}
      
      void Assign(const void* pSource, void* pDestination) const 
	{
	  (*static_cast<TDataType* >(pDestination) ) = (*static_cast<const TDataType* >(pSource) );
	}
      
      void AssignZero(void* pDestination) const 
	{
	  //(*static_cast<TDataType* >(pDestination) ) = mZero;
	  new (pDestination) TDataType(mZero);
	}
      
      void Delete(void* pSource) const
	{
	  delete static_cast<TDataType* >(pSource);
	}
      
      void Destruct(void* pSource) const
	{
	   static_cast<TDataType* >(pSource)->~TDataType();
	}
	
	void Print(const void* pSource, std::ostream& rOStream) const
	{
	  rOStream << Name() << " : " << *static_cast<const TDataType* >(pSource) ;
	}
	
	virtual void Save(Serializer& rSerializer, void* pData) const
	{
	  // I'm saving by the value, it can be done by the pointer to detect shared data. Pooyan.
	  rSerializer.save("Data",*static_cast<TDataType* >(pData));
	}
	
	virtual void Allocate(void** pData) const
	{
	  *pData = new TDataType;
	}
	
	virtual void Load(Serializer& rSerializer, void* pData) const
	{
	  rSerializer.load("Data",*static_cast<TDataType* >(pData));
	}
	
      static const VariableType& StaticObject(){return msStaticObject;}
        
      ///@}
      ///@name Access
      ///@{ 
      
      const TDataType& Zero() const
	{
	  return mZero;
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
	  buffer << Name() << " variable";
	  return buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Name() << " variable";
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
      VariableType& operator=(const VariableType& rOtherVariable)
	{
	  VariableData::operator=(rOtherVariable);
	  mZero = rOtherVariable.mZero;
	  return *this;
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
        
      static const VariableType msStaticObject;
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      TDataType mZero;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
    
      ///@} 
      ///@name Serialization
      ///@{ 

	friend class Serializer;
	
	virtual void save(Serializer& rSerializer) const
	{
	  KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, VariableData );
	  rSerializer.save("Zero",mZero);
	}
        
	virtual void load(Serializer& rSerializer)
	{
	  KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, VariableData );
	  rSerializer.load("Zero",mZero);
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
      
      /** Default constructor is un accessible due to the fact that
	each variable must have a name defined.*/
      Variable(){}
        
      ///@}    
        
    }; // Class Variable 

  ///@} 
  
  template<class TDataType>
  const Variable<TDataType> Variable<TDataType>::msStaticObject("NONE");

  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TDataType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Variable<TDataType>& rThis);

  /// output stream function
  template<class TDataType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Variable<TDataType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VARIABLE_H_INCLUDED  defined 


