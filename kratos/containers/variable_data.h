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


#if !defined(KRATOS_VARIABLE_DATA_H_INCLUDED )
#define  KRATOS_VARIABLE_DATA_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/counter.h"


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
  
  /// This class is the base of variables and variable's components which contains their common data.
  /** This class hold variables name and key and also adaptor type for variables components.
      It also has a counter to generate automatically key numbers for new variables.
  */
  class VariableData
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of VariableData
      KRATOS_CLASS_POINTER_DEFINITION(VariableData);

	  typedef std::size_t KeyType;

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Copy constructor
      VariableData(const VariableData& rOtherVariable) 
	: mName(rOtherVariable.mName), mKey(rOtherVariable.mKey){}

      /// Destructor.
      virtual ~VariableData(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /** This operator return the key. by this method user can use Variable
	  as argument for the places which variable key is needed. */
      operator size_t() const {return mKey;}
      
      ///@}
      ///@name Operations
      ///@{
      
      virtual void* Clone(const void* pSource) const {return 0;}
        
      virtual void* Copy(const void* pSource, void* pDestination) const {return 0;}
        
      virtual void Assign(const void* pSource, void* pDestination) const {}
        
      virtual void AssignZero(void* pDestination) const {}
        
      virtual void Destruct(void* pSource) const {}
      
      virtual void Delete(void* pSource) const {}
      
      virtual void Print(const void* pSource, std::ostream& rOStream) const {}
      
      ///@}
      ///@name Access
      ///@{ 
      
      KeyType Key() const {return mKey;}
    
	  /// NOTE: This function is for internal use and not 
	  /// to change arbitrary any variable's key
      void SetKey(KeyType NewKey) {mKey = NewKey;}
    
      const std::string& Name() const {return mName;}

      
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
	  buffer << mName << " variable data";
	  return buffer.str();
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << mName << " variable data";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream <<" #" << static_cast<unsigned int>(mKey);
	}
      
            
      ///@}      
      ///@name Friends
      ///@{
      
      friend bool operator==(const VariableData& rFirstVariable, const VariableData& rSecondVariable)
	{
	  return (rFirstVariable.mKey == rSecondVariable.mKey);
	}
            
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
        
      VariableData& operator=(const VariableData& rOtherVariable)
      {
	mName = rOtherVariable.mName;
	mKey = rOtherVariable.mKey;
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
      
      /// Constructor.
/*       VariableData(const std::string& NewName) : mName(NewName), mKey(gCounter++){} */
      //VariableData(const std::string& NewName) : mName(NewName), mKey(Counter<VariableData>::Increment()){}
      VariableData(const std::string& NewName) : mName(NewName), mKey(0){}
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      std::string mName;
        
      /** Key value of this variable. Each variable will be locate by this
	  value in each data structure. Variable constructor will initialize it. */
      KeyType mKey;

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
      
      /** default constructor is un accessible due to the fact that
	  each variable must have a name defined.*/
      VariableData();
        
      ///@}    
        
    }; // Class VariableData 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    VariableData& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const VariableData& rThis)
    {
      rThis.PrintInfo(rOStream);
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_VARIABLE_DATA_H_INCLUDED  defined 


