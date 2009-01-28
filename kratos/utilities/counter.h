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
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_COUNTER_H_INCLUDED )
#define  KRATOS_COUNTER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


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
  template<class TCountedType>
    class Counter
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Counter
      KRATOS_CLASS_POINTER_DEFINITION(Counter);

      typedef unsigned int SizeType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      // return counter 
      static SizeType GetCounter() 
	{ 
	  return Counter<TCountedType>::msCounter; 
	} 
      
      // return counter 
      static SizeType Increment() 
	{ 
	  return ++Counter<TCountedType>::msCounter; 
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
	  return std::string("Counter");
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << Info();
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	  rOStream << Counter();
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
      
      /// Default constructor.
      Counter() 
	{ 
	} 
            
      /// Copy constructor.
      Counter(Counter<TCountedType> const& rOther) 
	{ 
	} 
 
      /// Destructor.
      virtual ~Counter() 
	{
	} 
 
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
      /// Number of created objects.
      static SizeType msCounter;
        
      ///@} 
      ///@name Member Variables 
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
      
      /// Assignment operator.
      Counter& operator=(Counter const& rOther);

      ///@}    
        
    }; // Class Counter 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<class TCountedType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Counter<TCountedType>& rThis);

  /// output stream function
  template<class TCountedType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Counter<TCountedType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << " : ";
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  // initialize counter with zero 
  template<typename TCountedType>  
    typename Counter<TCountedType>::SizeType Counter<TCountedType>::msCounter = 0; 
  
}  // namespace Kratos.

#endif // KRATOS_COUNTER_H_INCLUDED  defined 


