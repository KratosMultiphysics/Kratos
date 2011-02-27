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


#if !defined(KRATOS_PROCESS_H_INCLUDED )
#define  KRATOS_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"


namespace Kratos
{

  ///@name Kratos Classes
  ///@{
  
  /// The base class for all processes in Kratos.
  /** The process is the base class for all processes and defines a simple interface for them.
      Execute method is used to execute the Process algorithms. While the parameters of this method 
	  can be very different from one Process to other there is no way to create enough overridden 
	  versions of it. For this reason this method takes no argument and all Process parameters must
	  be passed at construction time. The reason is that each constructor can take different set of 
	  argument without any dependency to other processes or the base Process class.
  */
  class Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Process
      KRATOS_CLASS_POINTER_DEFINITION(Process);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Process(){}

      /// Destructor.
      virtual ~Process(){}
      

      ///@}
      ///@name Operators 
      ///@{

	  /// This operator is provided to call the process as a function and simply calls the Execute method.
      void operator()()
	{
	  Execute();
	}
      
      
      ///@}
      ///@name Operations
      ///@{


	  /// Execute method is used to execute the Process algorithms.
      virtual void Execute(){}
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  return "Process";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "Process";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	{
	}
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      Process& operator=(Process const& rOther);

      /// Copy constructor.
      //Process(Process const& rOther);

        
      ///@}    
        
    }; // Class Process 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Process& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Process& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PROCESS_H_INCLUDED  defined 


