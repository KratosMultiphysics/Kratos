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


#if !defined(KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED )
#define  KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/ublas_interface.h"


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
  class GraphColoringProcess : public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of GraphColoringProcess
      KRATOS_CLASS_POINTER_DEFINITION(GraphColoringProcess);
      typedef std::size_t SizeType;
      typedef std::size_t IndexType;
      typedef matrix<int> GraphType;
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Constructor.
      GraphColoringProcess(int NumberOfPartitions, GraphType& rDomainsGraph, GraphType& rDomainsColoredGraph, int& rMaxColor):
	mNumberOfPartitions(NumberOfPartitions), mrMaxColor(rMaxColor), mrDomainsGraph(rDomainsGraph), mrDomainsColoredGraph(rDomainsColoredGraph)
	{}

      /// Destructor.
      virtual ~GraphColoringProcess(){}
      

      ///@}
      ///@name Operators 
      ///@{

      void operator()()
	{
	  Execute();
	}
      
      
      ///@}
      ///@name Operations
      ///@{

      virtual void Execute()
      {
	mrMaxColor = 0;
	// Initializing the coloered graph. -1 means no connection
	// TODO: I have to change this part and create this matrix using max color number
	mrDomainsColoredGraph = ScalarMatrix(mNumberOfPartitions, 2*mNumberOfPartitions, -1.00);
		    
	// Start coloring...
	for(SizeType i = 0 ; i < mrDomainsGraph.size1() ; i++) // for each domain
	  for(SizeType j = i + 1 ; j < mrDomainsGraph.size2() ; j++) // finding neighbor domains
	    if(mrDomainsGraph(i,j) != 0.00) // domain i has interface with domain j
	      for(SizeType color = 0 ; color < mrDomainsColoredGraph.size2() ; color++) // finding color
		if((mrDomainsColoredGraph(i,color) == -1.00) && (mrDomainsColoredGraph(j,color) == -1.00)) // the first unused color 
		  {
		    mrDomainsColoredGraph(i,color) = j;
		    mrDomainsColoredGraph(j,color) = i;
		    if(mrMaxColor < static_cast<int>(color + 1))
		      mrMaxColor = color + 1;
		    break;
		  }
      }
      
      
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
	  return "GraphColoringProcess";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "GraphColoringProcess";
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
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
      int mNumberOfPartitions;
      int& mrMaxColor;
      GraphType& mrDomainsGraph;
      GraphType& mrDomainsColoredGraph;
        
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
      GraphColoringProcess& operator=(GraphColoringProcess const& rOther);

      /// Copy constructor.
      //GraphColoringProcess(GraphColoringProcess const& rOther);

        
      ///@}    
        
    }; // Class GraphColoringProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    GraphColoringProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const GraphColoringProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_GRAPH_COLORING_PROCESS_H_INCLUDED  defined 


