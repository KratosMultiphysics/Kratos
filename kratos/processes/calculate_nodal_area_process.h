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
//   Date:                $Date: 2007-10-31 17:51:34 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
	typedef  ModelPart::NodesContainerType NodesContainerType;
	typedef  ModelPart::ElementsContainerType ElementsContainerType;
	
  
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
  class CalculateNodalAreaProcess 
	: public Process
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CalculateNodalAreaProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalAreaProcess);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      /// avg_elems ------ expected number of neighbour elements per node., 
      /// avg_nodes ------ expected number of neighbour Nodes 
      /// the better the guess for the quantities above the less memory occupied and the fastest the algorithm
      CalculateNodalAreaProcess(ModelPart& model_part, unsigned int domain_size)
		: mr_model_part(model_part), mdomain_size(domain_size)
	{
	}

      /// Destructor.
      virtual ~CalculateNodalAreaProcess()
	{
	}
      

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
		KRATOS_TRY
				
		//set to zero the nodal area
		for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); 
		in!=mr_model_part.NodesEnd(); in++)
		{
			in->FastGetSolutionStepValue(NODAL_AREA) = 0.00;
		}

		if(mdomain_size == 2)
		{
			double area = 0.0;
			for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin(); 
						 i!=mr_model_part.ElementsEnd(); i++)
			{	
					//calculating shape functions values
				Geometry< Node<3> >& geom = i->GetGeometry();
					
				area = GeometryUtils::CalculateVolume2D(geom);
				area *= 0.333333333333333333333333333;


				geom[0].FastGetSolutionStepValue(NODAL_AREA) += area;
				geom[1].FastGetSolutionStepValue(NODAL_AREA) += area;
				geom[2].FastGetSolutionStepValue(NODAL_AREA) += area;
			}
		}
		else if(mdomain_size == 3)
		{
			for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin(); 
						 i!=mr_model_part.ElementsEnd(); i++)
			{	
				double vol;
					//calculating shape functions values
				Geometry< Node<3> >& geom = i->GetGeometry();
					
				vol = GeometryUtils::CalculateVolume3D(geom);
				vol *= 0.25;
					
				geom[0].FastGetSolutionStepValue(NODAL_AREA) += vol;
				geom[1].FastGetSolutionStepValue(NODAL_AREA) += vol;
				geom[2].FastGetSolutionStepValue(NODAL_AREA) += vol;
				geom[3].FastGetSolutionStepValue(NODAL_AREA) += vol;
			}
		}

		mr_model_part.GetCommunicator().AssembleCurrentData(NODAL_AREA);

	      

		KRATOS_CATCH("");

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
	  return "CalculateNodalAreaProcess";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	  rOStream << "CalculateNodalAreaProcess";
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
	ModelPart& mr_model_part;
	unsigned int mdomain_size;
        
        
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
      CalculateNodalAreaProcess& operator=(CalculateNodalAreaProcess const& rOther);

      /// Copy constructor.
      //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);

        
      ///@}    
        
    }; // Class CalculateNodalAreaProcess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    CalculateNodalAreaProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const CalculateNodalAreaProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_CALCULATE_NODAL_AREA_PROCESS_H_INCLUDED  defined 


