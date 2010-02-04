/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.5 $ 
//
//  this process puts a IS_FLUID flag on the nodes of fluid elements 

#if !defined(KRATOS_MARK_FLUID_PROCESS_INCLUDED )
#define  KRATOS_MARK_FLUID_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
//#include "custom_utilities/geometry_utilities2D.h"
#include "custom_elements/updated_lagrangian_fluid.h"
#include "custom_elements/updated_lagrangian_fluid3D.h"
#include "custom_elements/updated_lagrangian_fluid_inc.h"
#include "custom_elements/updated_lagrangian_fluid3D_inc.h"


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
		Update the PRESSURE_FORCE on the nodes

		
	*/

	class MarkFluidProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(MarkFluidProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MarkFluidProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~MarkFluidProcess()
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

		ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
		double temp;

		for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
		{		
		in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
		}
		//set to 1 all the nodes surrounded by at least one fluid element
		for(ModelPart::ElementsContainerType::iterator ie = mr_model_part.ElementsBegin(); ie!=mr_model_part.ElementsEnd(); ie++)
		{		
		  ie->Calculate(IS_FLUID,temp,proc_info);	
		}	
			
/*		for(ModelPart::NodesContainerType::const_iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
			{		
				if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
					{
					in->FastGetSolutionStepValue(IS_FLUID) = 1.0;
					
					}
				else 
							in->FastGetSolutionStepValue(IS_FLUID) = 0.0;
			}*/
		KRATOS_CATCH("")
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
			return "MarkFluidProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MarkFluidProcess";
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
//		MarkFluidProcess& operator=(MarkFluidProcess const& rOther);

		/// Copy constructor.
//		MarkFluidProcess(MarkFluidProcess const& rOther);


		///@}    

	}; // Class MarkFluidProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MarkFluidProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MarkFluidProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_MARK_FLUID_PROCESS_INCLUDED  defined 


