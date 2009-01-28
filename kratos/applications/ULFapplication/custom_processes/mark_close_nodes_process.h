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
//   Revision:            $Revision: 1.3 $ 
//
//  

#if !defined(KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED )
#define  KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED



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

	class MarkCloseNodesProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(MarkCloseNodesProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MarkCloseNodesProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~MarkCloseNodesProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

		


		///@}
		///@name Operations
		///@{

		void MarkCloseNodes(const double admissible_distance_factor)
		{
		KRATOS_TRY
			
			double fact2 = admissible_distance_factor*admissible_distance_factor;
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
			{
				if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 0) //if it is not a wall node i can erase
				{
					double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
					hnode2 *= hnode2; //take the square

					//loop on neighbours and erase if they are too close
					for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
									i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
					{
						if(i->GetValue(ERASE_FLAG) == false) //we can erase the current node only if the neighb is not to be erased
						{
							double dx = i->X() - in->X();
							double dy = i->Y() - in->Y();
							double dz = i->Z() - in->Z();
							
							double dist2 = dx*dx + dy*dy + dz*dz;

							if(dist2 < fact2 *  hnode2)
								in->GetValue(ERASE_FLAG) = true;
						}
					}
				}
			}
		/*
		this was my old version. now Riccardos version is implemented
		for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
		{
			
				const double& X0 = in->X();		const double& Y0 = in->Y();
				KRATOS_WATCH("ENTERED MARKING CLOSE NODES FUCNTION!");
			
				for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
									i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
				{
				const double& X1 = i->X();	const double& Y1 = i->Y();
				const double& dist_sq = (X1-X0)*(X1-X0)+(Y1-Y0)*(Y1-Y0);
				//if (dist_sq<(i->GetValue(NODAL_H))*(i->GetValue(NODAL_H)) && in->GetId()>i->GetId())
				if (dist_sq<0.005 && in->GetId()>i->GetId())
					{
					if (i->FastGetSolutionStepValue(IS_STRUCTURE)==false)
						{
						i->GetValue(ERASE_FLAG)= true;
						KRATOS_WATCH("ERASING NODE!!!!!!");
						KRATOS_WATCH(in->GetId());
						}

					}
				
				}
				
				
		}
		*/	
			
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
			return "MarkCloseNodesProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MarkCloseNodesProcess";
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
//		MarkCloseNodesProcess& operator=(MarkCloseNodesProcess const& rOther);

		/// Copy constructor.
//		MarkCloseNodesProcess(MarkCloseNodesProcess const& rOther);


		///@}    

	}; // Class MarkCloseNodesProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MarkCloseNodesProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MarkCloseNodesProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_MARK_CLOSE_NODES_PROCESS_INCLUDED  defined 


