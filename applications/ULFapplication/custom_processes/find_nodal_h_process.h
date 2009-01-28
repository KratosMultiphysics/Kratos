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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_FIND_NODAL_H_PROCESS_INCLUDED )
#define  KRATOS_FIND_NODAL_H_PROCESS_INCLUDED



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
		calculate the nodal H for all the nodes depending on the min distance
		of the neighbouring nodes.

		lonely nodes are given the average value of the H
	*/

	class FindNodalHProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of FindNodalHProcess
		KRATOS_CLASS_POINTER_DEFINITION(FindNodalHProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		FindNodalHProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~FindNodalHProcess()
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
			ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
			//

			double havg = 0.00;
			double h_nodes = 0;

			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				if((in->GetValue(NEIGHBOUR_NODES)).size() != 0)
				{
					double xc = in->X();
					double yc = in->Y();
					double zc = in->Z();

					double h = 1000.0;
					for( WeakPointerVector< Node<3> >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
									i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
					{
						double x = i->X();
						double y = i->Y();
						double z = i->Z();
						double l = (x-xc)*(x-xc);
						l += (y-yc)*(y-yc);
						l += (z-zc)*(z-zc);

						//if(l>h) h = l;
						if(l<h) h = l;
					}
					h = sqrt(h);
					havg += h;
					h_nodes += 1;

					in->FastGetSolutionStepValue(NODAL_H) = h;
				} 
			}

			//calculate average h
			if(h_nodes == 0)
				KRATOS_ERROR(std::logic_error,"no node has neighbours!!!!","");
			havg /= h_nodes;

			//set the h to the average to the nodes without neighours
			for(ModelPart::NodesContainerType::iterator in = rNodes.begin(); in!=rNodes.end(); in++)
			{
				if((in->GetValue(NEIGHBOUR_NODES)).size() == 0)
				{
					in->FastGetSolutionStepValue(NODAL_H) = havg;
				} 
			}
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
			return "FindNodalHProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "FindNodalHProcess";
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
		double m_min_h;


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
		FindNodalHProcess& operator=(FindNodalHProcess const& rOther);

		/// Copy constructor.
		//FindNodalHProcess(FindNodalHProcess const& rOther);


		///@}    

	}; // Class FindNodalHProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		FindNodalHProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const FindNodalHProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_FIND_NODAL_H_PROCESS_INCLUDED  defined 


