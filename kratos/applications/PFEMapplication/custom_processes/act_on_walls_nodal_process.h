/*
==============================================================================
KratosPFEMApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
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
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_ACT_ON_WALLS_NODAL_PROCESS_H_INCLUDED )
#define  KRATOS_ACT_ON_WALLS_NODAL_PROCESS_H_INCLUDED



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
	//assign a size to each node depending on the distance fof the neighbouring nodes
	*/

	class ActOnWallsNodalProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of ActOnWallsNodalProcess
		KRATOS_CLASS_POINTER_DEFINITION(ActOnWallsNodalProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		ActOnWallsNodalProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~ActOnWallsNodalProcess()
		{
		}

		/// Copy constructor.
		//ActOnWallsNodalProcess(ActOnWallsNodalProcess const& rOther);


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

		//removes the velocity towards the wall for all the nodes which are getting to close
		virtual void Execute()
		{
			KRATOS_TRY
			
			array_1d<double,3> n, final_pos;
			array_1d<double,3> dist, corrective_disp, delta_disp_wall;

			for(ModelPart::NodesContainerType::iterator i = mr_model_part.NodesBegin(); 
				i!=mr_model_part.NodesEnd(); i++)
			{
			//	const WeakPointerVector<Node<3> >& rN = i->GetValue(NEIGHBOUR_NODES);
				
				//if the conditions is composed only by wall nodes remove a percentage of the velocity
				if(   i->FastGetSolutionStepValue(IS_STRUCTURE) == 1 
					&& (i->GetValue(NEIGHBOUR_NODES)).size() != 0)
				{
					
					noalias(n) = i->FastGetSolutionStepValue(NORMAL_TO_WALL);
					double A = inner_prod(n,n);
					A = sqrt(A);
					const array_1d<double,3>& disp_wall = i->FastGetSolutionStepValue(DISPLACEMENT);
					const array_1d<double,3>& disp_wall_old = i->FastGetSolutionStepValue(DISPLACEMENT,1);
					noalias(delta_disp_wall) = disp_wall;
					noalias(delta_disp_wall) -= disp_wall_old;
					if(A != 0)
					{
						n /= A;

						double x = i->X();
						double y = i->Y();
						double z = i->Z();
						double h = i->FastGetSolutionStepValue(NODAL_H);		

						WeakPointerVector<Node<3> >& rN = i->GetValue(NEIGHBOUR_NODES);
						for(WeakPointerVector<Node<3> >::iterator in = rN.begin(); in != rN.end(); in++)
						{
							if(in->FastGetSolutionStepValue(IS_STRUCTURE) == 0) //for nodes that are not of structure
							{
								//calculating the distance along the normal and orhtogonal to the normal
								array_1d<double,3>& disp = in->FastGetSolutionStepValue(DISPLACEMENT);
								array_1d<double,3>& olddisp = in->FastGetSolutionStepValue(DISPLACEMENT,1);
	
								array_1d<double,3> delta_disp = disp;
								noalias(delta_disp) -= olddisp;
	
								//calculating the end of the step position
								final_pos[0] = in->X0(); final_pos[1] = in->Y0(); final_pos[2] = in->Z0();
								noalias(final_pos) += disp;
								
								//calculating the distance at the end of the step
								noalias(dist) = final_pos;
								dist[0] = final_pos[0] - x;
								dist[1] = final_pos[1] - y;
								dist[2] = final_pos[2] - z;
	
								double dist_2 = inner_prod(dist,dist);
								double dist_n = inner_prod(dist,n);
								double dist_n2 = dist_n * dist_n;
								double dist_ort2 = dist_2 - dist_n2;
	
								//verify that the final position of the point is close to the normal in the wall node
								if(dist_ort2 < 0.5*h*h) 
								{
									//verify if the final position of the normal closer than 0.75h to the face (distance along the normal)
									if(dist_n > -0.5 * h)
									{
										double a = inner_prod(delta_disp,n);
										
										//if it is moving towards the outside put it back
										// to the allowed distance of 0.75h on the interior
										if(a > 0)
										{
											double coeff = - 0.5*h - dist_n;
											noalias(disp) += coeff * n;
	noalias(disp) += delta_disp_wall;	
										}
											
									}
										
								}
								
							}
						}
					}
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
			return "ActOnWallsNodalProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "ActOnWallsNodalProcess";
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
		//7/@{ 


		///@}    
		///@name Un accessible methods 
		///@{ 

		/// Assignment operator.
		ActOnWallsNodalProcess& operator=(ActOnWallsNodalProcess const& rOther);



		///@}    

	}; // Class ActOnWallsNodalProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		ActOnWallsNodalProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const ActOnWallsNodalProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_ACT_ON_WALLS_NODAL_PROCESS_H_INCLUDED  defined 


