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


#if !defined(KRATOS_MOVE_MESH_PROCESS_INCLUDED )
#define  KRATOS_MOVE_MESH_PROCESS_INCLUDED



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
	*/

	class MoveMeshProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of MoveMeshProcess
		KRATOS_CLASS_POINTER_DEFINITION(MoveMeshProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MoveMeshProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}

		/// Destructor.
		virtual ~MoveMeshProcess()
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
			KRATOS_TRY;

			double dt = mr_model_part.GetProcessInfo()[DELTA_TIME];

			//temporarily save the incremental displacement in the displacement variable
			for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				i != mr_model_part.NodesEnd() ; ++i)
			{
					const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
					i->X() = i->X0() + disp[0];
					i->Y() = i->Y0() + disp[1];
					i->Z() = i->Z0() + disp[2];

					array_1d<double,3>& mesh_vel = i->FastGetSolutionStepValue(MESH_VELOCITY);
					noalias(mesh_vel) = disp;
					noalias(mesh_vel) -= i->FastGetSolutionStepValue(DISPLACEMENT,1);
					mesh_vel /= dt;

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
			return "MoveMeshProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MoveMeshProcess";
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
		int mNumberOfLoops;


		///@} 
		///@name Private Operators
		///@{ 
		void SmoothingLoopOnNodes()
		{
			KRATOS_TRY
			//performs a laplacian smoothing over the incremental displacements 
			for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				i != mr_model_part.NodesEnd() ; ++i)
			{
				WeakPointerVector< Node<3> >& neighbours = (i->GetValue(NEIGHBOUR_NODES));
				//WeakPointerVector< Element >& el_neigh = (i->GetValue(NEIGHBOUR_ELEMENTS));


				double nneigh = neighbours.size();
				if( nneigh != 0.00 && //if it has neighoburs
					i->GetSolutionStepValue(IS_BOUNDARY) == 0  ) //and it is not on the boundary
				{
					array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
					for(WeakPointerVector< Node<3> >::iterator ii = neighbours.begin() ; 
						ii != neighbours.end() ; ++ii)
					{
						noalias(disp) += ii->FastGetSolutionStepValue(DISPLACEMENT);
					}
					disp/=nneigh;
				} 
			}
			KRATOS_CATCH("")
		}

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
		MoveMeshProcess& operator=(MoveMeshProcess const& rOther);

		/// Copy constructor.
		//MoveMeshProcess(MoveMeshProcess const& rOther);


		///@}    

	}; // Class MoveMeshProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MoveMeshProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MoveMeshProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_MOVE_MESH_PROCESS_INCLUDED  defined 


