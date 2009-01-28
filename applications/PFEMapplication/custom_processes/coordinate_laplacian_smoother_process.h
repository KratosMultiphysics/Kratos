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
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_COORDINATE_LAPLACIAN_SMOOTHER_PROCESS_INCLUDED )
#define  KRATOS_COORDINATE_LAPLACIAN_SMOOTHER_PROCESS_INCLUDED



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
	This process calculates and moves the mesh in a Lagrangian way
	using the formula
	d(n+1) = d(n) + v(n)*dt + 0.5*a(n)*dt*dt
	to do so it requires to recovery the material acceleration a(n) from the knowledge
	of the eulerian solution at the steps n and n+1
	a(n) = (v(n+1) - v(n))/Dt + dv/dx*(v-vm)
	*/

	class CoordinateLaplacianSmootherProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of CoordinateLaplacianSmootherProcess
		KRATOS_CLASS_POINTER_DEFINITION(CoordinateLaplacianSmootherProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		CoordinateLaplacianSmootherProcess(ModelPart& model_part, int number_of_loops = 20, double reduction_factor = 1.0)
			: mr_model_part(model_part),mreduction_factor(reduction_factor)
		{
			KRATOS_TRY
			mNumberOfLoops = number_of_loops;
			KRATOS_CATCH("")
		}

		/// Destructor.
		virtual ~CoordinateLaplacianSmootherProcess()
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

			double relaxation_param = 1.00/double(mNumberOfLoops);

			relaxation_param *= mreduction_factor;

			double Dt = mr_model_part.GetProcessInfo()[DELTA_TIME];

			for(unsigned int iloops=0; iloops<mNumberOfLoops; iloops++)
				SmoothingLoopOnNodes(relaxation_param,Dt);


			//calculate the mesh velocity and ripristinate the disp variable to the total value
			for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				i != mr_model_part.NodesEnd() ; ++i)
			{
				if(i->FastGetSolutionStepValue(IS_BOUNDARY) == 0)
				{
					array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);

					disp[0] = i->X() - i->X0();
					disp[1] = i->Y() - i->Y0();
					disp[2] = i->Z() - i->Z0();
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
			return "CoordinateLaplacianSmootherProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "CoordinateLaplacianSmootherProcess";
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
		unsigned int mNumberOfLoops;
		double mreduction_factor;


		///@} 
		///@name Private Operators
		///@{ 
		void SmoothingLoopOnNodes(double relaxation_param, double Dt)
		{
			KRATOS_TRY
		
			array_1d<double,3> new_disp;
			array_1d<double,3> zero = ZeroVector(3);
			array_1d<double,3> aux;
			double complement_param = 1.00 - relaxation_param;
			//performs a laplacian smoothing over the incremental displacements 
			for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				i != mr_model_part.NodesEnd() ; ++i)
			{
					WeakPointerVector< Node<3> >& neighbours = (i->GetValue(NEIGHBOUR_NODES));

					double nneigh = neighbours.size();
					if( nneigh != 0.00 && //if it has neighoburs
						i->GetSolutionStepValue(IS_BOUNDARY) == 0  &&
						i->GetSolutionStepValue(IS_STRUCTURE) == 0 ) //and it is not on the boundary
					{
						double xnew = 0.0;
						double ynew = 0.0;
						double znew = 0.0;
						double xold = i->X();
						double yold = i->Y();
						double zold = i->Z();

						for(WeakPointerVector< Node<3> >::iterator ii = neighbours.begin() ; 
							ii != neighbours.end() ; ++ii)
						{
							xnew += ii->X();
							ynew += ii->Y();
							znew += ii->Z();
						}
						xnew/=nneigh; ynew/=nneigh; znew/=nneigh;

						i->X() = complement_param * xold + relaxation_param*xnew ;
						i->Y() = complement_param * yold + relaxation_param*ynew ;
						i->Z() = complement_param * zold + relaxation_param*znew ;

// 						//verify that the displacement does not exceed the maximum allowable
// 						const array_1d<double,3>& vel = i->GetSolutionStepValue(VELOCITY);
// 						double v = inner_prod(vel,vel); v = sqrt(v);
// 
// 						double max_disp = v * Dt  ;
// // std::cout << v << " " << Dt << " " << relaxation_param << " " << max_disp << std::endl;
// 						aux[0] = i->X() - xold;
// 						aux[1] = i->Y() - yold;
// 						aux[2] = i->Z() - zold;
// 						double dist = inner_prod(aux,aux); dist = sqrt(v);
// // std::cout << dist << std::endl;
// 
// 						if(dist > max_disp && dist > 0.00)
// 						{
// 							double ratio = max_disp/dist;
// 							i->X() = xold + ratio * aux[0];
// 							i->Y() = yold + ratio * aux[1];
// 							i->Z() = zold + ratio * aux[2];
// // std::cout << ratio << " " << aux << std::endl;
// // std::cout << ratio << " " << ratio * aux << std::endl;
// 						}
						
						
						
					
						
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
		CoordinateLaplacianSmootherProcess& operator=(CoordinateLaplacianSmootherProcess const& rOther);


		///@}    

	}; // Class CoordinateLaplacianSmootherProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		CoordinateLaplacianSmootherProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const CoordinateLaplacianSmootherProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_COORDINATE_LAPLACIAN_SMOOTHER_PROCESS_INCLUDED  defined 


