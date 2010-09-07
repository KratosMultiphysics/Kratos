/*
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_SUBSCALE_ERRORESTIMATE_PROCESS_H_INCLUDED )
#define  KRATOS_SUBSCALE_ERRORESTIMATE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/model_part.h"

//tools for adaptivity
// #include "custom_utilities/adaptivity_utilities.h" 


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

	class SubscaleEstimatorProcess 
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of SubscaleEstimatorProcess
		KRATOS_CLASS_POINTER_DEFINITION(SubscaleEstimatorProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		SubscaleEstimatorProcess(ModelPart& model_part,int domain_size,  double limit_ratio)
			: mr_model_part(model_part),mlimit_ratio(limit_ratio),mdomain_size(domain_size)
		{
		}

		/// Destructor.
		virtual ~SubscaleEstimatorProcess()
		{
		}


		///@}
		///@name Operators 
		///@{



		///@}
		///@name Operations
		///@{

		//generate a list of new nodes refining accordingly to the REFINE_FLAG
		void Execute()
		{
			KRATOS_TRY

// 			double nu = mr_model_part.GetMesh().GetProperties(1)[VISCOSITY];
// 			double density = mr_model_part.GetMesh().GetProperties(1)[DENSITY];

			array_1d<double,3> aux;
			array_1d<double,3> u_tilda;
			double minimum_energy = 0.0001;
			double h;
			
// 			//reset flags
// 			for(ModelPart::ElementIterator ie = mr_model_part.ElementsBegin() ; 
// 				ie != mr_model_part.ElementsEnd() ; ++ie)
// 			{
// 				ie->GetValue(REFINE_FLAG) = 0;
// 
// 			}
			

			for(ModelPart::NodeIterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
			{
				const double density = in->FastGetSolutionStepValue(DENSITY);
				const double nu = in->FastGetSolutionStepValue(VISCOSITY);
				const array_1d<double,3>& v = in->FastGetSolutionStepValue(VELOCITY);
				const array_1d<double,3>& w = in->FastGetSolutionStepValue(MESH_VELOCITY);
				noalias(aux) = v; noalias(aux)-=w;
				const double nodal_area = in->FastGetSolutionStepValue(NODAL_MASS)/density;
				const double pressure = in->FastGetSolutionStepValue(PRESSURE);

				//calculate tau
				double norm_u = sqrt( aux[0]*aux[0] + aux[1]*aux[1] + aux[2]*aux[2] );
				if(mdomain_size == 2)
					h = sqrt(nodal_area);
				else
					h = pow(nodal_area,0.3333333333333);

				double tau = 1.00 / ( 4.0*nu/(h*h) + 2.0*norm_u/h );

//double aa = pow(4.0*nu/(h*h),2);
//double bb = pow(2.0*norm_u/(h),2);
//tau = sqrt(1.0/(aa+bb));

				//ESTIMATE SUBSCALE VELOCITY (projections already include density)
				noalias(u_tilda) = in->FastGetSolutionStepValue(PRESS_PROJ);
				noalias(u_tilda) += in->FastGetSolutionStepValue(CONV_PROJ);
				u_tilda *= tau;

				//calculate total fluid energy on the node
				double energy_tot = 0.5*density*inner_prod(v,v) + fabs(pressure); 
				double energy_tilda = 0.5*inner_prod(u_tilda,u_tilda);

				//calculate the ratio (error)
				double ratio = 0.0;
				//if the energy is 0 this comparison makes no sense
				if(energy_tot > minimum_energy) 
					ratio = energy_tilda/energy_tot;
					//ratio = energy_tilda;

				//saving the ratio as error
				in->FastGetSolutionStepValue(ERROR_RATIO) = ratio;
				

			}


// 			//decide which elements to refine
// 			for(ModelPart::ElementIterator ie = mr_model_part.ElementsBegin() ; 
// 				ie != mr_model_part.ElementsEnd() ; ++ie)
// 			{
// 				Geometry<Node<3> >& geom = ie->GetGeometry();
// 
// 				//calculate avg error
// 				double avg_err = 0.00;
// 				for(int i=0; i<geom.size(); i++)
// 				{
// 					avg_err += geom[i].FastGetSolutionStepValue(ERROR_RATIO);
// 				}
// 				avg_err /= double(geom.size());
// 
// 				if(avg_err > mlimit_ratio)
// 					ie->GetValue(REFINE_FLAG) = 1;
// 
// 			}


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
			return "SubscaleEstimatorProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "SubscaleEstimatorProcess";
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
		double 	mlimit_ratio;
		int mdomain_size;

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
// 		SubscaleEstimatorProcess& operator=(SubscaleEstimatorProcess const& rOther);

		/// Copy constructor.
// 		SubscaleEstimatorProcess(SubscaleEstimatorProcess const& rOther);


		///@}    

	}; // Class SubscaleEstimatorProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		SubscaleEstimatorProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const SubscaleEstimatorProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_SUBSCALE_ERRORESTIMATE_PROCESS_H_INCLUDED  defined 


