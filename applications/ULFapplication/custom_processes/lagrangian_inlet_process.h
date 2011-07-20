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


#if !defined(KRATOS_LAGRANGIAN_INLET_PROCESS_INCLUDED )
#define  KRATOS_LAGRANGIAN_INLET_PROCESS_INCLUDED



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
	//adds nodes to the inlet boundary ... every "insertion_time_step" time
	/** Detail class definition.
	*/

	class LagrangianInletProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of LagrangianInletProcess
		KRATOS_CLASS_POINTER_DEFINITION(LagrangianInletProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		LagrangianInletProcess(ModelPart& model_part, double insertion_time_step, array_1d<double,3> inlet_vel)
			: mr_model_part(model_part)
		{
			KRATOS_TRY
			minsertion_time = 0.00;
			minsertion_time_step = insertion_time_step;
			//by default we set the inlet vel to 0
			this->minlet_vel=inlet_vel;
			KRATOS_CATCH("")
		}

		/// Destructor.
		virtual ~LagrangianInletProcess()
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

			double time = mr_model_part.GetProcessInfo()[TIME];

			if(minsertion_time == 0.00)
			{
				minsertion_time_step = EstimateInsertionTime();
				minsertion_time = time + minsertion_time_step;
			}
			KRATOS_WATCH(time)
			KRATOS_WATCH(minsertion_time)
			if(time >= minsertion_time)
			{

				//detecting number of new nodes 
				int new_nodes_number = 0;
				for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
				{
					if(in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) == 1)
						new_nodes_number += 1;
				}
				KRATOS_WATCH("SHOULD ADD NODES")
				if(new_nodes_number != 0)
				{
					//allocating the memory needed
					int old_size = mr_model_part.Nodes().size();
					(mr_model_part.Nodes()).reserve(old_size + new_nodes_number);

					for(int i = 0; i<old_size; i++)
					{
						ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() + i; 

						if(in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) == 1)
						{
							//first of all set the node as normal fluid (not structure)
							in->FastGetSolutionStepValue(IS_STRUCTURE) = 0;
							in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
							in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) = 0;
							in->FastGetSolutionStepValue(IS_FREE_SURFACE) = 0;
							in->Free(DISPLACEMENT_X);
							in->Free(DISPLACEMENT_Y);
							in->Free(DISPLACEMENT_Z);
							
							//make a copy of the 
							//Node<3>::Pointer new_node = AdaptivityUtils::AddNode(in->X(),in->Y(),in->Z(),mr_model_part);
							Node<3>::DofsContainerType& reference_dofs = (mr_model_part.NodesBegin())->GetDofs();

							int id = mr_model_part.Nodes().size();
							id++;
							//KRATOS_WATCH(id)
							double init_pos_x=in->X()-in->FastGetSolutionStepValue(DISPLACEMENT_X);
							double init_pos_y=in->Y()-in->FastGetSolutionStepValue(DISPLACEMENT_Y);
							double init_pos_z=in->Z()-in->FastGetSolutionStepValue(DISPLACEMENT_Z);

							//KRATOS_WATCH(init_pos_x)
							//KRATOS_WATCH(init_pos_y)
							//KRATOS_WATCH(init_pos_z)

							Node<3>::Pointer new_node = mr_model_part.CreateNewNode(id, init_pos_x, init_pos_y, init_pos_z);//in->X(), in->Y(), in->Z());

		                                        new_node->SetBufferSize(mr_model_part.NodesBegin()->GetBufferSize() );	
							//KRATOS_WATCH(new_node->GetBufferSize())						

							for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
								{
									Node<3>::DofType& rDof = *iii;
									Node<3>::DofType::Pointer p_new_dof = new_node->pAddDof( rDof );
						
									(p_new_dof)->FreeDof();
									//(p_new_dof)->FixDof();
			//                                                (p_new_dof)->EquationId() = -1;
									KRATOS_WATCH(" ADDED NODES")

								}
							
							
							//new_node->X0() = init_pos_x;
							//new_node->Y0() = init_pos_y;
							//new_node->Z0() = init_pos_z;

							//array_1d<double,3>& dnew = new_node->FastGetSolutionStepValue(DISPLACEMENT);
							//const array_1d<double,3>& d = in->FastGetSolutionStepValue(DISPLACEMENT);

							//here we read the velocity of the Lagrangian inlet
							const array_1d<double,3>& inlet_vel = in->FastGetSolutionStepValue(VELOCITY,1);
							
							//array_1d<double,3> d=this->minlet_vel*dt;
							//KRATOS_WATCH(inlet_vel)
							array_1d<double,3> d=inlet_vel*dt;

							size_t sizee = 0;							
							new_node->FastGetSolutionStepValue(DISPLACEMENT)+=d;
							//KRATOS_WATCH(new_node->FastGetSolutionStepValue(DISPLACEMENT))
							
							//new_node->FastGetSolutionStepValue(VELOCITY) = this->minlet_vel;//in->FastGetSolutionStepValue(VELOCITY);
							//in->FastGetSolutionStepValue(VELOCITY) = this->minlet_vel;//in->FastGetSolutionStepValue(VELOCITY);
							new_node->FastGetSolutionStepValue(VELOCITY) = inlet_vel;//in->FastGetSolutionStepValue(VELOCITY);
							in->FastGetSolutionStepValue(VELOCITY) = inlet_vel;//in->FastGetSolutionStepValue(VELOCITY);
							
							new_node->FastGetSolutionStepValue(NODAL_AREA) = in->FastGetSolutionStepValue(NODAL_AREA);
							new_node->FastGetSolutionStepValue(PRESSURE,1) = in->FastGetSolutionStepValue(PRESSURE,1);
							new_node->FastGetSolutionStepValue(PRESSURE) = in->FastGetSolutionStepValue(PRESSURE);
							new_node->FastGetSolutionStepValue(BULK_MODULUS) = in->FastGetSolutionStepValue(BULK_MODULUS);
							new_node->FastGetSolutionStepValue(VISCOSITY) = in->FastGetSolutionStepValue(VISCOSITY);
							new_node->FastGetSolutionStepValue(DENSITY) = in->FastGetSolutionStepValue(DENSITY);
							new_node->FastGetSolutionStepValue(BODY_FORCE) = in->FastGetSolutionStepValue(BODY_FORCE);
							new_node->FastGetSolutionStepValue(IS_FLUID)=1.0;
							new_node->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
							//new_node->FastGetSolutionStepValue(ACCELERATION) = in->FastGetSolutionStepValue(PRESS_PROJ);
							new_node->Fix(DISPLACEMENT_X);
							new_node->Fix(DISPLACEMENT_Y);
							new_node->Fix(DISPLACEMENT_Z);

							//new_node->FastGetSolutionStepValue(IS_STRUCTURE) = 1;
							new_node->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) = 1;
							double nodal_h = in->FastGetSolutionStepValue(NODAL_H);
							new_node->FastGetSolutionStepValue(NODAL_H) = nodal_h;
							
						}
					}
					minsertion_time_step = EstimateInsertionTime();
					KRATOS_WATCH(minsertion_time_step);
				
					//update the time for the next insertion
					minsertion_time = time + minsertion_time_step;

				}		
			}	
			//if we do not insert in this step - we simply have to move
			else
			{			
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
				{
					if(in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) == 1)
						{
						array_1d<double,3> old_disp=in->FastGetSolutionStepValue(DISPLACEMENT,1);
						const array_1d<double,3>& inlet_vel = in->FastGetSolutionStepValue(VELOCITY,1);
						//KRATOS_WATCH(inlet_vel)
						//array_1d<double,3> inc_disp=this->minlet_vel*dt;
						array_1d<double,3> inc_disp=inlet_vel*dt;	
						in->FastGetSolutionStepValue(DISPLACEMENT)=inc_disp+old_disp;
						//in->FastGetSolutionStepValue(VELOCITY)=this->minlet_vel;
						in->FastGetSolutionStepValue(VELOCITY)=inlet_vel;
						in->Fix(DISPLACEMENT_X);
						in->Fix(DISPLACEMENT_Y);
						in->Fix(DISPLACEMENT_Z);
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
			return "LagrangianInletProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "LagrangianInletProcess";
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
		double minsertion_time_step;
		double minsertion_time;
		array_1d<double,3> minlet_vel;


		///@} 
		///@name Private Operators
		///@{ 
		double EstimateInsertionTime()
		{
			KRATOS_TRY

				double dt_estimate = 100.0;

			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin(); in!=mr_model_part.NodesEnd(); in++)
			{
				if(in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) == 1)
				{
					const array_1d<double,3> v = this->minlet_vel;//in->FastGetSolutionStepValue(VELOCITY);
					//KRATOS_WATCH("INLET VELOCITY")
					//KRATOS_WATCH(v)
					double nodal_h = in->FastGetSolutionStepValue(NODAL_H);

					//estimating the next time step
					double normv = norm_2(v);
					double dtcandidate = nodal_h / normv;
					KRATOS_WATCH(dtcandidate)
					if( dtcandidate < dt_estimate)
						dt_estimate = dtcandidate;
					dt_estimate*=1.1;
					KRATOS_WATCH(dt_estimate)
				}
			}

			return dt_estimate;

			KRATOS_CATCH("");
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
		LagrangianInletProcess& operator=(LagrangianInletProcess const& rOther);


		///@}    

	}; // Class LagrangianInletProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		LagrangianInletProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const LagrangianInletProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_INLET_PROCESS_INCLUDED  defined 


