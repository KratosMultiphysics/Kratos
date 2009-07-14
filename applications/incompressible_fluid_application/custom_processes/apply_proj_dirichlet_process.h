//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED )
#define  KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED


//This process loops over the destination model part elements, and fixes the velocity (applies Dirichlet) to the boundary nodes(on the fictitious side) of the elements, that belong to the interface



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h" 
#include "incompressible_fluid_application.h"

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

	class ApplyProjDirichletProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(ApplyProjDirichletProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		ApplyProjDirichletProcess()			
		{
		}

		/// Destructor.
		virtual ~ApplyProjDirichletProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

//		void operator()()
//		{
//			SaveStructure();
//		}


		///@}
		///@name Operations
		///@{

		void ApplyProjDirichlet(ModelPart& full_model_part, ModelPart& aux_conditions_model_part)
		{
		KRATOS_TRY
		int n_int;

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE,1);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=false;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=false;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=false;
				}

		}
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=true;
				}

		}

		
		
				
		
		
		//first we remove the Dirichlet conditions from the nodes that were defining the interface in the previous step:

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE,1);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE,1);
			
			if (n_int>0.0)
				{
								
				im->GetGeometry()[0].Free(VELOCITY_X);
				im->GetGeometry()[0].Free(VELOCITY_Y);

				im->GetGeometry()[1].Free(VELOCITY_X);
				im->GetGeometry()[1].Free(VELOCITY_Y);

				im->GetGeometry()[2].Free(VELOCITY_X);
				im->GetGeometry()[2].Free(VELOCITY_Y);

				im->GetGeometry()[0].Free(AUX_VEL_X);
				im->GetGeometry()[0].Free(AUX_VEL_Y);

				im->GetGeometry()[1].Free(AUX_VEL_X);
				im->GetGeometry()[1].Free(AUX_VEL_Y);

				im->GetGeometry()[2].Free(AUX_VEL_X);
				im->GetGeometry()[2].Free(AUX_VEL_Y);

				im->GetGeometry()[0].Free(PRESSURE);				
				im->GetGeometry()[1].Free(PRESSURE);
				im->GetGeometry()[2].Free(PRESSURE);
				}		

		}
		

		if (aux_conditions_model_part.Conditions().size()>0);
		{
			//and corerct the velocity of the elements that were "intersected" and stored and solved in aux_condition_model_part
			for(ModelPart::ConditionsContainerType::iterator ic = aux_conditions_model_part.ConditionsBegin() ; 
					ic != aux_conditions_model_part.ConditionsEnd() ; ++ic)
			{	  
			
				n_int=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
				n_int+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
				n_int+=ic->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
				//elements with n_int=3 lie inside the fictitious domain, and therefore are of no interest
				if (n_int==1 || n_int==2)
			
					{
					//and now fix the velocity of the IS_INTERFACE nodes (those are lying inside of the fictitious domain)
					for (int i=0;i<2;i++)
						{
						if (ic->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE)==1.0)
							{
							ic->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)=ic->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);
							//KRATOS_WATCH(im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL))
						
							//im->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE)=0.0;

							ic->GetGeometry()[i].Fix(VELOCITY_X);
							ic->GetGeometry()[i].Fix(VELOCITY_Y);
							//im->GetGeometry()[i].Fix(PRESSURE);
						
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
			return "ApplyProjDirichletProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "ApplyProjDirichletProcess";
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
		//ModelPart& mr_fluid_model_part;
		//ModelPart& mr_structure_model_part;
		
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
//		ApplyProjDirichletProcess& operator=(ApplyProjDirichletProcess const& rOther);

		/// Copy constructor.
//		ApplyProjDirichletProcess(ApplyProjDirichletProcess const& rOther);


		///@}    

	}; // Class ApplyProjDirichletProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,   
		ApplyProjDirichletProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const ApplyProjDirichletProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_PROJ_DIRICHLET_PROCESS_INCLUDED  defined 


