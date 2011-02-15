//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_SAVE_CONDITIONS_PROCESS_INCLUDED )
#define  KRATOS_SAVE_CONDITIONS_PROCESS_INCLUDED



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
//#include "custom_elements/updated_lagrangian_fluid.h"
//#include "custom_elements/updated_lagrangian_fluid3D.h"
//#include "custom_elements/updated_lagrangian_fluid_inc.h"
//#include "custom_elements/updated_lagrangian_fluid3D_inc.h"


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

	class SaveConditionsProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(SaveConditionsProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		SaveConditionsProcess()
		{
		}

		/// Destructor.
		virtual ~SaveConditionsProcess()
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

		void  SaveStructureConditions(ModelPart& fluid_model_part, ModelPart& face_model_part,  int domain_size)
		{
			KRATOS_TRY
			face_model_part.Conditions().clear();
			//number of structure nodes
			KRATOS_WATCH("SAVING STRUCTURE CONDITION")
		    	for(ModelPart::ConditionsContainerType::iterator ic = fluid_model_part.ConditionsBegin() ; 
					ic != fluid_model_part.ConditionsEnd() ; ++ic)
			{		
				//PointerVector<Element> struct_elements_list;
				//check number of structure nodes
				unsigned int n_struct=0;
				//unsigned int n_fluid = 0;
				for (unsigned int i=0;i<ic->GetGeometry().size();i++)
				{
					n_struct += int( ic->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE) );
					//n_fluid += int( ic->GetGeometry()[i].FastGetSolutionStepValue(IS_FLUID) );
					
				}
				
				if(n_struct==ic->GetGeometry().size())// && n_fluid != im->GetGeometry().size())
				{
					face_model_part.Conditions().push_back(*(ic.base()));	
					for (unsigned int i=0;i<ic->GetGeometry().size();i++)
					{
					ic->GetGeometry()[i].FastGetSolutionStepValue(IS_FLUID)=0.0; 
					}
								
				}
				
				//3D problem involving membrane moving in fluid (only for this)
				
				
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
			return "SaveConditionsProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "SaveConditionsProcess";
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
//		SaveConditionsProcess& operator=(SaveConditionsProcess const& rOther);

		/// Copy constructor.
//		SaveConditionsProcess(SaveConditionsProcess const& rOther);


		///@}    

	}; // Class SaveConditionsProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,   
		SaveConditionsProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const SaveConditionsProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_SAVE_STRUCTURE_CONDITIONS_PROCESS_INCLUDED  defined 


