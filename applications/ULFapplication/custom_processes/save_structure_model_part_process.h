//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-02 10:47:21 $
//   Revision:            $Revision: 1.8 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_SAVE_STRUCTURE_MODEL_PART_PROCESS_INCLUDED )
#define  KRATOS_SAVE_STRUCTURE_MODEL_PART_PROCESS_INCLUDED



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

	class SaveStructureModelPartProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(SaveStructureModelPartProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		SaveStructureModelPartProcess()
		{
		}

		/// Destructor.
		virtual ~SaveStructureModelPartProcess()
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

		void  SaveStructure(ModelPart& fluid_model_part, ModelPart& structure_model_part,  int domain_size)
		{
			KRATOS_TRY
				
			//number of structure nodes
			KRATOS_WATCH("SAVING STRUCTURE")
		    for(ModelPart::ElementsContainerType::iterator im = fluid_model_part.ElementsBegin() ; 
					im != fluid_model_part.ElementsEnd() ; ++im)
			{		
				//PointerVector<Element> struct_elements_list;
				//check number of structure nodes
				unsigned int n_struct=0;
				unsigned int n_fluid = 0;
				for (unsigned int i=0;i<im->GetGeometry().size();i++)
				{
					n_struct += int( im->GetGeometry()[i].FastGetSolutionStepValue(IS_STRUCTURE) );
					n_fluid += int( im->GetGeometry()[i].FastGetSolutionStepValue(IS_FLUID) );
					
				}
				
				if(n_struct==im->GetGeometry().size())// && n_fluid != im->GetGeometry().size())
				{
					structure_model_part.Elements().push_back(*(im.base()));					
				}
				
				//3D problem involving membrane moving in fluid (only for this)
				/*
				if (int(domain_size)==3)
				{
				//just to check FOR FLUID-MEMBRANE!!!! otherwise th eone above is correct
				
					if(n_struct==im->GetGeometry().size() && n_struct!= 4)
					{
					//KRATOS_WATCH("Domain size iz 3333333333333333333333333333333333333333333333333333333!!!!!!!")
					structure_model_part.Elements().push_back(*(im.base()));
				
					}
				}
				*/
				
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
			return "SaveStructureModelPartProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "SaveStructureModelPartProcess";
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
//		SaveStructureModelPartProcess& operator=(SaveStructureModelPartProcess const& rOther);

		/// Copy constructor.
//		SaveStructureModelPartProcess(SaveStructureModelPartProcess const& rOther);


		///@}    

	}; // Class SaveStructureModelPartProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		SaveStructureModelPartProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const SaveStructureModelPartProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_SAVE_STRUCTURE_MODEL_PART_PROCESS_INCLUDED  defined 


