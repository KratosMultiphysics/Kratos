//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2007-11-06 12:34:26 $
//   Revision:            $Revision: 1.4 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_MERGE_MODEL_PARTS_PROCESS_INCLUDED )
#define  KRATOS_MERGE_MODEL_PARTS_PROCESS_INCLUDED



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

	class MergeModelPartsProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(MergeModelPartsProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		MergeModelPartsProcess()
			//ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
			//: mr_fluid_model_part(fluid_model_part), mr_structure_model_part(structure_model_part), mr_combined_model_part(combined_model_part)
		{
		}

		/// Destructor.
		virtual ~MergeModelPartsProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

	//	void operator()()
	//	{
	//		MergeParts();
	//	}


		///@}
		///@name Operations
		///@{

		void MergeParts(ModelPart& fluid_model_part, ModelPart& structure_model_part, ModelPart& combined_model_part)
		{
			KRATOS_TRY
			combined_model_part.Elements().clear();
			combined_model_part.Nodes();
			combined_model_part.Nodes().clear();

			combined_model_part.Nodes()=fluid_model_part.Nodes();
			combined_model_part.Elements()=fluid_model_part.Elements();
			combined_model_part.Conditions()=fluid_model_part.Conditions();

			//sorting and renumbering the fluid elements
			unsigned int id=1;
			for(ModelPart::ElementsContainerType::iterator im = combined_model_part.ElementsBegin() ; 
				im != combined_model_part.ElementsEnd() ; ++im)
			{
                                im->SetId(id);
//				im->Id() = id;
				id++;
			}
	/*		KRATOS_WATCH("Last FLUID ID is");
			KRATOS_WATCH(id);

			KRATOS_WATCH("The size of fluid elements container is");
			KRATOS_WATCH(fluid_model_part.Elements().size());
	*/		
			//adding the structure elements to the combined part
			//combined_model_part.Elements().push_back(*( structure_model_part.ElementsBegin().base()));	
			id += 1000;
			for(ModelPart::ElementsContainerType::iterator im = structure_model_part.ElementsBegin() ; 
				im != structure_model_part.ElementsEnd() ; ++im)
			{		
                                im->SetId(id);
//				im->Id() = id;
				combined_model_part.Elements().push_back(*(im.base()));
				id++;
			}
			combined_model_part.Elements().Sort();

			
			//renumbering
			id=1;
			for(ModelPart::ConditionsContainerType::iterator im = combined_model_part.ConditionsBegin() ; 
				im != combined_model_part.ConditionsEnd() ; ++im)
			{
                            im->SetId(id);
//				im->Id() = id;
				id++;
			}
			//adding the structure conditions to the combined part
			id += 1000;
			for(ModelPart::ConditionsContainerType::iterator im = structure_model_part.ConditionsBegin() ; 
				im != structure_model_part.ConditionsEnd() ; ++im)
			{
                                im->SetId(id);
//				im->Id() = id;
				combined_model_part.Conditions().push_back(*(im.base()));
				id++;
			}
			combined_model_part.Conditions().Sort();

			//renumbering
			/*
			unsigned int new_id=1;
			for(ModelPart::ConditionsContainerType::iterator im = structure_model_part.ConditionsBegin() ; 
				im != structure_model_part.ConditionsEnd() ; ++im)
			{		
				im->Id() =new_id;
				new_id++;
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
			return "MergeModelPartsProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "MergeModelPartsProcess";
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
		//ModelPart& mr_combined_model_part;
		
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
//		MergeModelPartsProcess& operator=(MergeModelPartsProcess const& rOther);

		/// Copy constructor.
//		MergeModelPartsProcess(MergeModelPartsProcess const& rOther);


		///@}    

	}; // Class MergeModelPartsProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		MergeModelPartsProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const MergeModelPartsProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_MERGE_MODEL_PARTS_PROCESS_INCLUDED  defined 


