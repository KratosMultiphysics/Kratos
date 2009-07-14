//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_PSEUDO_LAG_PROCESS_INCLUDED )
#define  KRATOS_PSEUDO_LAG_PROCESS_INCLUDED


//This process permits one to "REMOVE" the elements that have all nodes "DISABLED" (flag) from the model part
// this reduced model part will be saved

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

	class PseudoLagPartProcess
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(PseudoLagPartProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		PseudoLagPartProcess()
		{
		}

		/// Destructor.
		virtual ~PseudoLagPartProcess()
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

		void SavePseudoLagPart(ModelPart& full_model_part, ModelPart& lagrangian_part, ModelPart& pseudo_lag_part)
		{
		KRATOS_TRY
		int n_int; 
		int n_disabled;
		//clear reduced_model_part
		pseudo_lag_part.Conditions().clear();
		pseudo_lag_part.Elements().clear();
		pseudo_lag_part.Nodes().clear();

		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				pseudo_lag_part.AddElement(*(im.base()));
				}

		}
		
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			
			double n_int=in->FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==1.0)
				{
				pseudo_lag_part.AddNode(*(in.base()));				
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
			return "PseudoLagPartProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "PseudoLagPartProcess";
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
//		PseudoLagPartProcess& operator=(PseudoLagPartProcessconst& rOther);

		/// Copy constructor.
//		PseudoLagPartProcess(PseudoLagPartProcessconst& rOther);


		///@}    

	}; // Class PseudoLagPartProcess

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,   
		PseudoLagPartProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const PseudoLagPartProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_PSEUDO_LAG_PROCESS_INCLUDED  defined 


