//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED )
#define  KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED


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

	class SubdomainDisableProcess
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(SubdomainDisableProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		SubdomainDisableProcess()
		{
		}

		/// Destructor.
		virtual ~SubdomainDisableProcess()
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

		void SaveReducedPart(ModelPart& full_model_part, ModelPart& reduced_model_part)
		{
		KRATOS_TRY
		int n_int;
		int n_disabled;
		//clear reduced_model_part
		reduced_model_part.Conditions().clear();
		reduced_model_part.Elements().clear();
		reduced_model_part.Nodes().clear();

		reduced_model_part.Conditions().reserve(full_model_part.Conditions().size());
		reduced_model_part.Elements().reserve(full_model_part.Elements().size());
		reduced_model_part.Nodes().reserve(full_model_part.Nodes().size());

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
		for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; 
				im != full_model_part.ElementsEnd() ; ++im)
		{	  
			
			/*n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
			n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
			
			if (n_int==3)
				{
				im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE)=true;
				im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE)=true;
				}*/

			n_disabled=im->GetGeometry()[0].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[1].FastGetSolutionStepValue(DISABLE);
			n_disabled+=im->GetGeometry()[2].FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled<3)
				{
				reduced_model_part.Elements().push_back(*(im.base()));
				//reduced_model_part.AddNode(im->GetGeometry()[0]);				
				}
			if (n_disabled>3)
				KRATOS_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes.... " , "");

			if (n_disabled==1 || n_disabled==2)
				{
				for (int i=0;i<3;i++)
					{

					if (im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)==1)
						{
						reduced_model_part.Nodes().push_back(im->GetGeometry()(i));				

						}
	
					}
				}



		}
	
		
		for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; 
				in != full_model_part.NodesEnd() ; ++in)
		{	  
			
			n_disabled=in->FastGetSolutionStepValue(DISABLE);
			
			if (n_disabled==0.0)
				{
				reduced_model_part.Nodes().push_back(*(in.base()));
				}
			

		}

		reduced_model_part.Nodes().Unique();
	
		
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
			return "SubdomainDisableProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "SubdomainDisableProcess";
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
//		SubdomainDisableProcess& operator=(SubdomainDisableProcessconst& rOther);

		/// Copy constructor.
//		SubdomainDisableProcess(SubdomainDisableProcessconst& rOther);


		///@}    

	}; // Class SubdomainDisableProcess

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream,   
		SubdomainDisableProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const SubdomainDisableProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_SUBDOMAIN_DISABLE_PROCESS_INCLUDED  defined 


