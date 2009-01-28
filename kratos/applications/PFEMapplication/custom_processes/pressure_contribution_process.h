//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-23 15:38:00 $
//   Revision:            $Revision: 1.1 $ 20 October 2006 CLEAN VERSION
//
//		
// THIS PROCESS is INVENTED in order to CALCULATE THE VOLUMETRIC force and
//STORE it node-wise. that is the term Integral(BdeltaSvolumetric).
// this term will be later added to the RHS vector as a prestress!!!
// this is done inside the pressure_contribution.cpp

#if !defined(KRATOS_PRESSURE_CONTRIBUTION_PROCESS_INCLUDED )
#define  KRATOS_PRESSURE_CONTRIBUTION_PROCESS_INCLUDED



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
#include "utilities/geometry_utilities.h"
#include "custom_elements/updated_lagrangian_fluid.h"
#include "PFEM_application.h"


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

	class PressureContributionProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PressureContributionProcess
		KRATOS_CLASS_POINTER_DEFINITION(PressureContributionProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		PressureContributionProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~PressureContributionProcess()
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
			KRATOS_TRY
					
			ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
			array_1d<double,3> dummy;
			//first initialize the pressure force to the old value
			
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ; 
					in != mr_model_part.NodesEnd() ; ++in)
			{   
				noalias(in->FastGetSolutionStepValue(PRESSURE_FORCE))
					= in->FastGetSolutionStepValue(PRESSURE_FORCE,1);
			}
			
			//and now add the increment due to the current solution step displacement
			for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
				im != mr_model_part.ElementsEnd() ; ++im)
			{  
				im->Calculate(PRESSURE_FORCE, dummy, proc_info);
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
			return "PressureContributionProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "PressureContributionProcess";
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
		double m_min_h;


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
//		PressureContributionProcess& operator=(PressureContributionProcess const& rOther);

		/// Copy constructor.
//		PressureContributionProcess(PressureContributionProcess const& rOther);


		///@}    

	}; // Class PressureContributionProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		PressureContributionProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const PressureContributionProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_PRESSURE_CONTRIBUTION_PROCESS_INCLUDED  defined 


