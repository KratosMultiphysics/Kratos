//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_CHECK_FSI_CONVERGENCE_INCLUDED )
#define  KRATOS_CHECK_FSI_CONVERGENCE_INCLUDED



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
		calculate the nodal H for all the nodes depending on the min distance
		of the neighbouring nodes.

		lonely nodes are given the average value of the H
	*/

	class CheckFSIConvergence 
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of CheckFSIConvergence
		KRATOS_CLASS_POINTER_DEFINITION(CheckFSIConvergence);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		CheckFSIConvergence()
		{
		}

		/// Destructor.
		virtual ~CheckFSIConvergence()
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

		bool CheckConvergence(PointerVector< Node<3> >& fluid_interface_nodes, double toll)
		{
			KRATOS_TRY

			double dp_norm = 0.00;
			double p_norm = 0.00;
			double dp,p;
			//verifies that the pressure variation is sufficiently small
			for(PointerVector< Node<3> >::iterator it = fluid_interface_nodes.begin();
												it != fluid_interface_nodes.end(); it++);
			{
				p = it->FastGetSolutionStepValue(PRESSURE);
				dp = p - it->FastGetSolutionStepValue(PRESSURE,1);

				dp_norm += dp*dp;
				p_norm += p*p
			}

			if( p_norm = 0.00)
				p_norm = 1.00;

			double ratio = sqrt(p_norm/dp_norm);
			std::cout << "FSI convergence ratio = " << ratio << std::endl;

			if(ratio < toll)
				return true;
			else
				return false;

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
			return "CheckFSIConvergence";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "CheckFSIConvergence";
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
		CheckFSIConvergence& operator=(CheckFSIConvergence const& rOther);

		/// Copy constructor.
		//CheckFSIConvergence(CheckFSIConvergence const& rOther);


		///@}    

	}; // Class CheckFSIConvergence 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		CheckFSIConvergence& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const CheckFSIConvergence& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_CHECK_FSI_CONVERGENCE_INCLUDED  defined 


