//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-23 15:38:00 $
//   Revision:            $Revision: 1.1 $
//
//		
// This process creates a pressue_force condition at each node

#if !defined(KRATOS_ASSIGN_PRESSURE_FORCE_PROCESS_INCLUDED )
#define  KRATOS_ASSIGN_PRESSURE_FORCE_PROCESS_INCLUDED



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
#include "custom_conditions/pressure_cond.h"
#include "custom_conditions/pressure_cond3d.h"


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

	class AssignPressureForceProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PressureContributionProcess
		KRATOS_CLASS_POINTER_DEFINITION(AssignPressureForceProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		AssignPressureForceProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~AssignPressureForceProcess()
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
			//first erase all the conditions
			mr_model_part.Conditions().clear();
			unsigned int dim = mr_model_part.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
									
			ModelPart::ConditionsContainerType& rConds = mr_model_part.Conditions();

			int end_id = 0;
			if (rConds.size()!=0)
				end_id = (rConds.end()-1)->Id();

			
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
			{
				if(in->GetValue(NEIGHBOUR_ELEMENTS).size() != 0)
				{
					Condition::NodesArrayType temp;
					temp.reserve(1);
		//			temp.push_back(*in); 
					temp.push_back(*(in.base())); 

					Geometry< Node<3> >::Pointer cond = 
								Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );

		//if(cond.size() != 1)
		//  KRATOS_WATCH("aaaaaaaaahhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh")

					int id = end_id+1;
						//in->Id();
							//
					//CHECK OUT IF THAT IS CORRECT!!!! (just the next line)
					//PressureCond::Pointer p_cond = PressureCond::Pointer(new PressureCond(id, cond) );
					if (dim==2)
					{
					Condition::Pointer p_cond = Condition::Pointer(new PressureCond(id, cond) );
								//Condition::Pointer(new Condition(id, cond) );
					//p_cond.
					end_id++;
					mr_model_part.Conditions().push_back(p_cond);
					}
					if (dim==3)
					{
					Condition::Pointer p_cond = Condition::Pointer(new PressureCond3D(id, cond) );
								//Condition::Pointer(new Condition(id, cond) );
					//p_cond.
					end_id++;
					mr_model_part.Conditions().push_back(p_cond);
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
			return "AssignPressureForceProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "AssignPressureForceProcess";
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
//		AssignPressureForceProcess& operator=(AssignPressureForceProcess const& rOther);

		/// Copy constructor.
//		AssignPressureForceProcess(AssignPressureForceProcess const& rOther);


		///@}    

	}; // Class AssignPressureForceProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		AssignPressureForceProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const AssignPressureForceProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_ASSIGN_PRESSURE_FORCE_PROCESS_INCLUDED  defined 


