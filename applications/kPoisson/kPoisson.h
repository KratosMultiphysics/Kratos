//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_POISSON_H_INCLUDED )
#define  KRATOS_POISSON_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/poisson_2d.h"

#include "includes/variables.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 

	// for Poisson application

	KRATOS_DEFINE_VARIABLE(double, DUMMY_UNKNOWN)
	KRATOS_DEFINE_VARIABLE(double, DUMMY_MATERIAL)
	KRATOS_DEFINE_VARIABLE(double, DUMMY_POINT_SOURCE)

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
	*/
	class KratosR1PoissonApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosR1PoissonApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosR1PoissonApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosR1PoissonApplication();

		/// Destructor.
		virtual ~KratosR1PoissonApplication(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{

		virtual void Register();



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
			return "KratosR1PoissonApplication";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
      	KRATOS_WATCH("in KratosR1PoissonApplication");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Elements:" << std::endl;
		KratosComponents<Element>().PrintData(rOStream);
		rOStream << std::endl;
		rOStream << "Conditions:" << std::endl;
		KratosComponents<Condition>().PrintData(rOStream);
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

		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
		///@name Member Variables 
		///@{ 
		const Poisson2D  mPoisson2D;

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
		KratosR1PoissonApplication& operator=(KratosR1PoissonApplication const& rOther);

		/// Copy constructor.
		KratosR1PoissonApplication(KratosR1PoissonApplication const& rOther);


		///@}    

	}; // Class KratosR1PoissonApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.


#endif // KRATOS_POISSON_H_INCLUDED  defined 


