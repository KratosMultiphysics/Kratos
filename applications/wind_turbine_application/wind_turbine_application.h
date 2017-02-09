//   
//   Project Name:        Kratos Wind Turbine Application
//   Last Modified by:    $Author: efusto $
//   Date:                $Date: 20120507 18:08:58 $
//   Revision:            $Revision: 0.1 $
//
//


#if !defined(KRATOS_WIND_TURBINE_APPLICATION_H_INCLUDED )
#define  KRATOS_WIND_TURBINE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

namespace Kratos
{
	///@name Kratos Globals
	///@{ 

        // Variables definition
        KRATOS_DEFINE_VARIABLE(int, AUX_ID);
        KRATOS_DEFINE_VARIABLE(int, AUX_INTERF_ID );
        KRATOS_DEFINE_VARIABLE(int, AUX_BASE_ID );
//	KRATOS_DEFINE_VARIABLE(double, AUX_MESH_VAR )
//	KRATOS_DEFINE_VARIABLE(double, IS_INTERFACE)
//	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


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
	class KratosWindTurbineApplication : public KratosApplication
	{
	
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of KratosWindTurbineApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosWindTurbineApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor
		KratosWindTurbineApplication();

		/// Destructor
		virtual ~KratosWindTurbineApplication(){}


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
			return "KratosWindTurbineApplication";
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
      	KRATOS_WATCH("in my application");
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
// 		const Elem2D   mElem2D; 
// 		const Elem3D   mElem3D; 


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
		KratosWindTurbineApplication& operator=(KratosWindTurbineApplication const& rOther);

		/// Copy constructor.
		KratosWindTurbineApplication(KratosWindTurbineApplication const& rOther);

		///@}    

	}; // Class KratosWindTurbineApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_WIND_TURBINE_APPLICATION_H_INCLUDED  defined 


