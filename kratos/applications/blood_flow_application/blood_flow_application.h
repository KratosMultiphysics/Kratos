//
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_BLOOD_FLOW_APPLICATION_H_INCLUDED )
#define  KRATOS_BLOOD_FLOW_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "includes/variables.h"

#include "custom_elements/artery_element.h"
#include "custom_conditions/artery_11_condition.h"
#include "custom_conditions/artery_1d_to_3d_condition.h"
#include "custom_conditions/artery_3d_to_1d_condition.h"
#include "custom_conditions/artery_12_condition.h"
#include "custom_conditions/artery_inlet_condition.h"
#include "custom_conditions/artery_inlet_condition_pressure.h"
#include "custom_conditions/artery_outlet_condition.h"
#include "custom_conditions/artery_outlet_condition_free.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	KRATOS_DEFINE_VARIABLE(double, FLOW )
    KRATOS_DEFINE_VARIABLE(double, PRESSURE_VENOUS)
    KRATOS_DEFINE_VARIABLE(double, TERMINAL_RESISTANCE)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(WORK)
	KRATOS_DEFINE_VARIABLE(double, BETA)
    KRATOS_DEFINE_VARIABLE(double, C0)
    KRATOS_DEFINE_VARIABLE(double, DYASTOLIC_PRESSURE)
    KRATOS_DEFINE_VARIABLE(double, SYSTOLIC_PRESSURE)
    KRATOS_DEFINE_VARIABLE(double, AVERAGE_PRESSURE)
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
	class KratosBloodFlowApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosBloodFlowApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosBloodFlowApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosBloodFlowApplication();

		/// Destructor.
		virtual ~KratosBloodFlowApplication(){}


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
			return "KratosBloodFlowApplication";
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

		///@} 
		///@name Member Variables 
		///@{ 
        const ArteryElement mArteryElement;
        const Artery11Condition mArtery11Condition;
        const Artery1Dto3DCondition mArtery1Dto3DCondition;
        const Artery3Dto1DCondition mArtery3Dto1DCondition;
        const Artery12Condition mArtery12Condition;
        const ArteryInletCondition mArteryInletCondition;
        const ArteryInletConditionPressure mArteryInletConditionPressure;
        const ArteryOutletCondition mArteryOutletCondition;
        const ArteryOutletFreeCondition mArteryOutletFreeCondition;

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
		KratosBloodFlowApplication& operator=(KratosBloodFlowApplication const& rOther);

		/// Copy constructor.
        KratosBloodFlowApplication(KratosBloodFlowApplication const& rOther);


		///@}    

	}; // Class KratosBloodFlowApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_BLOOD_FLOW_APPLICATION_H_INCLUDED  defined 


