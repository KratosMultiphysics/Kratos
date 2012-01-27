//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_THERMOMECHANICAL_APPLICATION_H_INCLUDED )
#define  KRATOS_THERMOMECHANICAL_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
#include "custom_elements/heat_contact_2d.h"
#include "custom_elements/heat_contact_3d.h"
#include "custom_elements/SUPG_conv_diff_2d.h"
#include "custom_elements/SUPG_conv_diff_3d.h"
#include "custom_elements/SUPG_conv_3d.h"
#include "custom_elements/thermal_face2d.h"
#include "custom_elements/thermal_face3d.h"
#include "custom_elements/environment_contact.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
        KRATOS_DEFINE_VARIABLE(int, NODE_PROPERTY_ID)
	KRATOS_DEFINE_VARIABLE(double,  AMBIENT_TEMPERATURE)	
	KRATOS_DEFINE_VARIABLE(double,  HTC)
        KRATOS_DEFINE_VARIABLE(int, REF_ID)
	
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
	class KratosThermoMechanicalApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosThermoMechanicalApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosThermoMechanicalApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosThermoMechanicalApplication();

		/// Destructor.
		virtual ~KratosThermoMechanicalApplication(){}


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
			return "KratosThermoMechanicalApplication";
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
 		const HeatContact2D  mHeatContact2D;
 		const HeatContact3D  mHeatContact3D;
		const ThermalFace2D  mThermalFace2D;
		const ThermalFace3D  mThermalFace3D;
		const EnvironmentContact  mEnvironmentContact;
		
 		const SUPGConvDiff2D  mSUPGConvDiff2D;		
 		const SUPGConvDiff3D  mSUPGConvDiff3D;	
 		const SUPGConv3D  mSUPGConv3D;		
		
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
		KratosThermoMechanicalApplication& operator=(KratosThermoMechanicalApplication const& rOther);

		/// Copy constructor.
		KratosThermoMechanicalApplication(KratosThermoMechanicalApplication const& rOther);


		///@}    

	}; // Class KratosThermoMechanicalApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_THERMOMECHANICAL_APPLICATION_H_INCLUDED  defined 


