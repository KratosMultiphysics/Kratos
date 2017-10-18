//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_FLUIDRVELAGRANGEMULTIPLIERS_APPLICATION_H_INCLUDED )
#define  KRATOS_FLUIDRVELAGRANGEMULTIPLIERS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"

#include "custom_conditions/lagrange_multiplier_mean_velocity_2d.h" //the condition
#include "custom_conditions/lagrange_multiplier_velocity_gradients_2d.h" //the condition2
#include "custom_conditions/inverse_tangent_velocity_periodic_condition_2d.h" //the condition3
#include "custom_conditions/inverse_normal_velocity_periodic_condition_2d.h" //the condition3

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS)
	
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_X_COMPONENT )
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_Y_COMPONENT )
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_Z_COMPONENT )
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_PRESSURE )
	
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_X_COMPONENT_ANTIPERIODIC )
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_Y_COMPONENT_ANTIPERIODIC )
	KRATOS_DEFINE_VARIABLE(int, NODE_PAIR_Z_COMPONENT_ANTIPERIODIC )	
	KRATOS_DEFINE_VARIABLE(double, LAGRANGE_MULTIPLIER_TANGENT_VELOCITY )
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_MULTIPLIER_NORMAL_VELOCITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(REYNOLDS_STRESS_2D)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_TILDA_TILDA_STRESS_2D)

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
	class KratosFluidRveLagrangeMultipliersApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosFluidRveLagrangeMultipliersApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosFluidRveLagrangeMultipliersApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosFluidRveLagrangeMultipliersApplication();

		/// Destructor.
		virtual ~KratosFluidRveLagrangeMultipliersApplication(){}


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
			return "KratosFluidRveLagrangeMultipliersApplication";
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
		const MeanVelocityLagrangeMultiplierCondition2D   mMeanVelocityLagrangeMultiplierCondition2D; 
		const VelocityGradientsLagrangeMultiplierCondition2D   mVelocityGradientsLagrangeMultiplierCondition2D; 
		const InverseTangentVelocityPeriodicCondition2D2N   mInverseTangentVelocityPeriodicCondition2D2N; 
		const InverseNormalVelocityPeriodicCondition2D2N   mInverseNormalVelocityPeriodicCondition2D2N; 

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
		KratosFluidRveLagrangeMultipliersApplication& operator=(KratosFluidRveLagrangeMultipliersApplication const& rOther);

		/// Copy constructor.
		KratosFluidRveLagrangeMultipliersApplication(KratosFluidRveLagrangeMultipliersApplication const& rOther);


		///@}    

	}; // Class KratosFluidRveLagrangeMultipliersApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_FLUIDRVELAGRANGEMULTIPLIERS_APPLICATION_H_INCLUDED  defined 


