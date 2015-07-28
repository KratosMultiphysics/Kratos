//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PFEM_FLUID_DYNAMIC_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_FLUID_DYNAMIC_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "utilities/body_normal_calculation_utils.h"
#include "processes/find_nodal_neighbours_process.h"
#include "containers/flags.h"
#include "includes/variables.h"
#include "pfem_fluid_dynamic_application_variables.h"




namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
//	KRATOS_DEFINE_VARIABLE(double, AUX_MESH_VAR )
//	KRATOS_DEFINE_VARIABLE(double, IS_INTERFACE)
//	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)
//	KRATOS_DEFINE_VARIABLE( double, NODAL_AREA )
//	KRATOS_DEFINE_VARIABLE( double, IS_STRUCTURE )
//	KRATOS_DEFINE_VARIABLE( double, IS_FLUID )
//	KRATOS_DEFINE_VARIABLE( double, IS_BOUNDARY )
//	KRATOS_DEFINE_VARIABLE( double, IS_FREE_SURFACE )

//	KRATOS_DEFINE_VARIABLE( double, NODAL_MASS )
//	KRATOS_DEFINE_VARIABLE( double, AUX_INDEX )
//	KRATOS_DEFINE_VARIABLE( double, EXTERNAL_PRESSURE )
//	KRATOS_DEFINE_VARIABLE( double, PRESSURE_OLD_IT )
//	KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_X_VALUE )
//	KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_Y_VALUE )
//	KRATOS_DEFINE_VARIABLE( double, IMPOSED_VELOCITY_Z_VALUE )
       
//	KRATOS_DEFINE_VARIABLE( double, DIVPROJ )
//      KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ADVPROJ )
//	KRATOS_DEFINE_VARIABLE( std::string, IDENTIFIER )  

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
	class KratosPFEMFluidDynamicApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosPFEMFluidDynamicApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPFEMFluidDynamicApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosPFEMFluidDynamicApplication();

		/// Destructor.
		virtual ~KratosPFEMFluidDynamicApplication(){}


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
			return "KratosPFEMFluidDynamicApplication";
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

    /// 2D instance of the VMS element


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
		KratosPFEMFluidDynamicApplication& operator=(KratosPFEMFluidDynamicApplication const& rOther);

		/// Copy constructor.
		KratosPFEMFluidDynamicApplication(KratosPFEMFluidDynamicApplication const& rOther);


		///@}    

	}; // Class KratosPFEMFluidDynamicApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_FLUID_DYNAMIC_APPLICATION_H_INCLUDED  defined 


