//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-05 09:39:14 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "includes/variables.h"
#include "includes/condition.h"
//#include "custom_elements/updated_lagrangian_fluid.h"
#include "custom_conditions/free_surface_cond2d.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
/*	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)
	KRATOS_DEFINE_VARIABLE(double, NODAL_H)
	KRATOS_DEFINE_VARIABLE(double, IS_STRUCTURE)
	KRATOS_DEFINE_VARIABLE(double, IS_FLUID)
	KRATOS_DEFINE_VARIABLE(double, IS_BOUNDARY)
	KRATOS_DEFINE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_DEFINE_VARIABLE(double, IS_FREE_SURFACE)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(NORMAL_TO_WALL)
*/
//	KRATOS_DEFINE_VARIABLE(double, PRESSURE_OLD_IT)
//	KRATOS_DEFINE_VARIABLE(double, NODAL_MASS)
//	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VAUX)

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
	class KratosPFEMApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosPFEMApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosPFEMApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosPFEMApplication();

		/// Destructor.
		virtual ~KratosPFEMApplication(){}


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
			return "KratosPFEMApplication";
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
      	KRATOS_WATCH("in KratosPFEMApplication");
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
		//const UpdatedLagrangianFluid mUpdatedLagrangianFluid2D;
		const FreeSurfaceCond2d mFreeSurfaceCond2d;



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
		KratosPFEMApplication& operator=(KratosPFEMApplication const& rOther);

		/// Copy constructor.
		KratosPFEMApplication(KratosPFEMApplication const& rOther);


		///@}    

	}; // Class KratosPFEMApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED  defined 


