//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmarti $
//   Date:                $Date: 2009-01-23 14:27:34 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_QCOMPRESSIBLE_FLUID_APPLICATION_H_INCLUDED )
#define KRATOS_QCOMPRESSIBLE_FLUID_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/Qfluid_3d.h"
#include "custom_elements/Qfluid_2d.h"


//#include "custom_conditions/fluid3d_neumann.h"

#include "includes/variables.h"
#include "includes/condition.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DESP)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)
	KRATOS_DEFINE_VARIABLE( double, MACH_NUMBER )
	KRATOS_DEFINE_VARIABLE( double, PRESSURE_COEFFICIENT )
	KRATOS_DEFINE_VARIABLE( double, PRESSURE_OLD_IT )
	KRATOS_DEFINE_VARIABLE( double, PRESSUREAUX_OLD_IT )
	KRATOS_DEFINE_VARIABLE( Vector, BDF_COEFFICIENTS )
	KRATOS_DEFINE_VARIABLE(double, NODAL_MASS)
        KRATOS_DEFINE_VARIABLE(double, NODAL_MASSAUX)
	KRATOS_DEFINE_VARIABLE(double, MASS)
	//KRATOS_DEFINE_VARIABLE(double, NODAL_MAUX)
	//KRATOS_DEFINE_VARIABLE(double, NODAL_PAUX)
	KRATOS_DEFINE_VARIABLE(double, NODAL_PRESS)
	KRATOS_DEFINE_VARIABLE(double, NODAL_PRESSAUX)
	KRATOS_DEFINE_VARIABLE(int, AUX_INDEX)
	KRATOS_DEFINE_VARIABLE(double, EXTERNAL_PRESSURE)
	KRATOS_DEFINE_VARIABLE(double, DIAMETER)
	KRATOS_DEFINE_VARIABLE(double, PERMEABILITY_INV)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
			
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
	class KratosQcompressibleFluidApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosR1IncompressibleFluidApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosQcompressibleFluidApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosQcompressibleFluidApplication();

		/// Destructor.
		virtual ~KratosQcompressibleFluidApplication(){}


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
			return "KratosQcompressibleFluidApplication";
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
      	KRATOS_WATCH("in KratosQcompressibleFluidApplication");
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


		///@} 
		///@name Private Operators
		///@{ 
		const QFluid3D  mQFluid3D; 
		const QFluid2D  mQFluid2D; 

	
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
		KratosQcompressibleFluidApplication& operator=(KratosQcompressibleFluidApplication const& rOther);

		/// Copy constructor.
		KratosQcompressibleFluidApplication(KratosQcompressibleFluidApplication const& rOther);


		///@}    

	}; // Class KratosR1IncompressibleFluidApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED  defined 


