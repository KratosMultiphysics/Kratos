//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-12-15 15:41:36 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/conv_diff_2d.h"
#include "custom_elements/conv_diff_3d.h"

#include "custom_elements/conv_diff_change_of_phase_2d.h"


#include "custom_conditions/thermal_face2D.h"
#include "custom_conditions/thermal_face3D.h"

#include "includes/variables.h"
#include "includes/condition.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
//	KRATOS_DEFINE_VARIABLE( Vector, BDF_COEFFICIENTS )
	//KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)
//	KRATOS_DEFINE_VARIABLE(int, AUX_INDEX)
	KRATOS_DEFINE_VARIABLE(double,  CONDUCTIVITY)
	KRATOS_DEFINE_VARIABLE(double,  SPECIFIC_HEAT)
	KRATOS_DEFINE_VARIABLE(double,  HEAT_FLUX)	
	KRATOS_DEFINE_VARIABLE(double,  TEMP_CONV_PROJ)	

	//Added by Pavel and Annelie	
	KRATOS_DEFINE_VARIABLE(double,  ENTHALPY)
	KRATOS_DEFINE_VARIABLE(double,  LATENT_HEAT)
	KRATOS_DEFINE_VARIABLE(double,  MELT_TEMPERATURE_1)
	KRATOS_DEFINE_VARIABLE(double,  MELT_TEMPERATURE_2)

	KRATOS_DEFINE_VARIABLE(double,  AMBIENT_TEMPERATURE)	
	KRATOS_DEFINE_VARIABLE(double,  EMISSIVITY)	
	KRATOS_DEFINE_VARIABLE(double,  CONVECTION_COEFFICIENT)	
	KRATOS_DEFINE_VARIABLE(double,  FACE_HEAT_FLUX)	
			
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
			

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
	class KratosConvectionDiffusionApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosConvectionDiffusionApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosConvectionDiffusionApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosConvectionDiffusionApplication();

		/// Destructor.
		virtual ~KratosConvectionDiffusionApplication(){}


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
			return "KratosConvectionDiffusionApplication";
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
      	KRATOS_WATCH("in KratosConvectionDiffusionApplication");
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

//		static const ConvDiff2D  msConvDiff2D; 
//		static const ConvDiff3D  msConvDiff3D; 
//		static const ConvDiff2DChangeOfPhase  msCConvDiff2DChangeOfPhase;
//		static const ThermalFace2D  msThermalFace2D;
//		static const ThermalFace3D  msThermalFace3D;


		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
		///@name Member Variables 
		///@{ 
		const ConvDiff2D  mConvDiff2D; 
		const ConvDiff3D  mConvDiff3D; 
		const ConvDiffChangeOfPhase2D  mConvDiffChangeOfPhase2D;
		const ThermalFace2D  mThermalFace2D;
		const ThermalFace3D  mThermalFace3D;


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
		KratosConvectionDiffusionApplication& operator=(KratosConvectionDiffusionApplication const& rOther);

		/// Copy constructor.
		KratosConvectionDiffusionApplication(KratosConvectionDiffusionApplication const& rOther);


		///@}    

	}; // Class KratosConvectionDiffusionApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED  defined 


