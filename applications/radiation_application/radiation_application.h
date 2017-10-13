//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


#if !defined(KRATOS_KRATOS_CONVECTION_DIFFUSIONr_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_CONVECTION_DIFFUSIONr_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/radiation_2d.h"
#include "custom_elements/radiation_3d.h"

//#include "custom_elements/conv_diff_change_of_phase_2d.h"


#include "custom_conditions/rad_face2D.h"
#include "custom_conditions/rad_face3D.h"

#include "includes/variables.h"
#include "includes/condition.h"

#include "includes/kratos_flags.h"
#include "includes/deprecated_variables.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	//KRATOS_DEFINE_VARIABLE(double,  TEMP_CONV_PROJ)	
	//KRATOS_DEFINE_VARIABLE(double,  EMISSIVITY)
			
	//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
			

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
	class KratosRadiationApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosConvectionDiffusionApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosRadiationApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosRadiationApplication();

		/// Destructor.
		virtual ~KratosRadiationApplication(){}


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
			return "KratosRadiationApplication";
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
      		KRATOS_WATCH("in KratosRadiationApplication");
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
		const Rad2D  mRad2D; 
		const Rad3D  mRad3D; 
		//const ConvDiffChangeOfPhase2D  mConvDiffChangeOfPhase2D;
		const RadFace2D  mRadFace2D;
		const RadFace3D  mRadFace3D;


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
		KratosRadiationApplication& operator=(KratosRadiationApplication const& rOther);

		/// Copy constructor.
		KratosRadiationApplication(KratosRadiationApplication const& rOther);


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


