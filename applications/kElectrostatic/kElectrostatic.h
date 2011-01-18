//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmora $ rrossi $
//   Date:                $Date: 2009-09-28 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_ELECTROSTATIC_H_INCLUDED )
#define  KRATOS_ELECTROSTATIC_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/electrostatic_2d.h"
#include "custom_elements/electrostatic_3d.h"
#include "custom_conditions/pointcharge2D.h"
#include "custom_conditions/pointcharge3D.h"
#include "custom_conditions/efield2D.h"
#include "custom_conditions/efield3D.h"

//#include "custom_utilities/custom_gid_io.h"

#include "includes/variables.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"


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

	KRATOS_DEFINE_VARIABLE(double,  AMBIENT_TEMPERATURE)	
	KRATOS_DEFINE_VARIABLE(double,  EMISSIVITY)	
	KRATOS_DEFINE_VARIABLE(double,  CONVECTION_COEFFICIENT)	
	KRATOS_DEFINE_VARIABLE(double,  FACE_HEAT_FLUX)	
	
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
			
	// for electromagnetic applications
	// for kElectrostatic application
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELECTRICAL_PERMITTIVITY)
	KRATOS_DEFINE_VARIABLE(double, ELECTROSTATIC_POTENTIAL)
	KRATOS_DEFINE_VARIABLE(double, ELECTROSTATIC_POINT_CHARGE)
	KRATOS_DEFINE_VARIABLE(double, ELECTROSTATIC_SURFACE_CHARGE)
	KRATOS_DEFINE_VARIABLE(double, ELECTROSTATIC_VOLUME_CHARGE)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_FIELD)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_DISPLACEMENT_FIELD)
	KRATOS_DEFINE_VARIABLE(double, INFINIT_COEFFICIENT)

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
	class KratosR1ElectrostaticApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosR1ElectrostaticApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosR1ElectrostaticApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosR1ElectrostaticApplication();

		/// Destructor.
		virtual ~KratosR1ElectrostaticApplication(){}


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
			return "KratosR1ElectrostaticApplication";
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
      	KRATOS_WATCH("in KratosR1ElectrostaticApplication");
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
//		static const ThermalFace2D  msThermalFace2D;
//		static const ThermalFace3D  msThermalFace3D;


		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
		///@name Member Variables 
		///@{ 
		const Electrostatic2D  mElectrostatic2D;
		const Electrostatic3D  mElectrostatic3D;
		const PointCharge2D  mPointCharge2D;
		const PointCharge3D  mPointCharge3D;
		const Efield2D  mEfield2D;
		const Efield3D  mEfield3D;

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
		KratosR1ElectrostaticApplication& operator=(KratosR1ElectrostaticApplication const& rOther);

		/// Copy constructor.
		KratosR1ElectrostaticApplication(KratosR1ElectrostaticApplication const& rOther);


		///@}    

	}; // Class KratosR1ElectrostaticApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.


#endif // KRATOS_ELECTROSTATIC_H_INCLUDED  defined 


