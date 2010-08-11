//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmora $ rrossi $
//   Date:                $Date: 2010-02-02 $
//   Revision:            $Revision: 1.4 $
//
//


#if !defined(KRATOS_MAGNETOSTATIC_H_INCLUDED )
#define  KRATOS_MAGNETOSTATIC_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// External includes 

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/magnetostatic_2d.h"
#include "custom_elements/magnetostatic_3d.h"
#include "custom_conditions/pointcurrent2D.h"
#include "custom_conditions/mfield2D.h"
//#include "custom_conditions/mfield3D.h"

#include "custom_elements/conductivity_3d.h"

#include "includes/variables.h"
#include "includes/condition.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
	///@name Kratos Globals
	///@{ 

	// Variables definition 
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
	// for kMagnetostatic application
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_PERMEABILITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(COERCIVITY)
	KRATOS_DEFINE_VARIABLE(double, MAGNETOSTATIC_POTENTIAL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VECTOR_POTENTIAL)
	KRATOS_DEFINE_VARIABLE(double, MAGNETOSTATIC_POINT_CURRENT)
	KRATOS_DEFINE_VARIABLE(double, MAGNETOSTATIC_SURFACE_CURRENT)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAGNETOSTATIC_VOLUME_CURRENT)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FIELD_INTENSITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(MAGNETIC_FLUX_DENSITY)
	KRATOS_DEFINE_VARIABLE(double, INFINIT_COEFFICIENT)

	KRATOS_DEFINE_VARIABLE(double, ELECTROSTATIC_POTENTIAL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELECTRICAL_CONDUCTIVITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_FIELD)

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
	class KratosR1MagnetostaticApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosR1MagnetostaticApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosR1MagnetostaticApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosR1MagnetostaticApplication();

		/// Destructor.
		virtual ~KratosR1MagnetostaticApplication(){}

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
			return "KratosR1MagnetostaticApplication";
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
      	KRATOS_WATCH("in KratosR1MagnetostaticApplication");
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
		const Magnetostatic2D  mMagnetostatic2D;
		const Magnetostatic3D  mMagnetostatic3D;
		const PointCurrent2D  mPointCurrent2D;
		const Mfield2D  mMfield2D;
		//const Mfield3D  mMfield3D;

		//const Conductivity2D mConductivity2D;
		const Conductivity3D mConductivity3D;

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
		KratosR1MagnetostaticApplication& operator=(KratosR1MagnetostaticApplication const& rOther);

		/// Copy constructor.
		KratosR1MagnetostaticApplication(KratosR1MagnetostaticApplication const& rOther);

		///@}    

	}; // Class KratosR1MagnetostaticApplication 

	///@} 

	///@name Type Definitions       
	///@{ 

	///@} 
	///@name Input and output 
	///@{ 

	///@} 

}  // namespace Kratos.


#endif // KRATOS_MAGNETOSTATIC_H_INCLUDED  defined 


