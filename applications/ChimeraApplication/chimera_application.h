//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================
#if !defined(KRATOS_CHIMERA_APPLICATION_H_INCLUDED )
#define  KRATOS_CHIMERA_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

#include "custom_elements/sksy_fluid_element.h"
#include "custom_conditions/sksy_fluid_condition.h"
#include "custom_conditions/chimera_fluid_coupling_condition.h"
#include "custom_conditions/chimera_thermal_coupling_condition.h"


namespace Kratos {

///@name Kratos Globals
///@{

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
class KratosChimeraApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosChimeraApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosChimeraApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosChimeraApplication();

	/// Destructor.
	virtual ~KratosChimeraApplication(){}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{

	virtual void Register() override;



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
	virtual std::string Info() const override{
		return "KratosChimeraApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const override{
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const override{
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

	// static const ApplicationCondition  msApplicationCondition;

	///@}
	///@name Member Variables
	///@{

	// Skew-symmetric fluid element
	const SkSyFluidElement<2> mSkSyFluidElement2D3N;
	const SkSyFluidElement<3> mSkSyFluidElement3D4N;

	const SkSyFluidCondition<2> mSkSyFluidCondition2D2N;
	const SkSyFluidCondition<3> mSkSyFluidCondition3D3N;

	// Old Chimera Neumann conditions
    const ChimeraFluidCouplingCondition<2> mChimeraFluidCouplingCondition2D;
    const ChimeraFluidCouplingCondition<3> mChimeraFluidCouplingCondition3D;
    const ChimeraThermalCouplingCondition<2> mChimeraThermalCouplingCondition2D;
    const ChimeraThermalCouplingCondition<3> mChimeraThermalCouplingCondition3D;


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
	KratosChimeraApplication& operator=(KratosChimeraApplication const& rOther);

	/// Copy constructor.
	KratosChimeraApplication(KratosChimeraApplication const& rOther);


	///@}

}; // Class KratosChimeraApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CHIMERA_APPLICATION_H_INCLUDED  defined
