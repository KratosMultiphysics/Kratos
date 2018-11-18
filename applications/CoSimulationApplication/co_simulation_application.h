//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED )
#define  KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


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
class KratosCoSimulationApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosCoSimulationApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosCoSimulationApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosCoSimulationApplication();

	/// Destructor.
	virtual ~KratosCoSimulationApplication() override {}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
	///@{
	///@{
  		void Register() override;
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
	virtual std::string Info() const override {
		return "KratosCoSimulationApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const override {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const override {
			rOStream << "in KratosCoSimulationApplication" << std::endl;
			KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
			KratosApplication::PrintData(rOStream);
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
	KratosCoSimulationApplication& operator=(KratosCoSimulationApplication const& rOther);

	/// Copy constructor.
	KratosCoSimulationApplication(KratosCoSimulationApplication const& rOther);


	///@}

}; // Class KratosCoSimulationApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CO_SIMULATION_APPLICATION_H_INCLUDED  defined
