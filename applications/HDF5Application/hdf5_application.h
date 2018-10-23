//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main authors:    Jordi Cotela
//


#if !defined(KRATOS_HDF5_APPLICATION_H_INCLUDED )
#define  KRATOS_HDF5_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"


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

/// This application provides an interface to write Kratos data using HDF5.
class KratosHDF5Application : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosHDF5Application
	KRATOS_CLASS_POINTER_DEFINITION(KratosHDF5Application);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosHDF5Application();

	/// Destructor.
	~KratosHDF5Application() override {}


	///@}
	///@name Operators
	///@{


	///@}
	///@name Operations
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
	std::string Info() const override {
		return "KratosHDF5Application";
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	void PrintData(std::ostream& rOStream) const override
	{
		rOStream << "In KratosHDF5Application" << std::endl;

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
	KratosHDF5Application& operator=(KratosHDF5Application const& rOther);

	/// Copy constructor.
	KratosHDF5Application(KratosHDF5Application const& rOther);


	///@}

}; // Class KratosHDF5Application

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_HDF5_APPLICATION_H_INCLUDED  defined
