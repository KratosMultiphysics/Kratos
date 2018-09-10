//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_EMPIRE_APPLICATION_H_INCLUDED )
#define  KRATOS_EMPIRE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"


namespace Kratos
{

	///@name Kratos Globals
	///@{

	// Variables definition
//	KRATOS_DEFINE_VARIABLE(double, AUX_MESH_VAR )
//	KRATOS_DEFINE_VARIABLE(double, IS_INTERFACE)
//	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


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
	class KratosEmpireApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{


		/// Pointer definition of KratosEmpireApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosEmpireApplication);

		///@}
		///@name Life Cycle
		///@{

		/// Default constructor.
		KratosEmpireApplication();

		/// Destructor.
		~KratosEmpireApplication() override {}


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
		std::string Info() const override
		{
			return "KratosEmpireApplication";
		}

		/// Print information about this object.
		void PrintInfo(std::ostream& rOStream) const override
		{
			rOStream << Info();
			PrintData(rOStream);
		}

		///// Print object's data.
      void PrintData(std::ostream& rOStream) const override
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
// 		const Elem2D   mElem2D;
// 		const Elem3D   mElem3D;


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
		KratosEmpireApplication& operator=(KratosEmpireApplication const& rOther);

		/// Copy constructor.
		KratosEmpireApplication(KratosEmpireApplication const& rOther);


		///@}

	}; // Class KratosEmpireApplication

	///@}


	///@name Type Definitions
	///@{


	///@}
	///@name Input and output
	///@{

	///@}


}  // namespace Kratos.

#endif // KRATOS_EMPIRE_APPLICATION_H_INCLUDED  defined


