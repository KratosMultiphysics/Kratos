//------------------------------------------------------------------
//           _   _            _                                    .
//   KRATOS | | | |_ __  __ _| |_                                  .
//          | |_| | '  \/ _` |  _|                                 .
//           \___/|_|_|_\__,_|\__| INTERFACE                       .                             
//			                                           .
//   License:(BSD)	  UmatApplication/license.txt              .
//   Main authors:        Lluis Monforte                           .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       LlMonforte  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_UMAT_APPLICATION_H_INCLUDED )
#define  KRATOS_UMAT_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

// ConstitutiveLaw interfaces
#include "custom_laws/small_strain_umat_3D_law.hpp"
#include "custom_laws/large_strain_umat_3D_law.hpp"

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
class KratosUmatApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosUmatApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosUmatApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosUmatApplication();

	/// Destructor.
	virtual ~KratosUmatApplication(){}


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
	virtual std::string Info() const {
		return "KratosUmatApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
  		KRATOS_WATCH("in UmatApplication");
  		KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
		rOStream << "Variables:" << std::endl;
		KratosComponents<VariableData>().PrintData(rOStream);
		rOStream << std::endl;
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

	const SmallStrainUmat3DLaw  mSmallStrainUmat3DLaw;
	const LargeStrainUmat3DLaw  mLargeStrainUmat3DLaw;
	
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
	KratosUmatApplication& operator=(KratosUmatApplication const& rOther);

	/// Copy constructor.
	KratosUmatApplication(KratosUmatApplication const& rOther);


	///@}

}; // Class KratosUmatApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_UMAT_APPLICATION_H_INCLUDED  defined
