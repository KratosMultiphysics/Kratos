//------------------------------------------------------------------
//           _   _            _                                    .
//   KRATOS | | | |_ __  __ _| |_                                  .
//          | |_| | '  \/ _` |  _|                                 .
//           \___/|_|_|_\__,_|\__| INTERFACE                       .                             
//			                                           .
//   License:(BSD)	  UmatApplication/license.txt              .
//   Main authors:        LlMonforte, JMCarbonell                  .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       JMCarbonell $
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

// Constitutive models
#include "custom_models/von_mises_umat_small_strain_model.hpp"
#include "custom_models/hypoplastic_umat_small_strain_model.hpp"
#include "custom_models/von_mises_umat_large_strain_model.hpp"

// Core applications
#include "umat_application_variables.h"

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
		return "KratosUmatApplication";
	}

	/// Print information about this object.
	void PrintInfo(std::ostream& rOStream) const override {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	void PrintData(std::ostream& rOStream) const override {
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

   const VonMisesSmallStrainUmatModel mVonMisesSmallStrainUmatModel;
   const HypoplasticSmallStrainUmatModel mHypoplasticSmallStrainUmatModel;
   const VonMisesLargeStrainUmatModel mVonMisesLargeStrainUmatModel;
	
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
