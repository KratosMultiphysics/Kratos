//-------------------------------------------------------------
//          ___         _           _   
//  KRATOS / __|___ _ _| |_ __ _ __| |_ 
//        | (__/ _ \ ' \  _/ _` / _|  _|
//         \___\___/_||_\__\__,_\__|\__|MECHANICS
//                                            
//  License:(BSD)    ContactMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   ...
//
//-------------------------------------------------------------
//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONTACT_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONTACT_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

// elements
#include "custom_elements/rigid_body_element.hpp"
#include "custom_elements/translatory_rigid_body_element.hpp"

// conditions
#include "custom_conditions/contact_domain_condition.hpp"
#include "custom_conditions/contact_domain_LM_2D_condition.hpp"
#include "custom_conditions/contact_domain_penalty_2D_condition.hpp"
#include "custom_conditions/axisym_contact_domain_LM_2D_condition.hpp"
#include "custom_conditions/axisym_contact_domain_penalty_2D_condition.hpp"

#include "custom_conditions/point_rigid_contact_condition.hpp"
#include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"
#include "custom_conditions/point_rigid_contact_penalty_2D_condition.hpp"
#include "custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.hpp"

// friction laws
#include "custom_friction/friction_law.hpp"
#include "custom_friction/coulomb_adhesion_friction_law.hpp"
#include "custom_friction/hardening_coulomb_friction_law.hpp"


// Core applications
#include "pfem_base_application.h"

#include "contact_mechanics_application_variables.h"

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
class KratosContactMechanicsApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{


	/// Pointer definition of KratosContactMechanicsApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosContactMechanicsApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosContactMechanicsApplication();

	/// Destructor.
	virtual ~KratosContactMechanicsApplication(){}


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
		return "KratosContactMechanicsApplication";
	}

	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const {
		rOStream << Info();
		PrintData(rOStream);
	}

	///// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const {
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

	//elements
	const RigidBodyElement                                             mRigidBodyElement;
	const TranslatoryRigidBodyElement                       mTranslatoryRigidBodyElement;

	//conditions
	const ContactDomainLM2DCondition                       mContactDomainLMCondition2D3N;
	const ContactDomainPenalty2DCondition             mContactDomainPenaltyCondition2D3N;

	const AxisymContactDomainLM2DCondition           mAxisymContactDomainLMCondition2D3N;
	const AxisymContactDomainPenalty2DCondition mAxisymContactDomainPenaltyCondition2D3N;

	const PointRigidContactPenalty2DCondition       mPointRigidContactPenalty2DCondition;
	const PointRigidContactPenalty3DCondition       mPointRigidContactPenalty3DCondition;
	const AxisymPointRigidContactPenalty2DCondition mAxisymPointRigidContactPenalty2DCondition;

	//friction laws
	const FrictionLaw                                                       mFrictionLaw;
	const CoulombAdhesionFrictionLaw                         mCoulombAdhesionFrictionLaw;
	const HardeningCoulombFrictionLaw                       mHardeningCoulombFrictionLaw;

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
	KratosContactMechanicsApplication& operator=(KratosContactMechanicsApplication const& rOther);

	/// Copy constructor.
	KratosContactMechanicsApplication(KratosContactMechanicsApplication const& rOther);


	///@}

}; // Class KratosContactMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_MECHANICS_APPLICATION_H_INCLUDED  defined
