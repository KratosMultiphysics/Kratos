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


#if !defined(KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"

#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "includes/model_part.h"
#include "solid_mechanics_application.h"
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

/* ELEMENTS */
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/meshless_membrane_element.h"
#include "custom_elements/meshless_laplace_element.h"
#include "custom_elements/meshless_shell_element.h"

//#include "custom_conditions/ContinuityConditionLagrange.h"
//#include "custom_conditions/ContinuityConditionPenalty.h"
//#include "custom_conditions/LoadCondition.h"
//#include "custom_conditions/SupportCondition.h"

/*CONDITIONS*/
#include "custom_conditions/meshless_support_rotation_condition.h"
#include "custom_conditions/meshless_load_condition.h"
#include "custom_conditions/meshless_lagrange_coupling_condition.h"
#include "custom_conditions/meshless_penalty_coupling_rotation_condition.h"

namespace Kratos {

///@name Kratos Globals
///@{

	//Variables definition
	KRATOS_DEFINE_VARIABLE(double, INTEGRATION_WEIGHT)
	KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)

	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)

	KRATOS_DEFINE_VARIABLE(Vector, TANGENTS)

	KRATOS_DEFINE_VARIABLE(double, PENALTY_FACTOR)

	KRATOS_DEFINE_VARIABLE(int, DISPLACEMENT_ROTATION_FIX)
	KRATOS_DEFINE_VARIABLE(int, LOAD_TYPE)
	//KRATOS_DEFINE_VARIABLE(std::string, CONDITION_TYPE_DEFINITION)
	KRATOS_DEFINE_VARIABLE(double, DISTRIBUTED_LOAD_FACTOR)

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
class KratosIGAStructuralMechanicsApplication : public KratosApplication {
public:
	///@name Type Definitions
	///@{
	
	/// Pointer definition of KratosIGAStructuralMechanicsApplication
	KRATOS_CLASS_POINTER_DEFINITION(KratosIGAStructuralMechanicsApplication);

	///@}
	///@name Life Cycle
	///@{

	/// Default constructor.
	KratosIGAStructuralMechanicsApplication();

	/// Destructor.
	virtual ~KratosIGAStructuralMechanicsApplication(){}
	
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
		return "KratosIGAStructuralMechanicsApplication";
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
	// Meshless Elements
	const MeshlessBaseElement  mMeshlessElement;
	const MeshlessMembraneElement  mMeshlessMembraneElement;
	const MeshlessLaplaceElement  mMeshlessLaplaceElement;
	const MeshlessShellElement  mMeshlessShellElement;

	//// Outdated Conditions
	//const LoadCondition mLoadCondition;
	//const SupportCondition mSupportCondition;
	//const ContinuityConditionLagrange mContinuityConditionLagrange;
	//const ContinuityConditionPenalty mContinuityConditionPenalty;

	// Meshless Conditions
	const MeshlessSupportRotationCondition mMeshlessSupportRotationCondition;
	const MeshlessLoadCondition mMeshlessLoadCondition;
	const MeshlessLagrangeCouplingCondition mMeshlessLagrangeCouplingCondition;
	const MeshlessPenaltyCouplingRotationCondition mMeshlessPenaltyCouplingRotationCondition;
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
	KratosIGAStructuralMechanicsApplication& operator=(KratosIGAStructuralMechanicsApplication const& rOther);

	/// Copy constructor.
	KratosIGAStructuralMechanicsApplication(KratosIGAStructuralMechanicsApplication const& rOther);


	///@}

}; // Class KratosIGAStructuralMechanicsApplication

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

}  // namespace Kratos.

#endif // KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined
