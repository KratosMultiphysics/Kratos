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
#include "structural_mechanics_application.h"
#include "structural_mechanics_application_variables.h"

/* ELEMENTS */
#include "custom_elements/meshless_base_element.h"
#include "custom_elements/meshless_base_surface_element.h"
#include "custom_elements/meshless_membrane_element.h"
#include "custom_elements/meshless_laplace_element.h"
#include "custom_elements/meshless_shell_element.h"
#include "custom_elements/meshless_shell_kl_element.h"

//#include "custom_conditions/ContinuityConditionLagrange.h"
//#include "custom_conditions/ContinuityConditionPenalty.h"
//#include "custom_conditions/LoadCondition.h"
//#include "custom_conditions/SupportCondition.h"

/*CONDITIONS*/
#include "custom_conditions/meshless_support_rotation_condition.h"
#include "custom_conditions/meshless_surface_support_condition.h"
#include "custom_conditions/meshless_load_condition.h"
#include "custom_conditions/meshless_lagrange_coupling_condition.h"
#include "custom_conditions/meshless_lagrange_coupling_condition_2.h"
#include "custom_conditions/meshless_penalty_coupling_rotation_condition.h"
#include "custom_conditions/meshless_penalty_coupling_crack_condition.h"

namespace Kratos {

///@name Kratos Globals
///@{

	//Variables definition
	KRATOS_DEFINE_VARIABLE(double, INTEGRATION_WEIGHT)
	KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES)
	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES)
	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)

	KRATOS_DEFINE_VARIABLE(Vector, SHAPE_FUNCTION_VALUES_SLAVE)
	KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)
	//KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER)
	//KRATOS_DEFINE_VARIABLE(Matrix, SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE)

	KRATOS_DEFINE_VARIABLE(Vector, TANGENTS)
	KRATOS_DEFINE_VARIABLE(Vector, TANGENTS_SLAVE)

	KRATOS_DEFINE_VARIABLE(double, PENALTY_FACTOR)

	KRATOS_DEFINE_VARIABLE(int, DISPLACEMENT_ROTATION_FIX)
	KRATOS_DEFINE_VARIABLE(int, LOAD_TYPE)
	//KRATOS_DEFINE_VARIABLE(std::string, CONDITION_TYPE_DEFINITION)
	KRATOS_DEFINE_VARIABLE(double, DISTRIBUTED_LOAD_FACTOR)

	// for damage constitutive law
	KRATOS_DEFINE_VARIABLE(Vector, GAP_INTERFACE)
	KRATOS_DEFINE_VARIABLE(double, CONVECTION_DEGRADATION)
	KRATOS_DEFINE_VARIABLE(int, EXPONENTIAL_DAMAGE)

	KRATOS_DEFINE_VARIABLE(double, DAMAGE_T_INTERFACE)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_C_INTERFACE)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_T)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_C)
	KRATOS_DEFINE_VARIABLE(double, FRACTURE_ENERGY_T)
	KRATOS_DEFINE_VARIABLE(double, FRACTURE_ENERGY_C)
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_T) /** @todo: to be removed*/
	KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS_C) /** @todo: to be removed*/
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRESS_T_0)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRESS_C_0)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRESS_C_P)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRESS_C_M) /** @todo: to be removed*/
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRESS_C_R)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRAIN_C_P)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_STRAIN_C_M) /** @todo: to be removed*/
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_COMPRESSIVE_LAW_C1)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_COMPRESSIVE_LAW_C2)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_COMPRESSIVE_LAW_C3)
	KRATOS_DEFINE_VARIABLE(double, BIAXIAL_COMPRESSION_MULTIPLIER)
	KRATOS_DEFINE_VARIABLE(double, SHEAR_COMPRESSION_REDUCTION)
	KRATOS_DEFINE_VARIABLE(double, DAMAGE_TENSILE_SURFACE_S1)
	KRATOS_DEFINE_VARIABLE(double, LUBLINER_SURFACE_PARAM_KC)
	KRATOS_DEFINE_VARIABLE(double, GENRANKINE_SURFACE_PARAM_A)
	KRATOS_DEFINE_VARIABLE(double, GENRANKINE_SURFACE_PARAM_B)
	KRATOS_DEFINE_VARIABLE(double, GENRANKINE_SURFACE_PARAM_C)
	KRATOS_DEFINE_VARIABLE(int, DAMAGE_SECANT_MATRIX)
	KRATOS_DEFINE_VARIABLE(int, DAMAGE_MODEL)
	KRATOS_DEFINE_VARIABLE(int, DAMAGE_TENSILE_MODEL)

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
	const MeshlessBaseElement     mMeshlessElement;
	const MeshlessMembraneElement mMeshlessMembraneElement;
	const MeshlessLaplaceElement  mMeshlessLaplaceElement;
	const MeshlessShellElement    mMeshlessShellElement;
	const MeshlessShellKLElement  mMeshlessShellKLElement;

	// Meshless Conditions
	const MeshlessSupportRotationCondition         mMeshlessSupportRotationCondition;
    const MeshlessSurfaceSupportCondition          mMeshlessSurfaceSupportCondition;
	const MeshlessLoadCondition                    mMeshlessLoadCondition;
	const MeshlessLagrangeCouplingCondition        mMeshlessLagrangeCouplingCondition;
    const MeshlessLagrangeCouplingCondition2       mMeshlessLagrangeCouplingCondition2;
	const MeshlessPenaltyCouplingRotationCondition mMeshlessPenaltyCouplingRotationCondition;
	const MeshlessPenaltyCouplingCrackCondition    mMeshlessPenaltyCouplingCrackCondition;
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
