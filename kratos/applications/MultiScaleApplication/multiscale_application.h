//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_MULTISCALE_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTISCALE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

#include "includes/kratos_application.h"

#include "includes/variables.h"
#include "includes/ublas_interface.h"

#include "custom_elements/small_displacement_interface_element.hpp"
#include "custom_elements/opt_triangle_element.hpp"
#include "custom_elements/eas_quad_element_v2.hpp"
#include "custom_elements/q4ristab_element.hpp"
#include "custom_elements/convdiff_interface_element.hpp"
#include "custom_elements/convdiff_element.hpp"

#include "custom_elements/ebst_element_2d3n.h"
#include "custom_elements/agq4_element.hpp"
#include "custom_conditions/periodic_condition_lm_2D2N.h"

namespace Kratos
{
///@name Kratos Globals
///@{

#define RVE_TEST_MOD_CONDITIONS

#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___ALWAYS_FINALIZE   0
#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_OR_ABORT 1
#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___ALWAYS_ABORT      2
#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___DO_NOTHING        3
#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___FINALIZE_ONLY     4
#define RVE_STRATEGY_FINALIZE_STEP_LEVEL___ABORT_ONLY        5

// for rve
KRATOS_DEFINE_VARIABLE( int, RVE_CONSTITUTIVE_LAW_FLAG )
KRATOS_DEFINE_VARIABLE( int, RVE_DAMAGE_SURFACE_FLAG )
KRATOS_DEFINE_VARIABLE(double, RVE_DAMAGE_SURFACE_LIMIT)
KRATOS_DEFINE_VARIABLE(std::string, RVE_CLAW_MAP_NAME)
//KRATOS_DEFINE_VARIABLE(TagStrainVectorMap, RVE_CLAW_MAP)

KRATOS_DEFINE_VARIABLE(Vector, FLUX_RVE)
KRATOS_DEFINE_VARIABLE(Vector, HEAT_FLUX_RVE)
KRATOS_DEFINE_VARIABLE(double, HOMOGENIZED_CONDUCTIVITY)
KRATOS_DEFINE_VARIABLE(double, HOMOGENIZED_CTE)
KRATOS_DEFINE_VARIABLE(Matrix, HOMOGENIZED_CONST_TENS)
KRATOS_DEFINE_VARIABLE(Matrix, INVERSE_HOMOGENIZED_CONST_TENS)
KRATOS_DEFINE_VARIABLE(double, RVE_NON_LINEAR_FLAG)
KRATOS_DEFINE_VARIABLE(Vector, ACTUAL_TAG)
KRATOS_DEFINE_VARIABLE(Matrix, RVE_GENERAL_STRESS_TENSOR)
KRATOS_DEFINE_VARIABLE(double, EQUIVALENT_DAMAGE)

KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_1111)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_1122)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_1133)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_2211)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_2222)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_2233)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_3311)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_3322)
KRATOS_DEFINE_VARIABLE(double, CONDUCTIVITY_3333)

KRATOS_DEFINE_VARIABLE(double, CLAW_LIMIT_DAMAGE)
KRATOS_DEFINE_VARIABLE(double, ASPERITIES_CONSTANT)
KRATOS_DEFINE_VARIABLE(double, BETA_THERMAL_COEFF)
KRATOS_DEFINE_VARIABLE(double, YOUNG_MAT1)
KRATOS_DEFINE_VARIABLE(double, CONDUCT_MAT1)
KRATOS_DEFINE_VARIABLE(double, POISSON_MAT1)
KRATOS_DEFINE_VARIABLE(double, YOUNG_MAT2)
KRATOS_DEFINE_VARIABLE(double, CONDUCT_MAT2)
KRATOS_DEFINE_VARIABLE(double, POISSON_MAT2)
KRATOS_DEFINE_VARIABLE(double, AIR_CONDUCTIVITY)

KRATOS_DEFINE_VARIABLE(double, RVE_FULL_TEMPERATURE)
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_REACTION)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( RVE_FULL_DISPLACEMENT )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( RVE_WPC_LAGRANGIAN_DOF )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( RVE_WPC_LAGRANGIAN_REACTION )
KRATOS_DEFINE_VARIABLE( double, RVE_WPR_LAGRANGIAN_DOF )
KRATOS_DEFINE_VARIABLE( double, RVE_WPR_LAGRANGIAN_REACTION )

// for lagrange multipliers
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( REACTION_DISPLACEMENT_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( REACTION_ROTATION_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_1 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_2 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_1 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_2 )
KRATOS_DEFINE_VARIABLE( double, DOUBLE_LAGRANGE_SCALE_FACTOR_ALPHA )
KRATOS_DEFINE_VARIABLE( double, DOUBLE_LAGRANGE_SCALE_FACTOR_BETA )

// for strategies
KRATOS_DEFINE_VARIABLE( int, STRATEGY_SOLUTION_STEP_SOLVED )
KRATOS_DEFINE_VARIABLE( int, STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL )
KRATOS_DEFINE_VARIABLE( double, CONSTITUTIVE_INTEGRATION_ERROR_CODE )
KRATOS_DEFINE_VARIABLE( int, ITERATION_CONVERGENCE_FLAG )
KRATOS_DEFINE_VARIABLE( double, SUGGESTED_TIME_STEP )

// for damage constitutive law
KRATOS_DEFINE_VARIABLE( Vector, GAP_INTERFACE )
KRATOS_DEFINE_VARIABLE( double, CONVECTION_DEGRADATION )
KRATOS_DEFINE_VARIABLE( int, EXPONENTIAL_DAMAGE )

KRATOS_DEFINE_VARIABLE( double, DAMAGE_T )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_C )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_T )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_C )
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_T ) /** @todo: to be removed*/
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_C ) /** @todo: to be removed*/
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRESS_T_0 )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRESS_C_0 )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRESS_C_P )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRESS_C_M ) /** @todo: to be removed*/
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRESS_C_R )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRAIN_C_P )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_STRAIN_C_M ) /** @todo: to be removed*/
KRATOS_DEFINE_VARIABLE( double, DAMAGE_COMPRESSIVE_LAW_C1 )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_COMPRESSIVE_LAW_C2 )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_COMPRESSIVE_LAW_C3 )
KRATOS_DEFINE_VARIABLE( double, BIAXIAL_COMPRESSION_MULTIPLIER )
KRATOS_DEFINE_VARIABLE( double, SHEAR_COMPRESSION_REDUCTION )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_TENSILE_SURFACE_S1 )
KRATOS_DEFINE_VARIABLE( double, LUBLINER_SURFACE_PARAM_KC )
KRATOS_DEFINE_VARIABLE( double, GENRANKINE_SURFACE_PARAM_A )
KRATOS_DEFINE_VARIABLE( double, GENRANKINE_SURFACE_PARAM_B )
KRATOS_DEFINE_VARIABLE( double, GENRANKINE_SURFACE_PARAM_C )
KRATOS_DEFINE_VARIABLE( int, DAMAGE_SECANT_MATRIX )
KRATOS_DEFINE_VARIABLE( int, DAMAGE_MODEL )
KRATOS_DEFINE_VARIABLE( int, DAMAGE_TENSILE_MODEL )

// for plots
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_2D_X )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_2D_Y )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_X )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_Y )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_Z )

// for interface constitutive law
KRATOS_DEFINE_VARIABLE( double, NORMAL_STIFFNESS )
KRATOS_DEFINE_VARIABLE( double, TANGENTIAL_STIFFNESS )
KRATOS_DEFINE_VARIABLE( double, NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_MODE_I )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_MODE_II )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_MODE_III )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_I )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_II )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_III )
KRATOS_DEFINE_VARIABLE( double, INITIAL_COHESION )
KRATOS_DEFINE_VARIABLE( double, INITIAL_FRICTION_ANGLE )
KRATOS_DEFINE_VARIABLE( double, RESIDUAL_FRICTION_ANGLE )
KRATOS_DEFINE_VARIABLE( double, INITIAL_DILATANCY_ANGLE )
KRATOS_DEFINE_VARIABLE( double, RESIDUAL_DILATANCY_ANGLE )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_TENSILE_LAW_S0 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_S0 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_SP )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_SR )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_EP )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_C1 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_C2 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_C3 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_COMPRESSIVE_LAW_C4 )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_PLASTIC_DAMAGE_FACTOR_T )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_PLASTIC_DAMAGE_FACTOR_C )
KRATOS_DEFINE_VARIABLE( double, INTERFACE_CAP_VALUE )
KRATOS_DEFINE_VARIABLE( Vector, INTERFACE_TRACTION )
KRATOS_DEFINE_VARIABLE( Vector, INTERFACE_DISPLACEMENT_JUMP )
KRATOS_DEFINE_VARIABLE( Vector, INTERFACE_PLASTIC_DISPLACEMENT_JUMP )
KRATOS_DEFINE_VARIABLE( double, YIELD_FUNCTION_VALUE )
KRATOS_DEFINE_VARIABLE( int, INTERFACE_REDUCED_INTEGRATION )

// for plastic constitutive law
KRATOS_DEFINE_VARIABLE( double, ISOTROPIC_HARDENING )
KRATOS_DEFINE_VARIABLE( double, KINEMATIC_HARDENING )
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_INFINITY )
KRATOS_DEFINE_VARIABLE( double, ISOTROPIC_HARDENING_EXPONENT )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_PLASTIC_STRAIN )
KRATOS_DEFINE_VARIABLE( Matrix, PLASTIC_STRAIN_TENSOR )

// for stabilized reduced integration
KRATOS_DEFINE_VARIABLE( double, RI_STABILIZATION )
KRATOS_DEFINE_VARIABLE( double, RI_STABILIZATION_RESIDUAL )

// for enhanced strain elements
KRATOS_DEFINE_VARIABLE( double, ENH_STRAIN_PARAM_1 )
KRATOS_DEFINE_VARIABLE( double, ENH_STRAIN_PARAM_2 )
KRATOS_DEFINE_VARIABLE( double, ENH_STRAIN_PARAM_3 )
KRATOS_DEFINE_VARIABLE( double, ENH_STRAIN_PARAM_4 )
KRATOS_DEFINE_VARIABLE( double, ENH_STRAIN_PARAM_5 )

// misc
KRATOS_DEFINE_VARIABLE( double, RANDOM_IMPERFECTION_FACTOR )
KRATOS_DEFINE_VARIABLE( Vector, DISCONTINUITY_DIRECTION )
KRATOS_DEFINE_VARIABLE( double, LAMBDA_OUTPUT )

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

/// MultiScale Application for KRATOS.
/**
 * This application is a Demo
 */
class KratosMultiScaleApplication : public KratosApplication
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMultiScaleApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultiScaleApplication);
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMultiScaleApplication();

    /// Destructor.
    virtual ~KratosMultiScaleApplication() {}

    ///@}

    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * Registers the structural application in the KRATOS kernel
     */
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
        return "KratosMultiScaleApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in KratosMultiScaleApplication application");
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

private:

    ///@name Member Variables
    ///@{

	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement2D4N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D6N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D8N;
	
	const OptTriangleElement mOptTriangleElement2D3N;
	const EASQuadElementV2 mEASQuadElementV22D4N;
	const Q4RIStabElement mQ4RIStabElement2D4N;

	const ConvDiffInterfaceElement mConvDiffInterfaceElement2D4N;
	const ConvDiffInterfaceElement mConvDiffInterfaceElement3D6N;
	const ConvDiffInterfaceElement mConvDiffInterfaceElement3D8N;

	const ConvDiffElement mConvDiffElement2D3N;
	const ConvDiffElement mConvDiffElement2D4N;
	const ConvDiffElement mConvDiffElement3D4N;
	const ConvDiffElement mConvDiffElement3D8N;
	
	const EBSTElement2D3N mEBSTElement2D3N;
	const AGQ4Element mAGQ4Element2D4N;

	const PeriodicConditionLM2D2N mPeriodicConditionLM2D2N;

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
    KratosMultiScaleApplication& operator=(KratosMultiScaleApplication const& rOther);

    /// Copy constructor.
    KratosMultiScaleApplication(KratosMultiScaleApplication const& rOther);

    ///@}

};

///@}

}
#endif // KRATOS_MULTISCALE_APPLICATION_H_INCLUDED  defined 
