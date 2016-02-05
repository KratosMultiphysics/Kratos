//--------------------------------------------------------------------
/*    |  /           |
      ' /   __| _` | __|  _ \   __|
      . \  |   (   | |   (   |\__ \ 
     _|\_\_|  \__,_|\__|\___/ ____/  
                   KRATOS  __|   _ \  |   |  _ \    
                         \__ \  (   | |   | | , )          
                         ____/ \___/ ___|_| ___/ MECHANICS 
			   
     License:		  SolidMechanicsApplication/license.txt
     Main authors:        Josep Maria Carbonell i Puigbo
                          ..                                        */
//--------------------------------------------------------------------
//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"

#include "containers/flags.h"

//conditions
#include "custom_conditions/point_torque_3D_condition.hpp"
#include "custom_conditions/point_load_2D_condition.hpp"

#include "custom_conditions/point_load_axisym_2D_condition.hpp"
#include "custom_conditions/point_load_3D_condition.hpp"

#include "custom_conditions/line_load_2D_condition.hpp"
#include "custom_conditions/line_load_axisym_2D_condition.hpp"
#include "custom_conditions/line_load_3D_condition.hpp"

#include "custom_conditions/surface_load_3D_condition.hpp"

//elements
// #include "custom_elements/small_displacement_beam_element_3D2N.hpp"
// #include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/small_displacement_element.hpp"
#include "custom_elements/axisym_small_displacement_element.hpp"

#include "custom_elements/total_lagrangian_element.hpp"
#include "custom_elements/updated_lagrangian_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_element.hpp"

#include "custom_elements/updated_lagrangian_U_P_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_P_element.hpp"
// #include "custom_elements/membrane_element.hpp"
// 
// #include "custom_elements/shell_thick_element_3D4N.hpp"
// #include "custom_elements/shell_thin_element_3D3N.hpp"

//flow rules
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"

//yield criteria
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"

//hardening laws
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"

//constitutive laws
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_plane_stress_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_simo_ju_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

// //cross sections
// #include "custom_utilities/shell_cross_section.hpp"

namespace Kratos
{
///@name Type Definitions
///@{
typedef array_1d<double,3> Vector3;
typedef array_1d<double,6> Vector6;
///@}

///@name Kratos Globals
///@{


//Define Variables

//for explicit schemes
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

//solution
KRATOS_DEFINE_VARIABLE(bool, COMPUTE_DYNAMIC_TANGENT )
KRATOS_DEFINE_VARIABLE(int, WRITE_ID )
KRATOS_DEFINE_VARIABLE(double, PREVIOUS_DELTA_TIME )
KRATOS_DEFINE_VARIABLE(double, NEWMARK_BETA )
KRATOS_DEFINE_VARIABLE(double, NEWMARK_GAMMA )
KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_ALPHA )
KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_BETA )

// //geometrical
// KRATOS_DEFINE_VARIABLE( double, AREA )
// KRATOS_DEFINE_VARIABLE( double, IX )
// KRATOS_DEFINE_VARIABLE( double, IY )
// KRATOS_DEFINE_VARIABLE( double, IZ )
// KRATOS_DEFINE_VARIABLE( double, CROSS_AREA )
// KRATOS_DEFINE_VARIABLE( double, MEAN_RADIUS )
// KRATOS_DEFINE_VARIABLE( int,    SECTION_SIDES )
// KRATOS_DEFINE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )

//constitutive law
KRATOS_DEFINE_VARIABLE(std::string, CONSTITUTIVE_LAW_NAME )
KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
KRATOS_DEFINE_VARIABLE(Matrix, CONSTITUTIVE_MATRIX )
KRATOS_DEFINE_VARIABLE(Matrix, DEFORMATION_GRADIENT )
KRATOS_DEFINE_VARIABLE(double, DETERMINANT_F )
KRATOS_DEFINE_VARIABLE(bool,   IMPLEX  )

// //cross section
// KRATOS_DEFINE_VARIABLE( ShellCrossSection::Pointer, SHELL_CROSS_SECTION )
// KRATOS_DEFINE_VARIABLE( int,          SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
// KRATOS_DEFINE_VARIABLE( double,	SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

//condition nodal load variables
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_LINE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_SURFACE_LOAD )

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( POINT_TORQUE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_TORQUE )

// //shell generalized variables
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_STRAIN )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_STRAIN_GLOBAL )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_CURVATURE )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_CURVATURE_GLOBAL )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_FORCE )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_FORCE_GLOBAL )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_MOMENT )
// KRATOS_DEFINE_VARIABLE( Matrix, SHELL_MOMENT_GLOBAL )

// //material orientation
KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DX )
KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DY )
KRATOS_DEFINE_VARIABLE( Vector3, MATERIAL_ORIENTATION_DZ )

// //othotropic/anisotropic constants
KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_X )
KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_Y )
KRATOS_DEFINE_VARIABLE( double, YOUNG_MODULUS_Z )
KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_XY )
KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_YZ )
KRATOS_DEFINE_VARIABLE( double, SHEAR_MODULUS_XZ )
KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_XY )
KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_YZ )
KRATOS_DEFINE_VARIABLE( double, POISSON_RATIO_XZ )

//material : hyperelastic_plastic
KRATOS_DEFINE_VARIABLE(double, NORM_ISOCHORIC_STRESS )
KRATOS_DEFINE_VARIABLE(double, PLASTIC_STRAIN )
KRATOS_DEFINE_VARIABLE(double, DELTA_PLASTIC_STRAIN )
KRATOS_DEFINE_VARIABLE(double, ISOTROPIC_HARDENING_MODULUS  )
KRATOS_DEFINE_VARIABLE(double, KINEMATIC_HARDENING_MODULUS  )
KRATOS_DEFINE_VARIABLE(double, HARDENING_EXPONENT )
KRATOS_DEFINE_VARIABLE(double, REFERENCE_HARDENING_MODULUS )
KRATOS_DEFINE_VARIABLE(double, INFINITY_HARDENING_MODULUS )

//material : isotropic damage
KRATOS_DEFINE_VARIABLE(double, DAMAGE_VARIABLE )
KRATOS_DEFINE_VARIABLE(double, DAMAGE_THRESHOLD )
KRATOS_DEFINE_VARIABLE(double, STRENGTH_RATIO )
KRATOS_DEFINE_VARIABLE(double, FRACTURE_ENERGY )

//thermal
KRATOS_DEFINE_VARIABLE(double, THERMAL_EXPANSION_COEFFICIENT )  
KRATOS_DEFINE_VARIABLE(double, REFERENCE_TEMPERATURE )
KRATOS_DEFINE_VARIABLE(double, PLASTIC_DISSIPATION )
KRATOS_DEFINE_VARIABLE(double, DELTA_PLASTIC_DISSIPATION )

//element
KRATOS_DEFINE_VARIABLE(Vector, RESIDUAL_VECTOR )
KRATOS_DEFINE_VARIABLE(Vector, EXTERNAL_FORCES_VECTOR )
KRATOS_DEFINE_VARIABLE(Vector, INTERNAL_FORCES_VECTOR )
KRATOS_DEFINE_VARIABLE(Vector, CONTACT_FORCES_VECTOR )

KRATOS_DEFINE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE(Vector, PK2_STRESS_VECTOR )

KRATOS_DEFINE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR )
KRATOS_DEFINE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
KRATOS_DEFINE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR )

KRATOS_DEFINE_VARIABLE(Matrix, MATERIAL_STIFFNESS_MATRIX )
KRATOS_DEFINE_VARIABLE(Matrix, GEOMETRIC_STIFFNESS_MATRIX )

KRATOS_DEFINE_VARIABLE(double, VON_MISES_STRESS )

//KRATOS_DEFINE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR );
//KRATOS_DEFINE_VARIABLE(Matrix, PK2_STRESS_TENSOR );
//KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR );

//nodal dofs
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION )
KRATOS_DEFINE_VARIABLE(double, PRESSURE_REACTION )


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
class KratosSolidMechanicsApplication : public KratosApplication
{
public:


    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosSolidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosSolidMechanicsApplication);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosSolidMechanicsApplication();

    /// Destructor.
    virtual ~KratosSolidMechanicsApplication() {}


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
    virtual std::string Info() const
    {
        return "KratosSolidMechanicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH( "in KratosSolidMechanicsApplication" )
        KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
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

//     //beams
// 
//     const SmallDisplacementBeamElement3D2N   mSmallDisplacementBeamElement3D2N;
// 
// 
//     //shells
// 
//     const IsotropicShellElement  mIsotropicShellElement3D3N;
//     const ShellThickElement3D4N  mShellThickElement3D4N;
//     const ShellThickElement3D4N  mShellThickCorotationalElement3D4N;
//     const ShellThinElement3D3N   mShellThinElement3D3N;
//     const ShellThinElement3D3N   mShellThinCorotationalElement3D3N;
	
    //solid

    //small displacement
    const SmallDisplacementElement mSmallDisplacementElement2D3N;
    const SmallDisplacementElement mSmallDisplacementElement2D4N;
    const SmallDisplacementElement mSmallDisplacementElement2D6N;
    const SmallDisplacementElement mSmallDisplacementElement2D8N;
    const SmallDisplacementElement mSmallDisplacementElement2D9N;

    const SmallDisplacementElement mSmallDisplacementElement3D4N;
    const SmallDisplacementElement mSmallDisplacementElement3D6N;
    const SmallDisplacementElement mSmallDisplacementElement3D8N;
    const SmallDisplacementElement mSmallDisplacementElement3D10N;
    const SmallDisplacementElement mSmallDisplacementElement3D15N;
    const SmallDisplacementElement mSmallDisplacementElement3D20N;
    const SmallDisplacementElement mSmallDisplacementElement3D27N;

    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D3N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D4N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D6N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D8N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D9N;

    //large displacement
    const LargeDisplacementElement     mLargeDisplacementElement;
    const LargeDisplacementUPElement mLargeDisplacementUPElement;

    //total lagrangian
    const TotalLagrangianElement mTotalLagrangianElement2D3N;
    const TotalLagrangianElement mTotalLagrangianElement2D4N;
    const TotalLagrangianElement mTotalLagrangianElement2D6N;
    const TotalLagrangianElement mTotalLagrangianElement2D8N;
    const TotalLagrangianElement mTotalLagrangianElement2D9N;

    const TotalLagrangianElement mTotalLagrangianElement3D4N;
    const TotalLagrangianElement mTotalLagrangianElement3D6N;
    const TotalLagrangianElement mTotalLagrangianElement3D8N;
    const TotalLagrangianElement mTotalLagrangianElement3D10N;
    const TotalLagrangianElement mTotalLagrangianElement3D15N;
    const TotalLagrangianElement mTotalLagrangianElement3D20N;
    const TotalLagrangianElement mTotalLagrangianElement3D27N;

 
    //updated lagrangian
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D3N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D4N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D6N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D8N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D9N;

    const UpdatedLagrangianElement mUpdatedLagrangianElement3D4N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D6N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D8N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D10N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D15N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D20N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D27N;

    const AxisymUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D3N;
    const AxisymUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D4N;
    const AxisymUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D6N;
    const AxisymUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D8N;
    const AxisymUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D9N;

    const UpdatedLagrangianUPElement             mUpdatedLagrangianUPElement2D3N;
    const AxisymUpdatedLagrangianUPElement mAxisymUpdatedLagrangianUPElement2D3N;
    const UpdatedLagrangianUPElement             mUpdatedLagrangianUPElement3D4N;
	
//     const MembraneElement mMembraneElement3D3N;
    
    //conditions
    const ForceLoadCondition                  mForceLoadCondition;

    const PointLoad2DCondition              mPointLoadCondition2D1N;
    const PointLoadAxisym2DCondition  mAxisymPointLoadCondition2D1N;
    const PointLoad3DCondition              mPointLoadCondition3D1N;
    const PointTorque3DCondition          mPointTorqueCondition3D1N;

    const LineLoad2DCondition              mLineLoadCondition2D2N;
    const LineLoad2DCondition              mLineLoadCondition2D3N;
    const LineLoadAxisym2DCondition  mAxisymLineLoadCondition2D2N;
    const LineLoadAxisym2DCondition  mAxisymLineLoadCondition2D3N;
    const LineLoad3DCondition              mLineLoadCondition3D2N;
    const LineLoad3DCondition              mLineLoadCondition3D3N;

    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D3N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D4N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D6N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D8N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D9N;


    //constitutive laws
    
    //Hyperelastic laws
    const HyperElastic3DLaw                       mHyperElastic3DLaw;
    const HyperElasticPlaneStrain2DLaw            mHyperElasticPlaneStrain2DLaw;
    const HyperElasticAxisym2DLaw                 mHyperElasticAxisym2DLaw;

    //Hyperelastic laws U-P
    const HyperElasticUP3DLaw                     mHyperElasticUP3DLaw;
    const HyperElasticUPPlaneStrain2DLaw          mHyperElasticUPPlaneStrain2DLaw;
    const HyperElasticUPAxisym2DLaw               mHyperElasticUPAxisym2DLaw;

    //Linear Elastic laws
    const LinearElastic3DLaw                      mLinearElastic3DLaw;
    const LinearElasticPlaneStrain2DLaw           mLinearElasticPlaneStrain2DLaw;
    const LinearElasticPlaneStress2DLaw           mLinearElasticPlaneStress2DLaw;
    const LinearElasticAxisym2DLaw                mLinearElasticAxisym2DLaw;

    //Hyperelastic Plastic laws
    const HyperElasticPlastic3DLaw                mHyperElasticPlastic3DLaw;
    const HyperElasticPlasticPlaneStrain2DLaw     mHyperElasticPlasticPlaneStrain2DLaw;
    const HyperElasticPlasticAxisym2DLaw          mHyperElasticPlasticAxisym2DLaw;    

    //Hyperelastic Plastic laws U-P
    const HyperElasticPlasticUP3DLaw              mHyperElasticPlasticUP3DLaw;
    const HyperElasticPlasticUPPlaneStrain2DLaw   mHyperElasticPlasticUPPlaneStrain2DLaw;
    const HyperElasticPlasticUPAxisym2DLaw        mHyperElasticPlasticUPAxisym2DLaw;    
 
    //Hyperelastic Plastic J2 specilization laws 
    const HyperElasticPlasticJ23DLaw              mHyperElasticPlasticJ23DLaw;
    const HyperElasticPlasticJ2PlaneStrain2DLaw   mHyperElasticPlasticJ2PlaneStrain2DLaw;
    const HyperElasticPlasticJ2Axisym2DLaw        mHyperElasticPlasticJ2Axisym2DLaw;

    //Hyperelastic Plastic J2 specilization laws U-P
    const HyperElasticPlasticUPJ23DLaw            mHyperElasticPlasticUPJ23DLaw;
    const HyperElasticPlasticUPJ2PlaneStrain2DLaw mHyperElasticPlasticUPJ2PlaneStrain2DLaw;
    const HyperElasticPlasticUPJ2Axisym2DLaw      mHyperElasticPlasticUPJ2Axisym2DLaw;
    
    //Linear Elastic Plastic Laws
    const LinearElasticPlastic3DLaw               mLinearElasticPlastic3DLaw;
    const LinearElasticPlasticPlaneStrain2DLaw    mLinearElasticPlasticPlaneStrain2DLaw;
    const LinearElasticPlasticPlaneStress2DLaw    mLinearElasticPlasticPlaneStress2DLaw;
    
    //Isotropic Damage Laws
    const IsotropicDamageSimoJu3DLaw              mIsotropicDamageSimoJu3DLaw;
    const IsotropicDamageSimoJuPlaneStrain2DLaw   mIsotropicDamageSimoJuPlaneStrain2DLaw;
    const IsotropicDamageSimoJuPlaneStress2DLaw   mIsotropicDamageSimoJuPlaneStress2DLaw;

    //Flow Rules
    const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticFlowRule;
    const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticFlowRule;
    const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
    
    //Yield Criteria
    const MisesHuberYieldCriterion                mMisesHuberYieldCriterion;
    const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
    
    //Hardening Laws
    const NonLinearIsotropicKinematicHardeningLaw mNonLinearIsotropicKinematicHardeningLaw;
    const LinearIsotropicKinematicHardeningLaw    mLinearIsotropicKinematicHardeningLaw;
    const ExponentialDamageHardeningLaw           mExponentialDamageHardeningLaw;

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
    KratosSolidMechanicsApplication& operator=(KratosSolidMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosSolidMechanicsApplication(KratosSolidMechanicsApplication const& rOther);


    ///@}

}; // Class KratosSolidMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED  defined 


