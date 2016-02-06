//--------------------------------------------------------------------
//    |  /           |                                               .
//    ' /   __| _` | __|  _ \   __|                                  .
//    . \  |   (   | |   (   |\__ \                                  .
//   _|\_\_|  \__,_|\__|\___/ ____/                                  .
//                 KRATOS  __|   _ \  |   |  _ \                     .
//                       \__ \  (   | |   | | , )                    .      
//                       |___/ \___/ ___|_| ___/ MECHANICS           .            
//			                                             .
//   License:(BSD)	  SolidMechanicsApplication/license.txt      .
//   Main authors:        Josep Maria Carbonell                      .
//                        ..                                         .
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
#include "includes/constitutive_law.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"

#include "containers/flags.h"

//conditions
#include "custom_conditions/point_load_2D_condition.hpp"
#include "custom_conditions/point_load_axisym_2D_condition.hpp"
#include "custom_conditions/point_load_3D_condition.hpp"

#include "custom_conditions/line_load_2D_condition.hpp"
#include "custom_conditions/line_load_axisym_2D_condition.hpp"
#include "custom_conditions/line_load_3D_condition.hpp"

#include "custom_conditions/surface_load_3D_condition.hpp"

//elements
#include "custom_elements/small_displacement_element.hpp"
#include "custom_elements/axisym_small_displacement_element.hpp"

#include "custom_elements/total_lagrangian_element.hpp"
#include "custom_elements/updated_lagrangian_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_element.hpp"

#include "custom_elements/updated_lagrangian_U_P_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_P_element.hpp"

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

#include "solid_mechanics_application_variables.h"

namespace Kratos
{
///@name Type Definitions
///@{
typedef array_1d<double,3> Vector3;
typedef array_1d<double,6> Vector6;
///@}

///@name Kratos Globals
///@{

//Application variables definition:  (see solid_mechanics_application_variables.h)

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


   ///@}
   ///@name Member Variables
   ///@{


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
	
   //conditions
   const ForceLoadCondition                  mForceLoadCondition;

   const PointLoad2DCondition              mPointLoadCondition2D1N;
   const PointLoadAxisym2DCondition  mAxisymPointLoadCondition2D1N;
   const PointLoad3DCondition              mPointLoadCondition3D1N;

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


