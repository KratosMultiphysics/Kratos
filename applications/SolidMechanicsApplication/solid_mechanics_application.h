//------------------------------------------------------------------
//           ___      _ _    _                                     .
//   KRATOS / __| ___| (_)__| |                                    .
//          \__ \/ _ \ | / _` |                                    .
//          |___/\___/_|_\__,_| MECHANICS                          .
//			                                           .
//   License:(BSD)	  SolidMechanicsApplication/license.txt    .
//   Main authors:        Josep Maria Carbonell                    .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED)
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

//elements

//solid elements
#include "custom_elements/solid_elements/linear_solid_element.hpp"

#include "custom_elements/solid_elements/small_displacement_element.hpp"
#include "custom_elements/solid_elements/small_displacement_bbar_element.hpp"
#include "custom_elements/solid_elements/axisymmetric_small_displacement_element.hpp"

#include "custom_elements/solid_elements/total_lagrangian_element.hpp"
#include "custom_elements/solid_elements/updated_lagrangian_element.hpp"
#include "custom_elements/solid_elements/axisymmetric_updated_lagrangian_element.hpp"

#include "custom_elements/solid_elements/updated_lagrangian_U_P_element.hpp"
#include "custom_elements/solid_elements/axisymmetric_updated_lagrangian_U_P_element.hpp"

#include "custom_elements/solid_elements/updated_lagrangian_V_element.hpp"
#include "custom_elements/solid_elements/updated_lagrangian_segregated_V_P_element.hpp"

//beam elements
#include "custom_elements/beam_elements/beam_element.hpp"
#include "custom_elements/beam_elements/small_displacement_beam_element.hpp"
#include "custom_elements/beam_elements/small_displacement_beam_element_3D2N.hpp"
#include "custom_elements/beam_elements/large_displacement_beam_element.hpp"
#include "custom_elements/beam_elements/large_displacement_beam_emc_element.hpp"
#include "custom_elements/beam_elements/large_displacement_beam_semc_element.hpp"
#include "custom_elements/beam_elements/geometrically_exact_rod_element.hpp"

//shell elements
#include "custom_elements/shell_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_elements/shell_thin_element_3D3N.hpp"

//thermal elements
#include "custom_elements/thermal_elements/thermal_element.hpp"
#include "custom_elements/thermal_elements/axisymmetric_thermal_element.hpp"

//conditions
#include "custom_conditions/load_conditions/axisymmetric_point_load_condition.hpp"
#include "custom_conditions/load_conditions/axisymmetric_line_load_condition.hpp"
#include "custom_conditions/load_conditions/surface_load_condition.hpp"

#include "custom_conditions/moment_conditions/point_moment_condition.hpp"
#include "custom_conditions/moment_conditions/line_moment_condition.hpp"
#include "custom_conditions/moment_conditions/surface_moment_condition.hpp"

#include "custom_conditions/elastic_conditions/axisymmetric_point_elastic_condition.hpp"
#include "custom_conditions/elastic_conditions/axisymmetric_line_elastic_condition.hpp"
#include "custom_conditions/elastic_conditions/surface_elastic_condition.hpp"

#include "custom_conditions/thermal_conditions/line_heat_flux_condition.hpp"

//flow rules
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/non_linear_rate_dependent_plastic_flow_rule.hpp"

//yield criteria
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/mises_huber_thermal_yield_criterion.hpp"

//hardening laws
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_thermal_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/johnson_cook_thermal_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/baker_johnson_cook_thermal_hardening_law.hpp"

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
#include "custom_constitutive/linear_elastic_orthotropic_3D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_J2_axisym_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_simo_ju_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.hpp"

#include "custom_constitutive/isotropic_damage_modified_mises_3D_law.hpp"
#include "custom_constitutive/isotropic_damage_modified_mises_plane_strain_2D_law.hpp"
#include "custom_constitutive/isotropic_damage_modified_mises_plane_stress_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_thermal_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_johnson_cook_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_baker_johnson_cook_plane_strain_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_johnson_cook_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_johnson_cook_axisym_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_thermal_U_P_baker_johnson_cook_plane_strain_2D_law.hpp"

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
 class KRATOS_API(SOLID_MECHANICS_APPLICATION) KratosSolidMechanicsApplication : public KratosApplication
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
   ~KratosSolidMechanicsApplication() override {}


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
     return "KratosSolidMechanicsApplication";
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
   const LinearSolidElement mLinearSolidElement2D3N;
   const LinearSolidElement mLinearSolidElement2D4N;
   const LinearSolidElement mLinearSolidElement3D4N;
   const LinearSolidElement mLinearSolidElement3D8N;


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

   //B-bar
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D3N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D4N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D6N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D8N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement2D9N;

   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D4N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D6N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D8N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D10N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D15N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D20N;
   const SmallDisplacementBbarElement mSmallDisplacementBbarElement3D27N;

   const AxisymmetricSmallDisplacementElement mAxisymSmallDisplacementElement2D3N;
   const AxisymmetricSmallDisplacementElement mAxisymSmallDisplacementElement2D4N;
   const AxisymmetricSmallDisplacementElement mAxisymSmallDisplacementElement2D6N;
   const AxisymmetricSmallDisplacementElement mAxisymSmallDisplacementElement2D8N;
   const AxisymmetricSmallDisplacementElement mAxisymSmallDisplacementElement2D9N;

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

   const AxisymmetricUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D3N;
   const AxisymmetricUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D4N;
   const AxisymmetricUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D6N;
   const AxisymmetricUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D8N;
   const AxisymmetricUpdatedLagrangianElement mAxisymUpdatedLagrangianElement2D9N;

   //velocity based elements
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement2D3N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement2D4N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement2D6N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement2D8N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement2D9N;

   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D4N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D6N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D8N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D10N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D15N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D20N;
   const UpdatedLagrangianVElement mUpdatedLagrangianVElement3D27N;

   //segregated VP elements
   const UpdatedLagrangianSegregatedVPElement mUpdatedLagrangianSegregatedVPElement2D3N;
   const UpdatedLagrangianSegregatedVPElement mUpdatedLagrangianSegregatedVPElement3D4N;

   //mixed elements UP
   const UpdatedLagrangianUPElement         mUpdatedLagrangianUPElement2D3N;
   const AxisymmetricUpdatedLagrangianUPElement mAxisymUpdatedLagrangianUPElement2D3N;
   const UpdatedLagrangianUPElement         mUpdatedLagrangianUPElement3D4N;

   //beams
   const SmallDisplacementBeamElement       mSmallDisplacementBeamElement3D2N;
   const LargeDisplacementBeamElement       mLargeDisplacementBeamElement3D2N;
   const LargeDisplacementBeamElement       mLargeDisplacementBeamElement3D3N;
   const LargeDisplacementBeamEMCElement    mLargeDisplacementBeamEMCElement3D2N;
   const LargeDisplacementBeamEMCElement    mLargeDisplacementBeamEMCElement3D3N;
   const LargeDisplacementBeamSEMCElement   mLargeDisplacementBeamSEMCElement3D2N;
   const GeometricallyExactRodElement       mGeometricallyExactRodElement3D2N;
   const LargeDisplacementBeamElement       mLargeDisplacementBeamElement2D2N;

   //shells
   const ShellThickElement3D4N              mShellThickElement3D4N;
   const ShellThickElement3D4N  mShellThickCorotationalElement3D4N;
   const ShellThinElement3D3N                mShellThinElement3D3N;
   const ShellThinElement3D3N    mShellThinCorotationalElement3D3N;

   //thermal
   const ThermalElement             mThermalElement2D3N;
   const ThermalElement             mThermalElement3D4N;

   const AxisymmetricThermalElement mAxisymThermalElement2D3N;


   //conditions
   const PointLoadCondition                    mPointLoadCondition3D1N;
   const PointLoadCondition                    mPointLoadCondition2D1N;
   const AxisymmetricPointLoadCondition  mAxisymPointLoadCondition2D1N;

   const LineLoadCondition                      mLineLoadCondition3D2N;
   const LineLoadCondition                      mLineLoadCondition3D3N;
   const LineLoadCondition                      mLineLoadCondition2D2N;
   const LineLoadCondition                      mLineLoadCondition2D3N;
   const AxisymmetricLineLoadCondition    mAxisymLineLoadCondition2D2N;
   const AxisymmetricLineLoadCondition    mAxisymLineLoadCondition2D3N;

   const SurfaceLoadCondition                mSurfaceLoadCondition3D3N;
   const SurfaceLoadCondition                mSurfaceLoadCondition3D4N;
   const SurfaceLoadCondition                mSurfaceLoadCondition3D6N;
   const SurfaceLoadCondition                mSurfaceLoadCondition3D8N;
   const SurfaceLoadCondition                mSurfaceLoadCondition3D9N;

   const PointMomentCondition                mPointMomentCondition3D1N;
   const PointMomentCondition                mPointMomentCondition2D1N;

   const LineMomentCondition                  mLineMomentCondition3D2N;
   const LineMomentCondition                  mLineMomentCondition3D3N;
   const LineMomentCondition                  mLineMomentCondition2D2N;
   const LineMomentCondition                  mLineMomentCondition2D3N;

   const SurfaceMomentCondition            mSurfaceMomentCondition3D3N;
   const SurfaceMomentCondition            mSurfaceMomentCondition3D4N;
   const SurfaceMomentCondition            mSurfaceMomentCondition3D6N;
   const SurfaceMomentCondition            mSurfaceMomentCondition3D8N;
   const SurfaceMomentCondition            mSurfaceMomentCondition3D9N;

   const PointElasticCondition                    mPointElasticCondition3D1N;
   const PointElasticCondition                    mPointElasticCondition2D1N;
   const AxisymmetricPointElasticCondition  mAxisymPointElasticCondition2D1N;

   const LineElasticCondition                      mLineElasticCondition3D2N;
   const LineElasticCondition                      mLineElasticCondition3D3N;
   const LineElasticCondition                      mLineElasticCondition2D2N;
   const LineElasticCondition                      mLineElasticCondition2D3N;
   const AxisymmetricLineElasticCondition    mAxisymLineElasticCondition2D2N;
   const AxisymmetricLineElasticCondition    mAxisymLineElasticCondition2D3N;

   const SurfaceElasticCondition                mSurfaceElasticCondition3D3N;
   const SurfaceElasticCondition                mSurfaceElasticCondition3D4N;
   const SurfaceElasticCondition                mSurfaceElasticCondition3D6N;
   const SurfaceElasticCondition                mSurfaceElasticCondition3D8N;
   const SurfaceElasticCondition                mSurfaceElasticCondition3D9N;

   const LineHeatFluxCondition                    mLineHeatFluxCondition2D2N;


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

   //Hyperelastic Plastic J2 specilization laws
   const HyperElasticPlasticJ23DLaw              mHyperElasticPlasticJ23DLaw;
   const HyperElasticPlasticJ2PlaneStrain2DLaw   mHyperElasticPlasticJ2PlaneStrain2DLaw;
   const HyperElasticPlasticJ2Axisym2DLaw        mHyperElasticPlasticJ2Axisym2DLaw;

   //Hyperelastic Plastic J2 specilization laws U-P
   const HyperElasticPlasticUPJ23DLaw            mHyperElasticPlasticUPJ23DLaw;
   const HyperElasticPlasticUPJ2PlaneStrain2DLaw mHyperElasticPlasticUPJ2PlaneStrain2DLaw;
   const HyperElasticPlasticUPJ2Axisym2DLaw      mHyperElasticPlasticUPJ2Axisym2DLaw;

   //Isotropic Damage Laws
   const IsotropicDamageSimoJu3DLaw              mIsotropicDamageSimoJu3DLaw;
   const IsotropicDamageSimoJuPlaneStrain2DLaw   mIsotropicDamageSimoJuPlaneStrain2DLaw;
   const IsotropicDamageSimoJuPlaneStress2DLaw   mIsotropicDamageSimoJuPlaneStress2DLaw;

   const IsotropicDamageModifiedMises3DLaw            mIsotropicDamageModifiedMises3DLaw;
   const IsotropicDamageModifiedMisesPlaneStrain2DLaw mIsotropicDamageModifiedMisesPlaneStrain2DLaw;
   const IsotropicDamageModifiedMisesPlaneStress2DLaw mIsotropicDamageModifiedMisesPlaneStress2DLaw;

   //Thermal Laws
   const HyperElasticPlasticThermalJ2PlaneStrain2DLaw mHyperElasticPlasticThermalJ2PlaneStrain2DLaw;
   const HyperElasticPlasticThermalJohnsonCookPlaneStrain2DLaw mHyperElasticPlasticThermalJohnsonCookPlaneStrain2DLaw;
   const HyperElasticPlasticThermalBakerJohnsonCookPlaneStrain2DLaw mHyperElasticPlasticThermalBakerJohnsonCookPlaneStrain2DLaw;

   const HyperElasticPlasticThermalUPJ23DLaw mHyperElasticPlasticThermalUPJ23DLaw;
   const HyperElasticPlasticThermalUPJ2PlaneStrain2DLaw mHyperElasticPlasticThermalUPJ2PlaneStrain2DLaw;
   const HyperElasticPlasticThermalUPJ2Axisym2DLaw mHyperElasticPlasticThermalUPJ2Axisym2DLaw;
   const HyperElasticPlasticThermalUPJohnsonCookPlaneStrain2DLaw mHyperElasticPlasticThermalUPJohnsonCookPlaneStrain2DLaw;
   const HyperElasticPlasticThermalUPJohnsonCookAxisym2DLaw mHyperElasticPlasticThermalUPJohnsonCookAxisym2DLaw;
   const HyperElasticPlasticThermalUPBakerJohnsonCookPlaneStrain2DLaw mHyperElasticPlasticThermalUPBakerJohnsonCookPlaneStrain2DLaw;


   //Flow Rules
   const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticFlowRule;
   const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticFlowRule;
   const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
   const NonLinearRateDependentPlasticFlowRule   mNonLinearRateDependentPlasticFlowRule;

   //Yield Criteria
   const MisesHuberYieldCriterion                mMisesHuberYieldCriterion;
   const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
   const ModifiedMisesYieldCriterion             mModifiedMisesYieldCriterion;
   const MisesHuberThermalYieldCriterion         mMisesHuberThermalYieldCriterion;

   //Hardening Laws
   const NonLinearIsotropicKinematicHardeningLaw mNonLinearIsotropicKinematicHardeningLaw;
   const LinearIsotropicKinematicHardeningLaw    mLinearIsotropicKinematicHardeningLaw;
   const ExponentialDamageHardeningLaw           mExponentialDamageHardeningLaw;
   const ModifiedExponentialDamageHardeningLaw   mModifiedExponentialDamageHardeningLaw;

   const NonLinearIsotropicKinematicThermalHardeningLaw mNonLinearIsotropicKinematicThermalHardeningLaw;
   const JohnsonCookThermalHardeningLaw                 mJohnsonCookThermalHardeningLaw;
   const BakerJohnsonCookThermalHardeningLaw            mBakerJohnsonCookThermalHardeningLaw;


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
