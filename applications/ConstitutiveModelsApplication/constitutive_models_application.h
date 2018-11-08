//--------------------------------------------------------------------
//    |  /           |                                               .
//    ' /   __| _` | __|  _ \   __|                                  .
//    . \  |   (   | |   (   |\__ \                                  .
//   _|\_\_|  \__,_|\__|\___/ ____/                                  .
//                        __  __      _           _      _           .
//           CONSTITUTIVE|  \/  |__ _| |_ ___ _ _(_)__ _| |          .
//                       | |\/| / _` |  _/ -_) '_| / _` | |          .
//                       |_|  |_\__,_|\__\___|_| |_\__,_|_|MODELS    .
//			                                             .
//   License:(BSD)	  ConstitutiveModelsApplication/license.txt  .
//   Main authors:        Josep Maria Carbonell                      .
//                        ..                                         .
//--------------------------------------------------------------------
//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_CONSTITUTIVE_MODELS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_MODELS_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

//#include "containers/flags.h"

//outfitted python laws
//#include "custom_python/python_outfitted_constitutive_law.hpp"

//general constitutive laws

//small strain laws
#include "custom_laws/small_strain_laws/small_strain_orthotropic_3D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_strain_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_plane_stress_2D_law.hpp"
#include "custom_laws/small_strain_laws/small_strain_axisymmetric_2D_law.hpp"

//large strain laws
#include "custom_laws/large_strain_laws/large_strain_plane_strain_2D_law.hpp"
#include "custom_laws/large_strain_laws/large_strain_axisymmetric_2D_law.hpp"

//strain rate laws
#include "custom_laws/strain_rate_laws/strain_rate_plane_strain_2D_law.hpp"
#include "custom_laws/strain_rate_laws/newtonian_plane_strain_2D_law.hpp"

//specialized large strain laws

//elasticity models
#include "custom_models/elasticity_models/linear_elastic_model.hpp"
#include "custom_models/elasticity_models/saint_venant_kirchhoff_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/neo_hookean_J_1_squared_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/isochoric_neo_hookean_lnJ_squared_model.hpp"
#include "custom_models/elasticity_models/incompressible_neo_hookean_model.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
#include "custom_models/elasticity_models/tamagnini_model.hpp"
#include "custom_models/elasticity_models/ogden_model.hpp"
//#include "custom_models/elasticity_models/isochoric_ogden_model.hpp"
#include "custom_models/elasticity_models/incompressible_hypo_elastic_model.hpp"

//plasticity models
#include "custom_models/plasticity_models/von_mises_linear_elastic_plasticity_model.hpp"
#include "custom_models/plasticity_models/von_mises_neo_hookean_plasticity_model.hpp"
#include "custom_models/plasticity_models/simo_J2_plasticity_model.hpp"
#include "custom_models/plasticity_models/simo_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/johnson_cook_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/baker_johnson_cook_J2_thermo_plasticity_model.hpp"
#include "custom_models/plasticity_models/cam_clay_model.hpp"
#include "custom_models/plasticity_models/gens_nova_model.hpp"
#include "custom_models/plasticity_models/simo_ju_exponential_damage_model.hpp"
//#include "custom_models/plasticity_models/simo_ju_modified_exponential_damage_model.hpp"

//yield criteria
#include "custom_models/plasticity_models/yield_surfaces/mises_huber_thermal_yield_surface.hpp"
#include "custom_models/plasticity_models/yield_surfaces/simo_ju_yield_surface.hpp"
#include "custom_models/plasticity_models/yield_surfaces/modified_mises_yield_surface.hpp"
#include "custom_models/plasticity_models/yield_surfaces/modified_cam_clay_yield_surface.hpp"
#include "custom_models/plasticity_models/yield_surfaces/gens_nova_yield_surface.hpp"

//hardening rules
#include "custom_models/plasticity_models/hardening_rules/simo_linear_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/simo_exponential_thermal_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/johnson_cook_thermal_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/baker_johnson_cook_thermal_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/exponential_damage_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/modified_exponential_damage_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/cam_clay_hardening_rule.hpp"
#include "custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.hpp"


//#include "constitutive_models_application_variables.h"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) KratosConstitutiveModelsApplication : public KratosApplication
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosConstitutiveModelsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConstitutiveModelsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConstitutiveModelsApplication();

    /// Destructor.
    ~KratosConstitutiveModelsApplication() override{}


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
    std::string Info() const override{
      return "KratosConstitutiveModelsApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const  override{
      rOStream << Info();
      PrintData(rOStream);
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const  override{
      KRATOS_WATCH("in KratosConstitutiveModelsApplication");
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


    ///@}
    ///@name Member Variables
    ///@{

    //outfitted python laws
    //const PythonOutfittedConstitutiveLaw           mPythonOutfittedConstitutiveLaw;

    //general constitutive laws

    //small strain laws
    const SmallStrain3DLaw                         mSmallStrain3DLaw;
    const SmallStrainOrthotropic3DLaw              mSmallStrainOrthotropic3DLaw;
    const SmallStrainPlaneStrain2DLaw              mSmallStrainPlaneStrain2DLaw;
    const SmallStrainPlaneStress2DLaw              mSmallStrainPlaneStress2DLaw;
    const SmallStrainAxisymmetric2DLaw             mSmallStrainAxisymmetric2DLaw;

    //large strain laws
    const LargeStrain3DLaw                         mLargeStrain3DLaw;
    const LargeStrainPlaneStrain2DLaw              mLargeStrainPlaneStrain2DLaw;
    const LargeStrainAxisymmetric2DLaw             mLargeStrainAxisymmetric2DLaw;

    //strain rate laws
    const StrainRate3DLaw                          mStrainRate3DLaw;
    const StrainRatePlaneStrain2DLaw               mStrainRatePlaneStrain2DLaw;
    const Newtonian3DLaw                           mNewtonian3DLaw;
    const NewtonianPlaneStrain2DLaw                mNewtonianPlaneStrain2DLaw;

    //general constitutive models

    //elasticity models
    const LinearElasticModel                       mLinearElasticModel;
    const SaintVenantKirchhoffModel                mSaintVenantKirchhoffModel;
    const NeoHookeanModel                          mNeoHookeanModel;
    const NeoHookeanLnJSquaredModel                mNeoHookeanLnJSquaredModel;
    const NeoHookeanJ_1SquaredModel                mNeoHookeanJ_1SquaredModel;
    const IsochoricNeoHookeanModel                 mIsochoricNeoHookeanModel;
    const IsochoricNeoHookeanLnJSquaredModel       mIsochoricNeoHookeanLnJSquaredModel;
    const IncompressibleNeoHookeanModel            mIncompressibleNeoHookeanModel;
    const BorjaModel                               mBorjaModel;
    const TamagniniModel                           mTamagniniModel;
    const OgdenModel                               mOgdenModel;
    const OgdenModel                               mIsochoricOgdenModel;
    const HypoElasticModel                         mHypoElasticModel;
    const IsochoricHypoElasticModel                mIsochoricHypoElasticModel;
    const IncompressibleHypoElasticModel           mIncompressibleHypoElasticModel;

    //plasticity models
    const VonMisesLinearElasticPlasticityModel     mVonMisesLinearElasticPlasticityModel;
    const VonMisesNeoHookeanPlasticityModel        mVonMisesNeoHookeanPlasticityModel;
    const SimoJ2PlasticityModel                    mSimoJ2PlasticityModel;
    const SimoJ2ThermoPlasticityModel              mSimoJ2ThermoPlasticityModel;
    const JohnsonCookJ2ThermoPlasticityModel       mJohnsonCookJ2ThermoPlasticityModel;
    const BakerJohnsonCookJ2ThermoPlasticityModel  mBakerJohnsonCookJ2ThermoPlasticityModel;
    const CamClayModel                             mCamClayModel;
    const GensNovaModel                            mGensNovaModel;
    const SimoJuExponentialDamageModel             mSimoJuExponentialDamageModel;
    const SimoJuExponentialDamageModel             mSimoJuModifiedExponentialDamageModel;

    //yield criteria
    const MisesHuberYieldSurface<HardeningRule>         mMisesHuberYieldSurface;
    const MisesHuberThermalYieldSurface<HardeningRule>  mMisesHuberThermalYieldSurface;
    const SimoJuYieldSurface<HardeningRule>             mSimoJuYieldSurface;
    const ModifiedMisesYieldSurface<HardeningRule>      mModifiedMisesYieldSurface;
    const ModifiedCamClayYieldSurface<HardeningRule>    mModifiedCamClayYieldSurface;
    const GensNovaYieldSurface<HardeningRule>           mGensNovaYieldSurface;
    
    //hardening rules
    const SimoExponentialHardeningRule              mSimoExponentialHardeningRule;
    const SimoLinearHardeningRule                   mSimoLinearHardeningRule;
    const SimoExponentialThermalHardeningRule       mSimoExponentialThermalHardeningRule;
    const JohnsonCookThermalHardeningRule           mJohnsonCookThermalHardeningRule;
    const BakerJohnsonCookThermalHardeningRule      mBakerJohnsonCookThermalHardeningRule;
    const ExponentialDamageHardeningRule            mExponentialDamageHardeningRule;
    const ModifiedExponentialDamageHardeningRule    mModifiedExponentialDamageHardeningRule;
    const CamClayHardeningRule                      mCamClayHardeningRule;
    const GensNovaHardeningRule                     mGensNovaHardeningRule;
      
       
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
    KratosConstitutiveModelsApplication& operator=(KratosConstitutiveModelsApplication const& rOther);

    /// Copy constructor.
    KratosConstitutiveModelsApplication(KratosConstitutiveModelsApplication const& rOther);


    ///@}

  }; // Class KratosConstitutiveModelsApplication

  ///@}


  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KRATOS_CONSTITUTIVE_MODELS_APPLICATION_H_INCLUDED  defined
