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
//   Date:                $Date:                December 2016 $
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

#include "containers/flags.h"

//outfitted python laws
/* #include "custom_python/python_outfitted_constitutive_law.hpp" */

/* //general constitutive laws */
/* //${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/elastic_3D_law.hpp */

/* //isotropic linear elastic laws */

/* //#include "custom_laws/linear_elastic_laws/linear_elastic_3D_law.hpp" */
/* #include "custom_laws/linear_elastic_laws/linear_elastic_plane_strain_2D_law.hpp" */
/* #include "custom_laws/linear_elastic_laws/linear_elastic_plane_stress_2D_law.hpp" */
/* #include "custom_laws/linear_elastic_laws/linear_elastic_axisymmetric_2D_law.hpp" */

/* //orthotropic linear elastic laws */
/* //#include "custom_laws/linear_elastic_laws/linear_elastic_orthotropic_3D_law.hpp" */

/* //isotropic hyperelastic laws */
/* //#include "custom_laws/linear_elastic_laws/hyperelastic_3D_law.hpp" */
/* #include "custom_laws/hyperelastic_laws/hyperelastic_plane_strain_2D_law.hpp" */
/* #include "custom_laws/hyperelastic_laws/hyperelastic_axisymmetric_2D_law.hpp" */

/* #include "custom_laws/hyperelastic_laws/hyperelastic_U_P_plane_strain_2D_law.hpp" */
/* #include "custom_laws/hyperelastic_laws/hyperelastic_U_P_axisymmetric_2D_law.hpp" */

/* //plasticity models */

/* //isotropic linear elastic plasticity laws */
/* //#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_3D_law.hpp" */
/* //#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_stress_2D_law.hpp" */
/* //#include "custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_axysimmetric_2D_law.hpp" */

/* //isotropic hyperelastic plasticity laws */
/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_axisymmetric_2D_law.hpp" */

/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_3D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_axisymmetric_2D_law.hpp" */

/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_3D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_J2_axisymmetric_2D_law.hpp" */

/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_3D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/hyperelastic_plastic_laws/incompressible_hyperelastic_plastic_J2_axisymmetric_2D_law.hpp" */

/* //isotropic linear elastic damage laws */
/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_3D_law.hpp" */
/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_stress_2D_law.hpp" */

/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_3D_law.hpp" */
/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_strain_2D_law.hpp" */
/* //#include "custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_stress_2D_law.hpp" */

/* //general constitutive models */

/* //hyperelastic models */
/* #include "custom_models/elasticity_models/hyperelastic_models/hyperelastic_model.hpp" */
/* #include "custom_models/elasticity_models/hyperelastic_models/neo_hookean_model.hpp" */

/* //hardening laws */
/* //#include "custom_models/plasticity_models/hardening_laws/hardening_law.hpp" */
/* #include "custom_models/plasticity_models/hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp" */
/* #include "custom_models/plasticity_models/hardening_laws/linear_isotropic_kinematic_hardening_law.hpp" */
/* #include "custom_models/plasticity_models/hardening_laws/exponential_damage_hardening_law.hpp" */
/* #include "custom_models/plasticity_models/hardening_laws/modified_exponential_damage_hardening_law.hpp" */

/* //yield criteria */
/* // #include "custom_models/plasticity_models/yield_criteria/yield_criterion.hpp" */
/* #include "custom_models/plasticity_models/yield_criteria/mises_huber_yield_criterion.hpp" */
/* #include "custom_models/plasticity_models/yield_criteria/simo_ju_yield_criterion.hpp" */
/* #include "custom_models/plasticity_models/yield_criteria/modified_mises_yield_criterion.hpp" */


#include "constitutive_models_application_variables.h"

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
  class KratosConstitutiveModelsApplication : public KratosApplication {
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
    virtual ~KratosConstitutiveModelsApplication(){}


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
      return "KratosConstitutiveModelsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {
      rOStream << Info();
      PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {
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

    /* //outfitted python laws */
    /* const PythonOutfittedConstitutiveLaw          mPythonOutfittedConstitutiveLaw; */
    
    /* //elasticity models */
    
    /* //isotropic linear elastic laws */
    /* const LinearElastic3DLaw                      mLinearElastic3DLaw; */
    /* const LinearElasticPlaneStrain2DLaw           mLinearElasticPlaneStrain2DLaw; */
    /* const LinearElasticPlaneStress2DLaw           mLinearElasticPlaneStress2DLaw; */
    /* const LinearElasticAxisymmetric2DLaw          mLinearElasticAxisymmetric2DLaw; */

    /* //orthotropic linear elastic laws */
    /* //const LinearElasticOrthotropic3DLaw           mLinearElasticOrthotropic3DLaw; */

    /* //isotropic hyperelastic laws */
    /* const HyperElastic3DLaw                       mHyperElastic3DLaw; */
    /* const HyperElasticPlaneStrain2DLaw            mHyperElasticPlaneStrain2DLaw; */
    /* const HyperElasticAxisymmetric2DLaw           mHyperElasticAxisymmetric2DLaw; */

    /* const HyperElasticUP3DLaw                     mHyperElasticUP3DLaw; */
    /* const HyperElasticUPPlaneStrain2DLaw          mHyperElasticUPPlaneStrain2DLaw; */
    /* const HyperElasticUPAxisymmetric2DLaw         mHyperElasticUPAxisymmetric2DLaw; */


    //plasticity models

    //isotropic linear elastic plasticity laws
    //const LinearElasticPlastic3DLaw               mLinearElasticPlastic3DLaw;
    //const LinearElasticPlasticPlaneStrain2DLaw    mLinearElasticPlasticPlaneStrain2DLaw;
    //const LinearElasticPlasticPlaneStress2DLaw    mLinearElasticPlasticPlaneStress2DLaw;
    //const LinearElasticPlasticPlaneStress2DLaw    mLinearElasticPlasticAxisymmetric2DLaw;
    
    //isotropic hyperelastic plasticity laws
    //const HyperElasticPlastic3DLaw                mHyperElasticPlastic3DLaw;
    //const HyperElasticPlasticPlaneStrain2DLaw     mHyperElasticPlasticPlaneStrain2DLaw;
    //const HyperElasticPlasticAxisymmetric2DLaw    mHyperElasticPlasticAxisymmetric2DLaw;    

    //const IncompressibleHyperElasticPlastic3DLaw              mIncompressibleHyperElasticPlastic3DLaw;
    //const IncompressibleHyperElasticPlasticPlaneStrain2DLaw   mIncompressibleHyperElasticPlasticPlaneStrain2DLaw;
    //const IncompressibleHyperElasticPlasticAxisymmetric2DLaw  mIncompressibleHyperElasticPlasticAxisymmetric2DLaw;  
 
    //const HyperElasticPlasticJ23DLaw              mHyperElasticPlasticJ23DLaw;
    //const HyperElasticPlasticJ2PlaneStrain2DLaw   mHyperElasticPlasticJ2PlaneStrain2DLaw;
    //const HyperElasticPlasticJ2Axisym2DLaw        mHyperElasticPlasticJ2Axisym2DLaw;

    //const IncompressibleHyperElasticPlasticUPJ23DLaw            mIncompressibleHyperElasticPlasticJ23DLaw;
    //const IncompressibleHyperElasticPlasticUPJ2PlaneStrain2DLaw mIncompressibleHyperElasticPlasticJ2PlaneStrain2DLaw;
    //const IncompressibleHyperElasticPlasticUPJ2Axisym2DLaw      mIncompressibleHyperElasticPlasticJ2Axisym2DLaw;
        

    //isotropic linear elastic damage laws
    //const IsotropicDamageSimoJu3DLaw              mIsotropicDamageSimoJu3DLaw;
    //const IsotropicDamageSimoJuPlaneStrain2DLaw   mIsotropicDamageSimoJuPlaneStrain2DLaw;
    //const IsotropicDamageSimoJuPlaneStress2DLaw   mIsotropicDamageSimoJuPlaneStress2DLaw;

    //const IsotropicDamageModifiedMises3DLaw            mIsotropicDamageModifiedMises3DLaw;
    //const IsotropicDamageModifiedMisesPlaneStrain2DLaw mIsotropicDamageModifiedMisesPlaneStrain2DLaw;
    //const IsotropicDamageModifiedMisesPlaneStress2DLaw mIsotropicDamageModifiedMisesPlaneStress2DLaw;

    /* //hyperelastic models */
    /* const HyperElasticModel                        mHyperElasticModel; */
    /* const NeoHookeanModel                          mNeoHookeanModel; */
    
    /* //hardening laws */
    /* const HardeningLaw                             mHardeningLaw; */
    /* const NonLinearIsotropicKinematicHardeningLaw  mNonLinearIsotropicKinematicHardeningLaw; */
    /* const LinearIsotropicKinematicHardeningLaw     mLinearIsotropicKinematicHardeningLaw; */
    /* const ExponentialDamageHardeningLaw            mExponentialDamageHardeningLaw; */
    /* const ModifiedExponentialDamageHardeningLaw    mModifiedExponentialDamageHardeningLaw; */
   
    /* //yield criteria */
    /* const MisesHuberYieldCriterion                 mMisesHuberYieldCriterion; */
    /* const SimoJuYieldCriterion                     mSimoJuYieldCriterion; */
    /* const ModifiedMisesYieldCriterion              mModifiedMisesYieldCriterion; */
    
    //plasticity models
    //const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticModel;
    //const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticModel;
    //const IsotropicDamageFlowRule                 mIsotropicDamageModel;
    
    

    
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
