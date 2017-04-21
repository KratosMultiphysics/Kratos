//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes


// External includes


// Project includes
#include "constitutive_models_application.h"


namespace Kratos {

  KratosConstitutiveModelsApplication::KratosConstitutiveModelsApplication() {}

  void KratosConstitutiveModelsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "             __  __      _           _      _         "<< std::endl;
    std::cout << "CONSTITUTIVE|  \\/  |__ _| |_ ___ _ _(_)__ _| |        "<< std::endl;
    std::cout << "            | |\\/| / _` |  _/ -_) '_| / _` | |        "<< std::endl;
    std::cout << "            |_|  |_\\__,_|\\__\\___|_| |_\\__,_|_|MODELS" << std::endl;
    std::cout << "Initializing KratosConstitutiveModelsApplication... " << std::endl;

    //Register Variables (variables created in constitutive_models_application_variables.cpp) 

    
      KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NAME )
      KRATOS_REGISTER_VARIABLE( IMPLEX )
	
      //hyperelasticity
      //KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS )
      //KRATOS_REGISTER_VARIABLE( POISSON_RATIO )
      KRATOS_REGISTER_VARIABLE( LAME_MU )
      KRATOS_REGISTER_VARIABLE( LAME_LAMBDA )
      KRATOS_REGISTER_VARIABLE( HYPERELASTIC_MODEL_PARAMETERS )

      //plasticity
      KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS )
      KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN )
      KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_STRAIN )
      KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
      KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_MODULUS )
      KRATOS_REGISTER_VARIABLE( HARDENING_EXPONENT )
      KRATOS_REGISTER_VARIABLE( REFERENCE_HARDENING_MODULUS )
      KRATOS_REGISTER_VARIABLE( INFINITY_HARDENING_MODULUS )

      //isotropic damage
      KRATOS_REGISTER_VARIABLE( DAMAGE_VARIABLE )
      KRATOS_REGISTER_VARIABLE( DAMAGE_THRESHOLD )
      KRATOS_REGISTER_VARIABLE( STRENGTH_RATIO )
      KRATOS_REGISTER_VARIABLE( FRACTURE_ENERGY )
      KRATOS_REGISTER_VARIABLE( RESIDUAL_STRENGTH )
      KRATOS_REGISTER_VARIABLE( SOFTENING_SLOPE )
      
      //thermal
      KRATOS_REGISTER_VARIABLE( REFERENCE_TEMPERATURE )
      KRATOS_REGISTER_VARIABLE( PLASTIC_DISSIPATION )
      KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_DISSIPATION )

      //othotropic/anisotropic constants
      KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
      KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
      KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
      KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
      KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
      KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
      KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
      KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
      KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DX )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DY )
      KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DZ )

      
      //KRATOS_REGISTER_VARIABLE( THERMAL_EXPANSION_COEFFICIENT )
      //KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS )
      //KRATOS_REGISTER_VARIABLE( BULK_MODULUS )


      //Register Constitutive Laws

      //outfitted python laws
      Serializer::Register( "PythonOutfittedConstitutiveLaw", mPythonOutfittedConstitutiveLaw );

      //general constitutive laws
      
      //elasticity laws
      
      //isotropic linear elastic laws
      Serializer::Register( "LinearElastic3DLaw", mLinearElastic3DLaw );
      Serializer::Register( "LinearElasticPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw );
      Serializer::Register( "LinearElasticPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw );
      Serializer::Register( "LinearElasticAxisymmetric2DLaw", mLinearElasticAxisymmetric2DLaw );

      //orthotropic linear elastic laws
      Serializer::Register( "LinearElasticOrthotropic3DLaw", mLinearElasticOrthotropic3DLaw );
      
      //isotropic hyperelastic laws
      Serializer::Register( "HyperElastic3DLaw", mHyperElastic3DLaw );
      Serializer::Register( "HyperElasticPlaneStrain2DLaw", mHyperElasticPlaneStrain2DLaw );
      Serializer::Register( "HyperElasticAxisymmetric2DLaw", mHyperElasticAxisymmetric2DLaw );

      Serializer::Register( "HyperElasticUP3DLaw", mHyperElasticUP3DLaw );
      Serializer::Register( "HyperElasticUPPlaneStrain2DLaw", mHyperElasticUPPlaneStrain2DLaw );
      Serializer::Register( "HyperElasticUPAxisymmetric2DLaw", mHyperElasticUPAxisymmetric2DLaw );

            
      //plasticity laws

      //isotropic linear elastic plasticity laws
      //Serializer::Register( "LinearElasticPlastic3DLaw", mLinearElasticPlastic3DLaw );
      //Serializer::Register( "LinearElasticPlasticPlaneStrain2DLaw", mLinearElasticPlasticPlaneStrain2DLaw );
      //Serializer::Register( "LinearElasticPlasticPlaneStress2DLaw", mLinearElasticPlasticPlaneStress2DLaw );

      //isotropic hyperelastic plasticity laws
      Serializer::Register( "HyperElasticPlastic3DLaw", mHyperElasticPlastic3DLaw );
      Serializer::Register( "HyperElasticPlasticPlaneStrain2DLaw", mHyperElasticPlasticPlaneStrain2DLaw );
      Serializer::Register( "HyperElasticPlasticAxisymmetric2DLaw", mHyperElasticPlasticAxisymmetric2DLaw );

      Serializer::Register( "HyperElasticPlasticUP3DLaw", mHyperElasticPlasticUP3DLaw );
      Serializer::Register( "HyperElasticPlasticUPPlaneStrain2DLaw", mHyperElasticPlasticUPPlaneStrain2DLaw );
      Serializer::Register( "HyperElasticPlasticUPAxisymmetric2DLaw", mHyperElasticPlasticUPAxisymmetric2DLaw );


      //isotropic linear elastic damage laws
      // Serializer::Register( "IsotropicDamageSimoJu3DLaw", mIsotropicDamageSimoJu3DLaw );
      // Serializer::Register( "IsotropicDamageSimoJuPlaneStrain2DLaw", mIsotropicDamageSimoJuPlaneStrain2DLaw );
      // Serializer::Register( "IsotropicDamageSimoJuPlaneStress2DLaw", mIsotropicDamageSimoJuPlaneStress2DLaw );

      // Serializer::Register( "IsotropicDamageModifiedMises3DLaw", mIsotropicDamageModifiedMises3DLaw );
      // Serializer::Register( "IsotropicDamageModifiedMisesPlaneStrain2DLaw", mIsotropicDamageModifiedMisesPlaneStrain2DLaw );
      // Serializer::Register( "IsotropicDamageModifiedMisesPlaneStress2DLaw", mIsotropicDamageModifiedMisesPlaneStress2DLaw );

      //elasticity models
      Serializer::Register( "LinearElasticModel", mLinearElasticModel );
      
      //hyperelastic models
      Serializer::Register( "SaintVenantKirchhoffModel", mSaintVenantKirchhoffModel );
      Serializer::Register( "NeoHookeanModel", mNeoHookeanModel );
      Serializer::Register( "CompressibleNeoHookeanModel", mCompressibleNeoHookeanModel );
      Serializer::Register( "IsochoricNeoHookeanModel", mIsochoricNeoHookeanModel );
      Serializer::Register( "IncompressibleNeoHookeanModel", mIncompressibleNeoHookeanModel );
      
      //plasticity model
      Serializer::Register( "NonLinearAssociativePlasticityModel", mNonLinearAssociativePlasticityModel );
      Serializer::Register( "VonMisesPlasticityModel", mVonMisesPlasticityModel );
          
      //yield criteria
      Serializer::Register( "MisesHuberYieldCriterion", mMisesHuberYieldCriterion );
      Serializer::Register( "SimoJuYieldCriterion", mSimoJuYieldCriterion );
      Serializer::Register( "ModifiedMisesYieldCriterion", mModifiedMisesYieldCriterion );
    
      //hardening laws
      Serializer::Register( "NonLinearIsotropicKinematicHardeningLaw", mNonLinearIsotropicKinematicHardeningLaw );
      Serializer::Register( "LinearIsotropicKinematicHardeningLaw", mLinearIsotropicKinematicHardeningLaw );
      Serializer::Register( "ExponentialDamageHardeningLaw", mExponentialDamageHardeningLaw );
      Serializer::Register( "ModifiedExponentialDamageHardeningLaw", mModifiedExponentialDamageHardeningLaw );


  }
}  // namespace Kratos.
