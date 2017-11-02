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

    std::cout << "            __  __      _           _      _         "<< std::endl;
    std::cout << "     KRATOS|  \\/  |__ _| |_ ___ _ _(_)__ _| |        "<< std::endl;
    std::cout << "           | |\\/| / _` |  _/ -_) '_| / _` | |        "<< std::endl;
    std::cout << "           |_|  |_\\__,_|\\__\\___|_| |_\\__,_|_|MODELS" << std::endl;
    std::cout << "Initializing KratosConstitutiveModelsApplication... " << std::endl;

    //Register Variables (variables created in constitutive_models_application_variables.cpp) 

    //specific constitutive models variables must be REGISTERED here
    
    //Register Constitutive Laws

    //outfitted python laws
    Serializer::Register( "PythonOutfittedConstitutiveLaw", mPythonOutfittedConstitutiveLaw );
    
    //general constitutive laws
    
    //elasticity laws
    
    //small strain laws
    Serializer::Register( "SmallStrain3DLaw", mSmallStrain3DLaw );
    Serializer::Register( "SmallStrainOrthotropic3DLaw", mSmallStrainOrthotropic3DLaw );
    Serializer::Register( "SmallStrainPlaneStrain2DLaw", mSmallStrainPlaneStrain2DLaw );
    Serializer::Register( "SmallStrainPlaneStress2DLaw", mSmallStrainPlaneStress2DLaw );
    Serializer::Register( "SmallStrainAxisymmetric2DLaw", mSmallStrainAxisymmetric2DLaw );
    
    //large strain laws
    Serializer::Register( "LargeStrain3DLaw", mLargeStrain3DLaw );
    Serializer::Register( "LargeStrainPlaneStrain2DLaw", mLargeStrainPlaneStrain2DLaw );
    Serializer::Register( "LargeStrainAxisymmetric2DLaw", mLargeStrainAxisymmetric2DLaw );
    
    //general constitutive models
    
    //elasticity models
    Serializer::Register( "LinearElasticModel", mLinearElasticModel );
    Serializer::Register( "SaintVenantKirchhoffModel", mSaintVenantKirchhoffModel );
    Serializer::Register( "NeoHookeanModel", mNeoHookeanModel );
    Serializer::Register( "NeoHookeanLnJSquaredModel", mNeoHookeanLnJSquaredModel );
    Serializer::Register( "NeoHookeanJ_1SquaredModel", mNeoHookeanJ_1SquaredModel );
    Serializer::Register( "IsochoricNeoHookeanModel", mIsochoricNeoHookeanModel );
    Serializer::Register( "IsochoricNeoHookeanLnJSquaredModel", mIsochoricNeoHookeanLnJSquaredModel );
    Serializer::Register( "IncompressibleNeoHookeanModel", mIncompressibleNeoHookeanModel );
    Serializer::Register( "BorjaModel", mBorjaModel );
    
    //plasticity models
    Serializer::Register( "VonMisesLinearElasticPlasticityModel", mVonMisesLinearElasticPlasticityModel );
    Serializer::Register( "VonMisesNeoHookeanPlasticityModel", mVonMisesNeoHookeanPlasticityModel );
    Serializer::Register( "SimoJ2PlasticityModel", mSimoJ2PlasticityModel );
    Serializer::Register( "CamClayModel", mCamClayModel );
    Serializer::Register( "SimoJ2ThermoPlasticityModel", mSimoJ2ThermoPlasticityModel );
	
    //yield criteria
    Serializer::Register( "MisesHuberYieldSurface", mMisesHuberYieldSurface );
    Serializer::Register( "MisesHuberThermalYieldSurface", mMisesHuberThermalYieldSurface );      
    Serializer::Register( "SimoJuYieldSurface", mSimoJuYieldSurface );
    Serializer::Register( "ModifiedMisesYieldSurface", mModifiedMisesYieldSurface );
    Serializer::Register( "ModifiedCamClaySurface", mModifiedCamClayYieldSurface );
    
    //hardening rules
    Serializer::Register( "SimoExponentialHardeningRule", mSimoExponentialHardeningRule );
    Serializer::Register( "SimoLinearHardeningRule", mSimoLinearHardeningRule );
    Serializer::Register( "SimoExponentialThermalHardeningRule", mSimoExponentialThermalHardeningRule );
    Serializer::Register( "JohnsonCookThermalHardeningRule", mJohnsonCookThermalHardeningRule );
    Serializer::Register( "BakerJohnsonCookThermalHardeningRule", mBakerJohnsonCookThermalHardeningRule );
    Serializer::Register( "ExponentialDamageHardeningRule", mExponentialDamageHardeningRule );
    Serializer::Register( "ModifiedExponentialDamageHardeningRule", mModifiedExponentialDamageHardeningRule );
    Serializer::Register( "CamClayHardeningRule", mCamClayHardeningRule );
    Serializer::Register( "SimoJuExponentialDamageModel", mSimoJuExponentialDamageModel );
    Serializer::Register( "SimoJuModifiedExponentialDamageModel", mSimoJuModifiedExponentialDamageModel );
      

  }
}  // namespace Kratos.
