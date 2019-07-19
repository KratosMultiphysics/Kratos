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

KratosConstitutiveModelsApplication::KratosConstitutiveModelsApplication()
    :KratosApplication("ConstitutiveModelsApplication")
{}

void KratosConstitutiveModelsApplication::Register() {
  // calling base class register to register Kratos components
  KratosApplication::Register();

  std::stringstream banner;

  banner << "            __  __      _           _      _          \n"
         << "    KRATOS |  \\/  |__ _| |_ ___ _ _(_)__ _| |         \n"
         << "           | |\\/| / _` |  _/ -_) '_| / _` | |         \n"
         << "           |_|  |_\\__,_|\\__\\___|_| |_\\__,_|_| MODELS\n"
         << "Initialize KratosConstitutiveModelsApplication... " << std::endl;

  // mpi initialization
  int mpi_is_initialized = 0;
  int rank = -1;

#ifdef KRATOS_MPI

  MPI_Initialized(&mpi_is_initialized);

  if (mpi_is_initialized)
  {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  }

#endif

  if (mpi_is_initialized)
  {
    if (rank == 0) KRATOS_INFO("") << banner.str();
  }
  else
  {
    KRATOS_INFO("") << banner.str();
  }

  //Register Variables (variables created in constitutive_models_application_variables.cpp)
  KRATOS_REGISTER_VARIABLE(TEMPERATURE_VARIABLE)
  KRATOS_REGISTER_VARIABLE(PRESSURE_VARIABLE)
  KRATOS_REGISTER_VARIABLE(PROPERTIES_LAYOUT)

   KRATOS_REGISTER_VARIABLE( DILATANCY_ANGLE )
   KRATOS_REGISTER_VARIABLE( RHOS )   
   KRATOS_REGISTER_VARIABLE( RHOT )   
   KRATOS_REGISTER_VARIABLE( KSIM )   
   KRATOS_REGISTER_VARIABLE( CHIS )   
   KRATOS_REGISTER_VARIABLE( CHIT )   
   KRATOS_REGISTER_VARIABLE( REFERENCE_PRESSURE )   

   KRATOS_REGISTER_VARIABLE( CASM_M )
   KRATOS_REGISTER_VARIABLE( ELASTIC_JACOBIAN )

   KRATOS_REGISTER_VARIABLE( SPACING_RATIO )
   KRATOS_REGISTER_VARIABLE( SHAPE_PARAMETER )

   KRATOS_REGISTER_VARIABLE( PS )   
   KRATOS_REGISTER_VARIABLE( PT )   
   KRATOS_REGISTER_VARIABLE( PM )  

   KRATOS_REGISTER_VARIABLE( COHESION )

   KRATOS_REGISTER_VARIABLE( PLASTIC_VOL_DEF )   
   KRATOS_REGISTER_VARIABLE( NONLOCAL_PLASTIC_VOL_DEF )   
   KRATOS_REGISTER_VARIABLE( PLASTIC_VOL_DEF_ABS )   
   KRATOS_REGISTER_VARIABLE( NONLOCAL_PLASTIC_VOL_DEF_ABS )   
   KRATOS_REGISTER_VARIABLE( PLASTIC_DEV_DEF )   
   KRATOS_REGISTER_VARIABLE( NONLOCAL_PLASTIC_DEV_DEF )   

   KRATOS_REGISTER_VARIABLE( MAX_IMPLEX_PLASTIC_MULTIPLIER )
   KRATOS_REGISTER_VARIABLE( NORMALIZE_FLOW_RULE )
  //specific constitutive models variables must be REGISTERED here

  //Register Constitutive Laws

  //outfitted python laws
  //Serializer::Register( "PythonOutfittedConstitutiveLaw", mPythonOutfittedConstitutiveLaw );

  //general constitutive laws

  //elasticity laws

  //small strain laws
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "SmallStrain3DLaw", mSmallStrain3DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "SmallStrainOrthotropic3DLaw", mSmallStrainOrthotropic3DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "SmallStrainPlaneStrain2DLaw", mSmallStrainPlaneStrain2DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "SmallStrainPlaneStress2DLaw", mSmallStrainPlaneStress2DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "SmallStrainAxisymmetric2DLaw", mSmallStrainAxisymmetric2DLaw );

  //large strain laws
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "LargeStrain3DLaw", mLargeStrain3DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "LargeStrainPlaneStrain2DLaw", mLargeStrainPlaneStrain2DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "LargeStrainAxisymmetric2DLaw", mLargeStrainAxisymmetric2DLaw );

  //strain rate laws
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "StrainRate3DLaw", mStrainRate3DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "StrainRatePlaneStrain2DLaw", mStrainRatePlaneStrain2DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "NewtonianFluid3DLaw", mNewtonianFluid3DLaw );
  KRATOS_REGISTER_CONSTITUTIVE_LAW( "NewtonianFluidPlaneStrain2DLaw", mNewtonianFluidPlaneStrain2DLaw );

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
  Serializer::Register( "TamagniniModel", mTamagniniModel );
  Serializer::Register( "OgdenModel", mOgdenModel );
  Serializer::Register( "IsochoricOgdenModel", mIsochoricOgdenModel );
  Serializer::Register( "HypoElasticModel", mHypoElasticModel );
  Serializer::Register( "IsochoricHypoElasticModel", mIsochoricHypoElasticModel );
  Serializer::Register( "IncompressibleHypoElasticModel", mIncompressibleHypoElasticModel );
  Serializer::Register( "HenckyLinearModel", mHenckyLinearModel );

  //plasticity models
  Serializer::Register( "VonMisesLinearElasticPlasticityModel", mVonMisesLinearElasticPlasticityModel );
  Serializer::Register( "VonMisesNeoHookeanPlasticityModel", mVonMisesNeoHookeanPlasticityModel );
  Serializer::Register( "SimoJ2PlasticityModel", mSimoJ2PlasticityModel );
  Serializer::Register( "CamClayModel", mCamClayModel );
  Serializer::Register( "NonlocalCamClayModel", mNonlocalCamClayModel );
  Serializer::Register( "GensNovaModel", mGensNovaModel );
  Serializer::Register( "V2GensNovaModel", mV2GensNovaModel );
  Serializer::Register( "NonlocalV2GensNovaModel", mNonlocalV2GensNovaModel );
  Serializer::Register( "SimoJ2ThermoPlasticityModel", mSimoJ2ThermoPlasticityModel );
  Serializer::Register( "MohrCoulombV1Model", mMohrCoulombV1Model );
  Serializer::Register( "MohrCoulombNonAssociativeModel", mMohrCoulombNonAssociativeModel );
  Serializer::Register( "TrescaModel", mTrescaModel );
  Serializer::Register( "CasmAssociatedSoilModel", mCasmAssociatedSoilModel );
  Serializer::Register( "CasmMCCSoilModel", mCasmMCCSoilModel );
  Serializer::Register( "CasmNubiaSoilModel", mCasmNubiaSoilModel );
  Serializer::Register( "NonlocalCasmMCCSoilModel", mNonlocalCasmMCCSoilModel );
  Serializer::Register( "NonlocalCasmNubiaSoilModel", mNonlocalCasmNubiaSoilModel );

  //yield criteria
  Serializer::Register( "MisesHuberYieldSurface", mMisesHuberYieldSurface );
  Serializer::Register( "MisesHuberThermalYieldSurface", mMisesHuberThermalYieldSurface );
  Serializer::Register( "SimoJuYieldSurface", mSimoJuYieldSurface );
  Serializer::Register( "ModifiedMisesYieldSurface", mModifiedMisesYieldSurface );
  Serializer::Register( "ModifiedCamClayYieldSurface", mModifiedCamClayYieldSurface );
  Serializer::Register( "MohrCoulombV1YieldSurface", mMohrCoulombV1YieldSurface );
  Serializer::Register( "MohrCoulombNonAssociativeYieldSurface", mMohrCoulombNonAssociativeYieldSurface );
  Serializer::Register( "TrescaYieldSurface", mTrescaYieldSurface );

  //hardening rules
  Serializer::Register( "SimoExponentialHardeningRule", mSimoExponentialHardeningRule );
  Serializer::Register( "SimoLinearHardeningRule", mSimoLinearHardeningRule );
  Serializer::Register( "SimoExponentialThermalHardeningRule", mSimoExponentialThermalHardeningRule );
  Serializer::Register( "JohnsonCookThermalHardeningRule", mJohnsonCookThermalHardeningRule );
  Serializer::Register( "BakerJohnsonCookThermalHardeningRule", mBakerJohnsonCookThermalHardeningRule );
  Serializer::Register( "ExponentialDamageHardeningRule", mExponentialDamageHardeningRule );
  Serializer::Register( "ModifiedExponentialDamageHardeningRule", mModifiedExponentialDamageHardeningRule );
  Serializer::Register( "CamClayHardeningRule", mCamClayHardeningRule );
  Serializer::Register( "GensNovaHardeningRule", mGensNovaHardeningRule);
  Serializer::Register( "MohrCoulombV1HardeningRule", mMohrCoulombV1HardeningRule );


  Serializer::Register( "SimoJuExponentialDamageModel", mSimoJuExponentialDamageModel );
  Serializer::Register( "SimoJuModifiedExponentialDamageModel", mSimoJuModifiedExponentialDamageModel );


}
}  // namespace Kratos.
