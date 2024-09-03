// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "constitutive_laws_application.h"
#include "constitutive_laws_application_variables.h"


namespace Kratos {

KratosConstitutiveLawsApplication::KratosConstitutiveLawsApplication():
    KratosApplication("ConstitutiveLawsApplication")
    {}

void KratosConstitutiveLawsApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosConstitutiveLawsApplication..." << std::endl;

    // Damage and plasticity
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticOrthotropic2DLaw", mLinearElasticOrthotropic2DLaw);
    // Register hyper elastic laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenant3DLaw", mHyperElasticIsotropicKirchhoff3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenantPlaneStress2DLaw", mHyperElasticIsotropicKirchhoffPlaneStress2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenantPlaneStrain2DLaw", mHyperElasticIsotropicKirchhoffPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElastic3DLaw", mHyperElasticIsotropicNeoHookean3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticPlaneStrain2DLaw", mHyperElasticIsotropicNeoHookeanPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2PlasticityPlaneStrain2DLaw", mSmallStrainJ2PlasticityPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2Plasticity3DLaw", mSmallStrainJ2Plasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainPlasticity3DAthLaw", mSmallStrainPlasticity3DAth);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamagePlaneStrain2DLaw", mSmallStrainIsotropicDamagePlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DLaw", mSmallStrainIsotropicDamage3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageAthira3DLaw", mSmallStrainIsotropicDamageAthira3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMazarsLaw", mSmallStrainIsotropicDamage3DMazars);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageImplex3DLaw", mSmallStrainIsotropicDamageImplex3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageTractionOnly3DLaw", mSmallStrainIsotropicDamageTractionOnly3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageTractionOnlyImplex3DLaw", mSmallStrainIsotropicDamageTractionOnlyImplex3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("TrussPlasticityConstitutiveLaw", mTrussPlasticityConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicOgden1D", mHyperElasticIsotropicOgden1D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicHenky1D", mHyperElasticIsotropicHenky1D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticPlaneStressUncoupledShear2DLaw", mElasticIsotropicPlaneStressUncoupledShear);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticityFactory", mSmallStrainIsotropicPlasticityFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticityFactory", mSmallStrainKinematicPlasticityFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageFactory", mSmallStrainIsotropicDamageFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedKelvin3D", mViscousGeneralizedKelvin3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedMaxwell3D", mViscousGeneralizedMaxwell3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("GenericSmallStrainViscoplasticity3D", mGenericSmallStrainViscoplasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("PlasticityIsotropicKinematicJ2Law", mPlasticityIsotropicKinematicJ2);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("WrinklingLinear2DLaw", mWrinklingLinear2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MultiLinearElastic1DLaw", mMultiLinearElastic1DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("MultiLinearIsotropicPlaneStress2D", mMultiLinearIsotropicPlaneStress2D);

    // Custom Constitutive laws
    // Serial-Parallel Rule Of Mixtures
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SerialParallelRuleOfMixturesLaw", mSerialParallelRuleOfMixturesLaw);

    // Anisotropic law
    KRATOS_REGISTER_CONSTITUTIVE_LAW("GenericAnisotropic3DLaw", mGenericAnisotropic3DLaw);

    /// Plasticity

    /* Small strain */
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesVonMises", mSmallStrainIsotropicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager", mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesTresca", mSmallStrainIsotropicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaVonMises", mSmallStrainIsotropicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaDruckerPrager", mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaTresca", mSmallStrainIsotropicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerVonMises", mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager", mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerTresca", mSmallStrainIsotropicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb", mSmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombVonMises", mSmallStrainIsotropicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb", mSmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager", mSmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombTresca", mSmallStrainIsotropicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaMohrCoulomb", mSmallStrainIsotropicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb", mSmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesVonMises", mSmallStrainKinematicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesDruckerPrager", mSmallStrainKinematicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesTresca", mSmallStrainKinematicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaVonMises", mSmallStrainKinematicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaDruckerPrager", mSmallStrainKinematicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaTresca", mSmallStrainKinematicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerVonMises", mSmallStrainKinematicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager", mSmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerTresca", mSmallStrainKinematicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesMohrCoulomb", mSmallStrainKinematicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombVonMises", mSmallStrainKinematicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb", mSmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager", mSmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombTresca", mSmallStrainKinematicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaMohrCoulomb", mSmallStrainKinematicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb", mSmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb);

    //Plastic Damage Model
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises", mSmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager", mSmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager);


    /* Finite strain */
    // Isotropic plasticity
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DVonMisesVonMises", mFiniteStrainIsotropicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb", mFiniteStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DVonMisesDruckerPrager", mFiniteStrainIsotropicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DVonMisesTresca", mFiniteStrainIsotropicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises", mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager", mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DModifiedMohrCoulombTresca", mFiniteStrainIsotropicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DTrescaVonMises", mFiniteStrainIsotropicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb", mFiniteStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DTrescaDruckerPrager", mFiniteStrainIsotropicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DTrescaTresca", mFiniteStrainIsotropicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DDruckerPragerVonMises", mFiniteStrainIsotropicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb", mFiniteStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DDruckerPragerDruckerPrager", mFiniteStrainIsotropicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DDruckerPragerTresca", mFiniteStrainIsotropicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DVonMisesMohrCoulomb", mFiniteStrainIsotropicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DMohrCoulombVonMises", mFiniteStrainIsotropicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb", mFiniteStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DMohrCoulombDruckerPrager", mFiniteStrainIsotropicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DMohrCoulombTresca", mFiniteStrainIsotropicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DTrescaMohrCoulomb", mFiniteStrainIsotropicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb", mFiniteStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb);
    // Kinematic plasticity
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DVonMisesVonMises", mFiniteStrainKinematicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb", mFiniteStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DVonMisesDruckerPrager", mFiniteStrainKinematicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DVonMisesTresca", mFiniteStrainKinematicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DModifiedMohrCoulombVonMises", mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager", mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DModifiedMohrCoulombTresca", mFiniteStrainKinematicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DTrescaVonMises", mFiniteStrainKinematicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb", mFiniteStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DTrescaDruckerPrager", mFiniteStrainKinematicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DTrescaTresca", mFiniteStrainKinematicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DDruckerPragerVonMises", mFiniteStrainKinematicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb", mFiniteStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DDruckerPragerDruckerPrager", mFiniteStrainKinematicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DDruckerPragerTresca", mFiniteStrainKinematicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DVonMisesMohrCoulomb", mFiniteStrainKinematicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DMohrCoulombVonMises", mFiniteStrainKinematicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DMohrCoulombMohrCoulomb", mFiniteStrainKinematicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DMohrCoulombDruckerPrager", mFiniteStrainKinematicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DMohrCoulombTresca", mFiniteStrainKinematicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DTrescaMohrCoulomb", mFiniteStrainKinematicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("FiniteStrainKinematicPlasticity3DDruckerPragerMohrCoulomb", mFiniteStrainKinematicPlasticity3DDruckerPragerMohrCoulomb);


    /// Damage
    /* Small strain */
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesVonMises", mSmallStrainIsotropicDamage3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesDruckerPrager", mSmallStrainIsotropicDamage3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesTresca", mSmallStrainIsotropicDamage3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca", mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaVonMises", mSmallStrainIsotropicDamage3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaDruckerPrager", mSmallStrainIsotropicDamage3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaTresca", mSmallStrainIsotropicDamage3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerVonMises", mSmallStrainIsotropicDamage3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager", mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerTresca", mSmallStrainIsotropicDamage3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineVonMises", mSmallStrainIsotropicDamage3DRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineDruckerPrager", mSmallStrainIsotropicDamage3DRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineTresca", mSmallStrainIsotropicDamage3DRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuVonMises", mSmallStrainIsotropicDamage3DSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuDruckerPrager", mSmallStrainIsotropicDamage3DSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuTresca", mSmallStrainIsotropicDamage3DSimoJuTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesMohrCoulomb", mSmallStrainIsotropicDamage3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombVonMises", mSmallStrainIsotropicDamage3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb", mSmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombTresca", mSmallStrainIsotropicDamage3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaMohrCoulomb", mSmallStrainIsotropicDamage3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb", mSmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineMohrCoulomb", mSmallStrainIsotropicDamage3DRankineMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuMohrCoulomb", mSmallStrainIsotropicDamage3DSimoJuMohrCoulomb);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesVonMises", mSmallStrainIsotropicDamage2DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesDruckerPrager", mSmallStrainIsotropicDamage2DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesTresca", mSmallStrainIsotropicDamage2DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises", mSmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombTresca", mSmallStrainIsotropicDamage2DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaVonMises", mSmallStrainIsotropicDamage2DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaDruckerPrager", mSmallStrainIsotropicDamage2DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaTresca", mSmallStrainIsotropicDamage2DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerVonMises", mSmallStrainIsotropicDamage2DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerDruckerPrager", mSmallStrainIsotropicDamage2DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerTresca", mSmallStrainIsotropicDamage2DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineVonMises", mSmallStrainIsotropicDamage2DRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineDruckerPrager", mSmallStrainIsotropicDamage2DRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineTresca", mSmallStrainIsotropicDamage2DRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuVonMises", mSmallStrainIsotropicDamage2DSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuDruckerPrager", mSmallStrainIsotropicDamage2DSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuTresca", mSmallStrainIsotropicDamage2DSimoJuTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesMohrCoulomb", mSmallStrainIsotropicDamage2DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombVonMises", mSmallStrainIsotropicDamage2DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb", mSmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage2DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombTresca", mSmallStrainIsotropicDamage2DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaMohrCoulomb", mSmallStrainIsotropicDamage2DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb", mSmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineMohrCoulomb", mSmallStrainIsotropicDamage2DRankineMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuMohrCoulomb", mSmallStrainIsotropicDamage2DSimoJuMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NonLocalElasticIsotropicDamage", mNonLocalElasticIsotropicDamage);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticAnisotropicDamage", mElasticAnisotropicDamage);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticIsotropicDamage3DNonLocalEquivalentStrain", mElasticIsotropicDamage3DNonLocalEquivalentStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticAnisotropicDamage3DNonLocalEquivalentStrain", mElasticAnisotropicDamage3DNonLocalEquivalentStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticAnisotropicDamage3DTwoNLEquivalentStrains", mElasticAnisotropicDamage3DTwoNLEquivalentStrains);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC", mElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC);

    // HCF (High Cycle Fatigue)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesVonMises", mSmallStrainHighCycleFatigue3DLawVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager", mSmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesTresca", mSmallStrainHighCycleFatigue3DLawVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaVonMises", mSmallStrainHighCycleFatigue3DLawTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaDruckerPrager", mSmallStrainHighCycleFatigue3DLawTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaTresca", mSmallStrainHighCycleFatigue3DLawTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerVonMises", mSmallStrainHighCycleFatigue3DLawDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager", mSmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerTresca", mSmallStrainHighCycleFatigue3DLawDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineVonMises", mSmallStrainHighCycleFatigue3DLawRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineDruckerPrager", mSmallStrainHighCycleFatigue3DLawRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineTresca", mSmallStrainHighCycleFatigue3DLawRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuVonMises", mSmallStrainHighCycleFatigue3DLawSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager", mSmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuTresca", mSmallStrainHighCycleFatigue3DLawSimoJuTresca);

    // d+d- laws (3D)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineRankine3D", mSmallStrainDplusDminusDamageRankineRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineSimoJu3D", mSmallStrainDplusDminusDamageRankineSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineVonMises3D", mSmallStrainDplusDminusDamageRankineVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineTresca3D", mSmallStrainDplusDminusDamageRankineTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineDruckerPrager3D", mSmallStrainDplusDminusDamageRankineDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuRankine3D", mSmallStrainDplusDminusDamageSimoJuRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuSimoJu3D", mSmallStrainDplusDminusDamageSimoJuSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuVonMises3D", mSmallStrainDplusDminusDamageSimoJuVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuTresca3D", mSmallStrainDplusDminusDamageSimoJuTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuDruckerPrager3D", mSmallStrainDplusDminusDamageSimoJuDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesRankine3D", mSmallStrainDplusDminusDamageVonMisesRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesSimoJu3D", mSmallStrainDplusDminusDamageVonMisesSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesVonMises3D", mSmallStrainDplusDminusDamageVonMisesVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesTresca3D", mSmallStrainDplusDminusDamageVonMisesTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesDruckerPrager3D", mSmallStrainDplusDminusDamageVonMisesDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaRankine3D", mSmallStrainDplusDminusDamageTrescaRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaSimoJu3D", mSmallStrainDplusDminusDamageTrescaSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaVonMises3D", mSmallStrainDplusDminusDamageTrescaVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaTresca3D", mSmallStrainDplusDminusDamageTrescaTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaDruckerPrager3D", mSmallStrainDplusDminusDamageTrescaDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerRankine3D", mSmallStrainDplusDminusDamageDruckerPragerRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerSimoJu3D", mSmallStrainDplusDminusDamageDruckerPragerSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerVonMises3D", mSmallStrainDplusDminusDamageDruckerPragerVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerTresca3D", mSmallStrainDplusDminusDamageDruckerPragerTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D", mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D", mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombRankine3D", mSmallStrainDplusDminusDamageMohrCoulombRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombSimoJu3D", mSmallStrainDplusDminusDamageMohrCoulombSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombVonMises3D", mSmallStrainDplusDminusDamageMohrCoulombVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombTresca3D", mSmallStrainDplusDminusDamageMohrCoulombTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D", mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineMohrCoulomb3D", mSmallStrainDplusDminusDamageRankineMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuMohrCoulomb3D", mSmallStrainDplusDminusDamageSimoJuMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesMohrCoulomb3D", mSmallStrainDplusDminusDamageVonMisesMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaMohrCoulomb3D", mSmallStrainDplusDminusDamageTrescaMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D", mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D);

    // d+d- laws (2D)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombRankine2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombTresca2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineRankine2D", mSmallStrainDplusDminusDamageRankineRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineSimoJu2D", mSmallStrainDplusDminusDamageRankineSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineVonMises2D", mSmallStrainDplusDminusDamageRankineVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineTresca2D", mSmallStrainDplusDminusDamageRankineTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineDruckerPrager2D", mSmallStrainDplusDminusDamageRankineDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuRankine2D", mSmallStrainDplusDminusDamageSimoJuRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuSimoJu2D", mSmallStrainDplusDminusDamageSimoJuSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuVonMises2D", mSmallStrainDplusDminusDamageSimoJuVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuTresca2D", mSmallStrainDplusDminusDamageSimoJuTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuDruckerPrager2D", mSmallStrainDplusDminusDamageSimoJuDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesRankine2D", mSmallStrainDplusDminusDamageVonMisesRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesSimoJu2D", mSmallStrainDplusDminusDamageVonMisesSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesVonMises2D", mSmallStrainDplusDminusDamageVonMisesVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesTresca2D", mSmallStrainDplusDminusDamageVonMisesTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesDruckerPrager2D", mSmallStrainDplusDminusDamageVonMisesDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaRankine2D", mSmallStrainDplusDminusDamageTrescaRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaSimoJu2D", mSmallStrainDplusDminusDamageTrescaSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaVonMises2D", mSmallStrainDplusDminusDamageTrescaVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaTresca2D", mSmallStrainDplusDminusDamageTrescaTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaDruckerPrager2D", mSmallStrainDplusDminusDamageTrescaDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerRankine2D", mSmallStrainDplusDminusDamageDruckerPragerRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerSimoJu2D", mSmallStrainDplusDminusDamageDruckerPragerSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerVonMises2D", mSmallStrainDplusDminusDamageDruckerPragerVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerTresca2D", mSmallStrainDplusDminusDamageDruckerPragerTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerDruckerPrager2D", mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D", mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombRankine2D", mSmallStrainDplusDminusDamageMohrCoulombRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombSimoJu2D", mSmallStrainDplusDminusDamageMohrCoulombSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombVonMises2D", mSmallStrainDplusDminusDamageMohrCoulombVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombTresca2D", mSmallStrainDplusDminusDamageMohrCoulombTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D", mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineMohrCoulomb2D", mSmallStrainDplusDminusDamageRankineMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuMohrCoulomb2D", mSmallStrainDplusDminusDamageSimoJuMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesMohrCoulomb2D", mSmallStrainDplusDminusDamageVonMisesMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaMohrCoulomb2D", mSmallStrainDplusDminusDamageTrescaMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D", mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("DamageDPlusDMinusPlaneStressMasonry2DLaw", mDamageDPlusDMinusPlaneStressMasonry2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("DamageDPlusDMinusMasonry3DLaw", mDamageDPlusDMinusMasonry3DLaw);

    // Orthotropic Damage
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageRankine3D", mSmallStrainOrthotropicDamageRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageVonMises3D", mSmallStrainOrthotropicDamageVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageDruckerPrager3D", mSmallStrainOrthotropicDamageDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageTresca3D", mSmallStrainOrthotropicDamageTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageMohrCoulomb3D", mSmallStrainOrthotropicDamageMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageModifiedMohrCoulomb3D", mSmallStrainOrthotropicDamageModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageSimoJu3D", mSmallStrainOrthotropicDamageSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageRankine2D", mSmallStrainOrthotropicDamageRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageVonMises2D", mSmallStrainOrthotropicDamageVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageDruckerPrager2D", mSmallStrainOrthotropicDamageDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageTresca2D", mSmallStrainOrthotropicDamageTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageMohrCoulomb2D", mSmallStrainOrthotropicDamageMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageModifiedMohrCoulomb2D", mSmallStrainOrthotropicDamageModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageSimoJu2D", mSmallStrainOrthotropicDamageSimoJu2D);

    // Rules of mixtures
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ParallelRuleOfMixturesLaw2D", mParallelRuleOfMixturesLaw2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ParallelRuleOfMixturesLaw3D", mParallelRuleOfMixturesLaw3D);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("AssociativePlasticDamageModel3DVonMisesVonMises", mAssociativePlasticDamageModel3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("AssociativePlasticDamageModel3DDruckerPragerDruckerPrager", mAssociativePlasticDamageModel3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("AssociativePlasticDamageModel3DModifiedMohrCoulombModifiedMohrCoulomb", mAssociativePlasticDamageModel3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("AssociativePlasticDamageModel3DRankineRankine", mAssociativePlasticDamageModel3DRankineRankine);

    // Constitutive laws variables

    //Faituge variables
    KRATOS_REGISTER_VARIABLE(HIGH_CYCLE_FATIGUE_COEFFICIENTS)
    KRATOS_REGISTER_VARIABLE(STRESS_LIMITS)
    KRATOS_REGISTER_VARIABLE(HARDENING_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(FATIGUE_REDUCTION_FACTOR)

    KRATOS_REGISTER_VARIABLE(LOCAL_NUMBER_OF_CYCLES)
    KRATOS_REGISTER_VARIABLE(WOHLER_STRESS)
    KRATOS_REGISTER_VARIABLE(REVERSION_FACTOR_RELATIVE_ERROR)
    KRATOS_REGISTER_VARIABLE(MAX_STRESS_RELATIVE_ERROR)
    KRATOS_REGISTER_VARIABLE(MAX_STRESS)
    KRATOS_REGISTER_VARIABLE(THRESHOLD_STRESS)
    KRATOS_REGISTER_VARIABLE(CYCLE_INDICATOR)
    KRATOS_REGISTER_VARIABLE(CYCLES_TO_FAILURE)
    KRATOS_REGISTER_VARIABLE(TIME_INCREMENT)
    KRATOS_REGISTER_VARIABLE(DAMAGE_ACTIVATION)
    KRATOS_REGISTER_VARIABLE(PREVIOUS_CYCLE);
    KRATOS_REGISTER_VARIABLE(CYCLE_PERIOD)
    KRATOS_REGISTER_VARIABLE(ADVANCE_STRATEGY_APPLIED);

    // Constitutive laws variables
    KRATOS_REGISTER_VARIABLE(LAYER_EULER_ANGLES)
    KRATOS_REGISTER_VARIABLE(INELASTIC_FLAG)
    KRATOS_REGISTER_VARIABLE(INFINITY_YIELD_STRESS)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_TENSION)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(DILATANCY_ANGLE)
    KRATOS_REGISTER_VARIABLE(SOFTENING_TYPE)
    KRATOS_REGISTER_VARIABLE(SOFTENING_TYPE_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(HARDENING_CURVE)
    KRATOS_REGISTER_VARIABLE(HARDENING_MODULUS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(HARDENING_MODULUS_TENSION)
    KRATOS_REGISTER_VARIABLE(MAX_NUMBER_NL_CL_ITERATIONS)
    KRATOS_REGISTER_VARIABLE(VISCOUS_PARAMETER)
    KRATOS_REGISTER_VARIABLE(DELAY_TIME)
    KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS)
    KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS_POSITION)
    KRATOS_REGISTER_VARIABLE(PLASTIC_DISSIPATION_LIMIT_LINEAR_SOFTENING)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS)
    KRATOS_REGISTER_VARIABLE(FRICTION_ANGLE)
    KRATOS_REGISTER_VARIABLE(COHESION)
    KRATOS_REGISTER_VARIABLE(DISSIPATION)
    KRATOS_REGISTER_VARIABLE(DAMAGE)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MATRIX)
    KRATOS_REGISTER_VARIABLE(DAMAGE_FIBER)
    KRATOS_REGISTER_VARIABLE(THRESHOLD)
    KRATOS_REGISTER_VARIABLE(INTEGRATED_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_TENSOR)
    KRATOS_REGISTER_VARIABLE(CURVE_FITTING_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(TANGENCY_REGION2)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_INDICATORS)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(KINEMATIC_PLASTICITY_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(KINEMATIC_HARDENING_TYPE)
    KRATOS_REGISTER_VARIABLE(CONSIDER_PERTURBATION_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(TANGENT_OPERATOR_ESTIMATION)
    KRATOS_REGISTER_VARIABLE(TENSION_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(COMPRESSION_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(TENSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(COMPRESSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_TENSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_COMPRESSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_TENSOR_FIBER)
    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_TENSOR_MATRIX)
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_TENSOR_MATRIX)
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_TENSOR_FIBER)
    KRATOS_REGISTER_VARIABLE(EXPONENTIAL_SATURATION_YIELD_STRESS)
    KRATOS_REGISTER_VARIABLE(ACCUMULATED_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(BACK_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(BACK_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(MULTI_LINEAR_ELASTICITY_MODULI)
    KRATOS_REGISTER_VARIABLE(MULTI_LINEAR_ELASTICITY_STRAINS)
    KRATOS_REGISTER_VARIABLE(STRAIN_DAMAGE_CURVE)
    KRATOS_REGISTER_VARIABLE(STRESS_DAMAGE_CURVE)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_STRESS_VECTOR_PLASTICITY_POINT_CURVE)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR_PLASTICITY_POINT_CURVE)
    KRATOS_REGISTER_VARIABLE(DAMAGE_VECTOR)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_STRAINS_TC)
    KRATOS_REGISTER_VARIABLE(DAMAGE_TENSOR)

    // Some variables related with SP
    KRATOS_REGISTER_VARIABLE(SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE)

    // D+D- Damage Constitutive laws variables
    KRATOS_REGISTER_VARIABLE(DAMAGE_TENSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(THRESHOLD_TENSION)
    KRATOS_REGISTER_VARIABLE(THRESHOLD_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS_TENSION)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_DAMAGE_PROCESS)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EULER_ANGLES)
    KRATOS_REGISTER_VARIABLE(OGDEN_BETA_1)
    KRATOS_REGISTER_VARIABLE(OGDEN_BETA_2)
    KRATOS_REGISTER_VARIABLE(DAMAGE_THRESHOLD_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_THRESHOLD_TENSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL_PARAMETER_BETA1_TENSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL_PARAMETER_BETA2_TENSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION)

    // D+D- Damage Constitutive laws variables, additional Masonry 2D & 3D
    KRATOS_REGISTER_VARIABLE(DAMAGE_ONSET_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(BIAXIAL_COMPRESSION_MULTIPLIER)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_TENSION)
    KRATOS_REGISTER_VARIABLE(RESIDUAL_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C1)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C2)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C3)
    KRATOS_REGISTER_VARIABLE(YIELD_STRAIN_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(SHEAR_COMPRESSION_REDUCTOR)
    KRATOS_REGISTER_VARIABLE(TRIAXIAL_COMPRESSION_COEFFICIENT)
	KRATOS_REGISTER_VARIABLE(INTEGRATION_IMPLEX)
    KRATOS_REGISTER_VARIABLE(TENSION_YIELD_MODEL)
    KRATOS_REGISTER_VARIABLE(PLASTIC_DAMAGE_PROPORTION)

    // For anisotropy + orthotrophy
    // The ratios between the yield strength in the isotropic space and the anisotropic space
    // at each direction in local coordinates ratio
    KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(ISOTROPIC_ANISOTROPIC_YIELD_RATIO);
    KRATOS_REGISTER_SYMMETRIC_3D_TENSOR_VARIABLE_WITH_COMPONENTS(ORTHOTROPIC_ELASTIC_CONSTANTS);

}
}  // namespace Kratos.
