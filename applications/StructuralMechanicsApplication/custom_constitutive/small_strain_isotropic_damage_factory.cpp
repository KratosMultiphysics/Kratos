// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/small_strain_isotropic_damage_factory.h"
#include "custom_constitutive/generic_small_strain_isotropic_damage.h"

namespace Kratos
{
ConstitutiveLaw::Pointer SmallStrainIsotropicDamageFactory::Create(Kratos::Parameters NewParameters) const
{
    const std::string law_type = NewParameters.Has("law_type") ? NewParameters["law_type"].GetString() : "3D";
    const std::string& yield = NewParameters["yield_surface"].GetString();
    const std::string& potential = NewParameters["plastic_potential"].GetString();


    if (law_type == "3D") {
        if (yield == "VonMises") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "Tresca") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "DruckerPrager") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "MohrCoulomb") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "Rankine") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (yield == "SimoJu") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "This yield surface is not available: " << yield << std::endl;
        }
    } else {
        if (yield == "VonMises") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "Tresca") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "DruckerPrager") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "MohrCoulomb") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "Rankine") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else if (yield == "SimoJu") {
            return GenericSmallStrainIsotropicDamage<GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<3>>>>().Clone();
        } else {
            KRATOS_ERROR << "This yield surface is not available: " << yield << std::endl;
        }
    }
}
} // namespace Kratos
