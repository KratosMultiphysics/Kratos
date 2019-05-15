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
#include "custom_constitutive/small_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/generic_small_strain_isotropic_plasticity.h"

namespace Kratos
{
ConstitutiveLaw::Pointer SmallStrainIsotropicPlasticityFactory::Create(Kratos::Parameters NewParameters) const
{
    const std::string law_type = NewParameters.Has("law_type") ? NewParameters["law_type"].GetString() : "3D";
    const std::string& yield = NewParameters["yield_surface"].GetString();
    const std::string& potential = NewParameters["plastic_potential"].GetString();

    KRATOS_ERROR_IF(yield == "SimoJu") << "SimoJu yield surface not available in plasticity "   << yield << std::endl;
    KRATOS_ERROR_IF(yield == "Rankine") << "Rankine yield surface not available in plasticity " << yield << std::endl;

    if (yield == "VonMises") {
        if (potential == "VonMises") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "ModifiedMohrCoulomb") {
        if (potential == "VonMises") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "Tresca") {
        if (potential == "VonMises") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "DruckerPrager") {
        if (potential == "VonMises") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "MohrCoulomb") {
        if (potential == "VonMises") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainIsotropicPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else {
        KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
    }
}

} // namespace Kratos
