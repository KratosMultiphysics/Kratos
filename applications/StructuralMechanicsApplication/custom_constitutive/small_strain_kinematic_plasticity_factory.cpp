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
#include "custom_constitutive/small_strain_kinematic_plasticity_factory.h"
#include "custom_constitutive/generic_small_strain_kinematic_plasticity.h"

namespace Kratos
{
ConstitutiveLaw::Pointer SmallStrainKinematicPlasticityFactory::Create(Kratos::Parameters NewParameters) const
{
    const std::string law_type = NewParameters.Has("law_type") ? NewParameters["law_type"].GetString() : "3D";
    const std::string& yield = NewParameters["yield_surface"].GetString();
    const std::string& potential = NewParameters["plastic_potential"].GetString();

    KRATOS_ERROR_IF(yield == "SimoJu") << "SimoJu yield surface not available in plasticity "   << yield << std::endl;
    KRATOS_ERROR_IF(yield == "Rankine") << "Rankine yield surface not available in plasticity " << yield << std::endl;

    if (yield == "VonMises") {
        if (potential == "VonMises") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "ModifiedMohrCoulomb") {
        if (potential == "VonMises") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "Tresca") {
        if (potential == "VonMises") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "DruckerPrager") {
        if (potential == "VonMises") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "ModifiedMohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else if (yield == "MohrCoulomb") {
        if (potential == "VonMises") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
        } else if (potential == "Tresca") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
        } else if (potential == "DruckerPrager") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
        } else if (potential == "MohrCoulomb") {
            return GenericSmallStrainKinematicPlasticity<GenericConstitutiveLawIntegratorKinematicPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    } else {
        KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
    }
}

} // namespace Kratos
