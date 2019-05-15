// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//  Collaborator:    Alejandro Cornejo & Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "custom_constitutive/finite_strain_isotropic_plasticity_factory.h"
#include "custom_constitutive/generic_finite_strain_isotropic_plasticity.h"

namespace Kratos
{
ConstitutiveLaw::Pointer FiniteStrainIsotropicPlasticityFactory::Create(Kratos::Parameters NewParameters) const
{
    const std::string law_type = NewParameters.Has("law_type") ? NewParameters["law_type"].GetString() : "3D";
    const std::string& elastic_behaviour = NewParameters["elastic_behaviour"].GetString();
    const std::string& yield = NewParameters["yield_surface"].GetString();
    const std::string& potential = NewParameters["plastic_potential"].GetString();

    KRATOS_ERROR_IF(yield == "SimoJu") << "SimoJu yield surface not available in plasticity " << yield << std::endl;
    KRATOS_ERROR_IF(yield == "Rankine") << "Rankine yield surface not available in plasticity " << yield << std::endl;

    if (elastic_behaviour == "HyperElasticIsotropicNeoHookean3D") {
        if (yield == "VonMises") {
            if (potential == "VonMises") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
            } else if (potential == "ModifiedMohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
            } else if (potential == "Tresca") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
            } else if (potential == "DruckerPrager") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
            } else if (potential == "MohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
            } else {
                KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
            }
        } else if (yield == "ModifiedMohrCoulomb") {
            if (potential == "VonMises") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
            } else if (potential == "ModifiedMohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
            } else if (potential == "Tresca") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
            } else if (potential == "DruckerPrager") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
            } else {
                KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
            }
        } else if (yield == "Tresca") {
            if (potential == "VonMises") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
            } else if (potential == "ModifiedMohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
            } else if (potential == "Tresca") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
            } else if (potential == "DruckerPrager") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
            } else if (potential == "MohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
            } else {
                KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
            }
        } else if (yield == "DruckerPrager") {
            if (potential == "VonMises") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
            } else if (potential == "ModifiedMohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>().Clone();
            } else if (potential == "Tresca") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
            } else if (potential == "DruckerPrager") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
            } else if (potential == "MohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
            } else {
                KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
            }
        } else if (yield == "MohrCoulomb") {
            if (potential == "VonMises") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>().Clone();
            } else if (potential == "Tresca") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>().Clone();
            } else if (potential == "DruckerPrager") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>().Clone();
            } else if (potential == "MohrCoulomb") {
                return GenericFiniteStrainIsotropicPlasticity<HyperElasticIsotropicNeoHookean3D, GenericFiniteStrainConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>().Clone();
            } else {
                KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
            }
        } else {
            KRATOS_ERROR << "The combination of the yield surface <" << yield << "> with the plastic potential <" << potential << "> is not available" << std::endl;
        }
    }

}

} // namespace Kratos
