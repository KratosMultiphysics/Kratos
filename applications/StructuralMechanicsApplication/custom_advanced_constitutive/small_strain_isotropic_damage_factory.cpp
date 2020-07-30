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
#include "custom_advanced_constitutive/small_strain_isotropic_damage_factory.h"
#include "custom_advanced_constitutive/generic_small_strain_isotropic_damage.h"

namespace Kratos
{
ConstitutiveLaw::Pointer SmallStrainIsotropicDamageFactory::Create(Kratos::Parameters NewParameters) const
{
    const std::string law_type = NewParameters.Has("law_type") ? NewParameters["law_type"].GetString() : "3D";
    const std::string& yield = NewParameters["yield_surface"].GetString();
    const std::string& potential = NewParameters["plastic_potential"].GetString();

    const std::string& name = "SmallStrainIsotropicDamage" + law_type + yield + potential;
    return KratosComponents<ConstitutiveLaw>::Get(name).Clone();
}
} // namespace Kratos
