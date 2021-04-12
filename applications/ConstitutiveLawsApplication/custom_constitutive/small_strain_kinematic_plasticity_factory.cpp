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

    KRATOS_ERROR_IF(yield == "SimoJu") << "SimoJu yield surface not available in plasticity " << yield << std::endl;
    KRATOS_ERROR_IF(yield == "Rankine") << "Rankine yield surface not available in plasticity " << yield << std::endl;

    const std::string name = "SmallStrainIsotropicKinematicPlasticity" + law_type + yield + potential;
    return KratosComponents<ConstitutiveLaw>::Get(name).Clone();
}

} // namespace Kratos
