// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main author:     Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "generic_finite_strain_kinematic_plasticity_plane_stress.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_kinematic_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/plane_stress_von_mises_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/plane_stress_von_mises_plastic_potential.h"

namespace Kratos
{

template class GenericFiniteStrainKinematicPlasticityPlaneStress<GenericConstitutiveLawIntegratorKinematicPlasticity<PlaneStressVonMisesYieldSurface<PlaneStressVonMisesPlasticPotential<3>>>>;

} // namespace Kratos
