// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_parallel_bond_2D_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage_parallel_bond_2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage_parallel_bond_2D(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage_parallel_bond_2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = radius * other_radius / radius_sum;
        calculation_area = 2.0 * equiv_radius;
        KRATOS_CATCH("")
    }

} // namespace Kratos
