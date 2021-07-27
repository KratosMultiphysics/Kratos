// System includes
#include <string>
#include <iostream>
#include <cmath>

// Project includes
#include "DEM_KDEM_with_damage_parallel_bond_Hertz_2D_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM_with_damage_parallel_bond_Hertz_2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM_with_damage_parallel_bond_Hertz_2D(*this));
        return p_clone;
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz_2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond_Hertz_2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz_2D::SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEM_with_damage_parallel_bond_Hertz_2D to Properties " << pProp->Id() <<" with given parameters"<< std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        TransferParametersToProperties(parameters, pProp);
        this->Check(pProp);
    }

    void DEM_KDEM_with_damage_parallel_bond_Hertz_2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = radius * other_radius / radius_sum; //double equiv_radius = 0.5 * radius_sum;
        calculation_area = 2.0 * equiv_radius * 1.0; // This last 1.0 is simply the unitary thickness of 2D, in meters.
        KRATOS_CATCH("")
    }

} // namespace Kratos
