
#include <iostream>
#include "dem_kdem_2d_cl.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEM2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEM2D(*this));
        return p_clone;
    }

    void DEM_KDEM2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_KDEM2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEM2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = radius * other_radius / radius_sum;
        calculation_area = 2.0 * equiv_radius;
        KRATOS_CATCH("")
    }

} /* namespace Kratos.*/
