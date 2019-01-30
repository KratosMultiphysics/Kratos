
#include <iostream>
#include "custom_elements/spheric_continuum_particle.h"
#include "DEM_KDEM_fabric_2D_CL.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEMFabric2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEMFabric2D(*this));
        return p_clone;
    }

    void DEM_KDEMFabric2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_KDEMFabric2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_KDEMFabric2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        KRATOS_TRY
        double radius_sum = radius + other_radius;
        double equiv_radius = radius * other_radius / radius_sum;
        calculation_area = 2.0 * equiv_radius;
        KRATOS_CATCH("")
    }

} // namespace Kratos
