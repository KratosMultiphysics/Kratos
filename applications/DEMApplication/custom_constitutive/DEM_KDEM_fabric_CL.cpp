
#include <iostream>
#include "DEM_KDEM_fabric_CL.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_KDEMFabric::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_KDEMFabric(*this));
        return p_clone;
    }

    void DEM_KDEMFabric::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        KRATOS_INFO("DEM") << "Assigning DEM_KDEMFabric to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }


    void DEM_KDEMFabric::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force,
                                          double calculation_area, BoundedMatrix<double, 3, 3>* mSymmStressTensor, SphericContinuumParticle* element1,
                                          SphericContinuumParticle* element2, const ProcessInfo& r_process_info, const int i_neighbor_count, const double indentation) {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

} // namespace Kratos
