#include "DEM_D_Linear_HighStiffness_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_D_Linear_HighStiffness_CL::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_HighStiffness_CL to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    double DEM_D_Linear_HighStiffness_CL::CalculateNormalForce(const double indentation) {
        return DEM_D_Linear_viscous_Coulomb::CalculateNormalForce(indentation) * 100;
    }

} // namespace Kratos
