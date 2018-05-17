#include "DEM_D_Linear_HighStiffness_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

        DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_HighStiffness::Clone() const {

        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_HighStiffness(*this));
        return p_clone;
    }

    void DEM_D_Linear_HighStiffness::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_HighStiffness_CL to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    double DEM_D_Linear_HighStiffness::CalculateNormalForce(const double indentation) {
        return DEM_D_Linear_viscous_Coulomb::CalculateNormalForce(indentation) * 100;
    }

} // namespace Kratos
