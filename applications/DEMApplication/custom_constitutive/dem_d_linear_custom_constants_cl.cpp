// Authors: S. Latorre (CIMNE)
// Date: April 2016

#include "dem_d_linear_custom_constants_cl.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_Custom_Constants::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_Custom_Constants(*this));
        return p_clone;
    }

    void DEM_D_Linear_Custom_Constants::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_Custom_Constants to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_D_Linear_Custom_Constants::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {

        //Normal and Tangent elastic constants
        mKn = 0.5 * (element1->GetParticleKNormal() + element2->GetParticleKNormal());
        mKt = 0.5 * (element1->GetParticleKTangential() + element2->GetParticleKTangential());
    }

    void DEM_D_Linear_Custom_Constants::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {

        //Normal and Tangent elastic constants
        mKn = element->GetParticleKNormal();
        mKt = element->GetParticleKTangential();
    }

} // namespace Kratos
