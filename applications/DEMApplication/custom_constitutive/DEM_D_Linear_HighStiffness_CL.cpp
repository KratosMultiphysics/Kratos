#include "DEM_D_Linear_HighStiffness_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

        DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_HighStiffness::Clone() const {

        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_HighStiffness(*this));
        return p_clone;
    }

    void DEM_D_Linear_HighStiffness::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_HighStiffness to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_D_Linear_HighStiffness::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        DEM_D_Linear_viscous_Coulomb::InitializeContact(element1, element2, indentation);
        const double kn_augmenter = 5.0;
        mKn *= kn_augmenter;
    }

    void DEM_D_Linear_HighStiffness::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        DEM_D_Linear_viscous_Coulomb::InitializeContactWithFEM(element, wall, indentation, ini_delta);
        const double kn_augmenterFEM = 5.0;
        mKn *= kn_augmenterFEM;
    }
} // namespace Kratos
