#include "DEM_D_Linear_HighStiffness_2D_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

        DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_HighStiffness_2D::Clone() const {

        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_HighStiffness_2D(*this));
        return p_clone;
    }

    void DEM_D_Linear_HighStiffness_2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Linear_HighStiffness_2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_D_Linear_HighStiffness_2D::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        DEM_D_Linear_viscous_Coulomb2D::InitializeContact(element1, element2, indentation);
        const double kn_augmenter = 1.0;
        mKn *= kn_augmenter;
    }

    void DEM_D_Linear_HighStiffness_2D::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        DEM_D_Linear_viscous_Coulomb2D::InitializeContactWithFEM(element, wall, indentation, ini_delta);
        const double kn_augmenterFEM = 1.0;
        mKn *= kn_augmenterFEM;
    }
} // namespace Kratos
