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

    void DEM_D_Linear_HighStiffness_2D::Check(Properties::Pointer pProp) const {

        DEM_D_Linear_viscous_Coulomb2D::Check(pProp);
        
        if (!pProp->Has(STIFFNESS_FACTOR)) {
            KRATOS_WARNING("DEM") << std::endl;
            KRATOS_WARNING("DEM") << "WARNING: Variable STIFFNESS_FACTOR should be present in the Properties when using DEM_D_Linear_HighStiffness_2D. A default value of 5.0 was assigned." << std::endl;
            KRATOS_WARNING("DEM") << std::endl;
            pProp->GetValue(STIFFNESS_FACTOR) = 5.0;
        }
    }

    void DEM_D_Linear_HighStiffness_2D::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        
        DEM_D_Linear_viscous_Coulomb2D::InitializeContact(element1, element2, indentation);
        Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
        mKn *= properties_of_this_contact[STIFFNESS_FACTOR];
    }

    void DEM_D_Linear_HighStiffness_2D::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        
        DEM_D_Linear_viscous_Coulomb2D::InitializeContactWithFEM(element, wall, indentation, ini_delta);
        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
        mKn *= properties_of_this_contact[STIFFNESS_FACTOR];
    }
} // namespace Kratos
