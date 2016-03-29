// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "DEM_application.h"
#include "DEM_Dempack_2D_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    void DEM_Dempack2D::Initialize() {
    KRATOS_TRY  
        mHistoryMaxInd              = 0.0; //maximum indentation achieved
        mHistoryMaxForce            = 0.0; //maximum force achieved
        mHistoryDamage              = 0.0; //cumulated_damage
        mHistoryDegradation         = 1.0; //degradation factor for G reducing in Dempack;
        mHistoryDisp                = 0.0; //displacement;
        mHistoryShearFlag           = 0.0; //superado el limite de cortante;  
    KRATOS_CATCH("")  
    }

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack2D(*this));
        return p_clone;
    }

    void DEM_Dempack2D::SetConstitutiveLawInProperties(Properties::Pointer pProp) const {
        std::cout << "Assigning DEM_Dempack2D to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_Dempack2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {
        
        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = 2*rmin;
    }

} /* namespace Kratos.*/
