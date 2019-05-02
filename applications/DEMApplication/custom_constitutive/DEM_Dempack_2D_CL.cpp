// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "DEM_Dempack_2D_CL.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack2D::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack2D(*this));
        return p_clone;
    }

    void DEM_Dempack2D::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_Dempack2D to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_Dempack2D::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = 2*rmin;
    }

} /* namespace Kratos.*/
