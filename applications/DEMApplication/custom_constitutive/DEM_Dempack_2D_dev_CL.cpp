// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "DEM_Dempack_2D_dev_CL.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::Pointer DEM_Dempack2D_dev::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEM_Dempack2D_dev(*this));
        return p_clone;
    }

    void DEM_Dempack2D_dev::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_Dempack2D_dev to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

    void DEM_Dempack2D_dev::CalculateContactArea(double radius, double other_radius, double& calculation_area) {

        double rmin = radius;
        if (other_radius < radius) rmin = other_radius;
        calculation_area = 2*rmin;
    }

} /* namespace Kratos.*/
