
#include <iostream>
#include "DEM_D_Hertz_viscous_Coulomb_nestle_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Hertz_viscous_Coulomb_nestle::Clone() const {
        
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Hertz_viscous_Coulomb_nestle(*this));
        return p_clone;
    }

    void DEM_D_Hertz_viscous_Coulomb_nestle::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) const {
        
        if (verbose) std::cout << "\nAssigning DEM_D_Hertz_viscous_Coulomb_nestle to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
    }

} //namespace Kratos
