
#if !defined(DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED)
#define DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED

#include "DEM_D_Hertz_viscous_Coulomb_CL.h"

namespace Kratos {
    
    class SphericParticle;

    class DEM_D_Hertz_viscous_Coulomb_nestle : public DEM_D_Hertz_viscous_Coulomb {
    
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Hertz_viscous_Coulomb_nestle);

        DEM_D_Hertz_viscous_Coulomb_nestle() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        ~DEM_D_Hertz_viscous_Coulomb_nestle() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;
    };
} // namespace Kratos

#endif // DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED defined
