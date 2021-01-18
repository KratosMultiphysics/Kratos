
#if !defined(DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED)
#define DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED

#include "DEM_D_Hertz_viscous_Coulomb_CL.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Hertz_viscous_Coulomb_Nestle : public DEM_D_Hertz_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Hertz_viscous_Coulomb_Nestle);

        DEM_D_Hertz_viscous_Coulomb_Nestle() {}

        ~DEM_D_Hertz_viscous_Coulomb_Nestle() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void CalculateViscoDampingForce(double LocalRelVel[3], double ViscoDampingLocalContactForce[3], SphericParticle* const element1, SphericParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3], double ViscoDampingLocalContactForce[3], SphericParticle* const element, Condition* const wall);
    };
} // namespace Kratos

#endif // DEM_D_HERTZ_VISCOUS_COULOMB_NESTLE_H_INCLUDED defined
