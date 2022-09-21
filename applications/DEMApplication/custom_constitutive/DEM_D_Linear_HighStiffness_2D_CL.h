#if !defined(DEM_D_LINEAR_HIGHSTIFFNESS_2D_CL_H_INCLUDED)
#define DEM_D_LINEAR_HIGHSTIFFNESS_2D_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_D_Linear_viscous_Coulomb_2D_CL.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_HighStiffness_2D : public DEM_D_Linear_viscous_Coulomb2D {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_HighStiffness_2D);

        DEM_D_Linear_HighStiffness_2D() {}

        ~DEM_D_Linear_HighStiffness_2D() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void Check(Properties::Pointer pProp) const override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;
    
    protected:

    };

} /* namespace Kratos.*/

#endif /* DEM_D_LINEAR_HIGHSTIFFNESS_2D_CL_H_INCLUDED  defined */
