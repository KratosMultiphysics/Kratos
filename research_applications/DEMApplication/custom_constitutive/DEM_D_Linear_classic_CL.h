/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

#if !defined (DEM_D_LINEAR_CLASSIC_CL_H_INCLUDE)
#define DEM_D_LINEAR_CLASSIC_CL_H_INCLUDE

#include <string>
#include <iostream>
#include "DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_classic : public DEM_D_Linear_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_classic);

        DEM_D_Linear_classic() {}

        ~DEM_D_Linear_classic() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

        double GetTangentialStiffness() override;
        
    protected:

    };
} // namespace Kratos

#endif // DEM_D_LINEAR_CLASSIC_CL_H_INCLUDE  defined