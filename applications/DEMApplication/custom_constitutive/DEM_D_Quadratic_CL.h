/////////////////////////////////////////////////
// Author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

//TODO: // Here this Quadratic model is established for Parallel Bond Model only.
        // It CAN NOT BE USED AS AN INDEPENDENT CONTACT MODEL UNTIL IT IS COMPLETED.

#if !defined(DEM_D_QUADRATIC_CL_H_INCLUDED)
#define DEM_D_QUADRATIC_CL_H_INCLUDED

//project include
#include <string>
#include <iostream>
#include "DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos{

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Quadratic : public DEM_D_Linear_viscous_Coulomb {

    public:

        //using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Quadratic);

        DEM_D_Quadratic () {}

        ~DEM_D_Quadratic () {}

        std::string GetTypeOfLaw() override;

        void Check(Properties::Pointer pProp) const override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        double CalculateNormalForce(const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

        double GetTangentialStiffness() override;

    private:

    }; //CLASS DEM_D_QUADRATIC

} //namespace Kratos

#endif //DEM_D_QUADRATIC_CL_H_INCLUDED defined