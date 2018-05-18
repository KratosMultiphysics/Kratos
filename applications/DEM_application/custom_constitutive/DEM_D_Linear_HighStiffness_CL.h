#if !defined(DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED)
#define DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_HighStiffness_CL : public DEM_D_Linear_viscous_Coulomb {

    public:

        using DEM_D_Linear_viscous_Coulomb::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_HighStiffness_CL);

        DEM_D_Linear_HighStiffness_CL() {}

        ~DEM_D_Linear_HighStiffness_CL() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        DEM_D_Linear_viscous_Coulomb::Pointer Clone() const override;

        double CalculateNormalForce(const double indentation) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEM_D_Linear_viscous_Coulomb)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEM_D_Linear_viscous_Coulomb)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_D_LINEAR_HIGHSTIFFNESS_CL_H_INCLUDED  defined */
