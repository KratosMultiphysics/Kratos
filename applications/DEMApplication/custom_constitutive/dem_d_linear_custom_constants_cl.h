// Authors: S. Latorre (CIMNE)
// Date: April 2016

#if !defined(DEM_D_LINEAR_CUSTOM_CONSTANTS_CL_H_INCLUDED)
#define DEM_D_LINEAR_CUSTOM_CONSTANTS_CL_H_INCLUDED

#include "DEM_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_Custom_Constants : public DEM_D_Linear_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_Custom_Constants);

        DEM_D_Linear_Custom_Constants() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        ~DEM_D_Linear_Custom_Constants() {}

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;
        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
            //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
            //rSerializer.load("MyMemberName",myMember);
        }
    };
} // namespace Kratos
#endif // DEM_D_LINEAR_CUSTOM_CONSTANTS_CL_H_INCLUDED defined
