
#if !defined(DEM_KDEM_2D_CL_H_INCLUDED)
#define DEM_KDEM_2D_CL_H_INCLUDED

#include "DEM_KDEM_CL.h"

namespace Kratos {

    class DEM_KDEM2D : public DEM_KDEM {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM2D);

        DEM_KDEM2D() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_KDEM2D() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const;
        void CalculateContactArea(double radius, double other_radius, double& calculation_area);

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
            //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
            //rSerializer.load("MyMemberName",myMember);
        }
    };
} // namespace Kratos
#endif // DEM_KDEM_2D_H_INCLUDED defined
