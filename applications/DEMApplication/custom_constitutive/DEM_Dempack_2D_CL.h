
#if !defined(DEM_DEMPACK_2D_CL_H_INCLUDED)
#define  DEM_DEMPACK_2D_CL_H_INCLUDED

/* Project includes */
#include "DEM_Dempack_CL.h"
//#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_Dempack2D : public DEM_Dempack {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Dempack2D);

        DEM_Dempack2D() {}

        //DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

        double mHistoryMaxInd;
        double mHistoryMaxForce;
        double mHistoryDamage;
        double mHistoryDegradation;
        double mHistoryDisp;
        double mHistoryShearFlag;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        ~DEM_Dempack2D() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;
        void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_DEMPACK_2D_H_INCLUDED  defined */
