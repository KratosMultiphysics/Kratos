
#if !defined(DEM_EXPONENTIALHC_CL_H_INCLUDED)
#define  DEM_EXPONENTIALHC_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"


namespace Kratos {

    class DEM_ExponentialHC : public DEMContinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_ExponentialHC);

        DEM_ExponentialHC() {}

        double mHistoryMaxInd;
        double mHistoryMaxForce;
        double mHistoryDamage;
        double mHistoryDegradation;
        double mGamma1;
        double mGamma2;
        double mGamma3;
        double mMaxDef;

        void Initialize();

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_ExponentialHC() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const;

        void CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps);

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

} /* namespace Kratos.*/
#endif /* DEM_EXPONENTIALHC_H_INCLUDED  defined */

