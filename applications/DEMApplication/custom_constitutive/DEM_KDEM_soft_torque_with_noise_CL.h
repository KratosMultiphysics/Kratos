
#if !defined(DEM_KDEM_SOFT_TORQUE_WITH_NOISE_H_INCLUDED)
#define DEM_KDEM_SOFT_TORQUE_WITH_NOISE_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_soft_torque_CL.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_soft_torque_with_noise : public DEM_KDEM_soft_torque {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_soft_torque_with_noise);

        DEM_KDEM_soft_torque_with_noise() {
        }

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;
        void Check(Properties::Pointer pProp) const override;
        void Initialize(SphericContinuumParticle* owner_sphere) override;

        ~DEM_KDEM_soft_torque_with_noise() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

    protected:

        double GetTauZero(SphericContinuumParticle* element1) override;

        double GetInternalFricc(SphericContinuumParticle* element1) override;

    private:
        double rand_normal(const double mean, const double stddev);

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override{
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override{
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_KDEM_SOFT_TORQUE_WITH_NOISE_H_INCLUDED  defined */
