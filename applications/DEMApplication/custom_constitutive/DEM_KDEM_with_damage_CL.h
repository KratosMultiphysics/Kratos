#if !defined(DEM_KDEM_WITH_DAMAGE_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_soft_torque_CL.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage : public DEM_KDEM_soft_torque {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage);

        DEM_KDEM_with_damage() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        ~DEM_KDEM_with_damage() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void Check(Properties::Pointer pProp) const override;

        void Initialize(SphericContinuumParticle* element) override;

        virtual void CalculateTangentialForces(double OldLocalElasticContactForce[3],
            double LocalElasticContactForce[3],
            double LocalElasticExtraContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            const double kt_el,
            const double equiv_shear,
            double& contact_sigma,
            double& contact_tau,
            double indentation,
            double calculation_area,
            double& failure_criterion_state,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            bool& sliding,
            int search_control,
            DenseVector<int>& search_control_vector,
            const ProcessInfo& r_process_info) override;

        double mKtUpdated = 1e30;
        double mDamage = 0.0;

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

} // namespace Kratos
#endif /* DEM_KDEM_WITH_DAMAGE_H_INCLUDED  defined */
