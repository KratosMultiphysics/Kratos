#if !defined(DEM_KDEM_WITH_DAMAGE_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_soft_torque_CL.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage : public DEM_KDEM_soft_torque {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage);

        DEM_KDEM_with_damage() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        ~DEM_KDEM_with_damage() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void Check(Properties::Pointer pProp) const override;

        void Initialize(SphericContinuumParticle* element) override;

        void CalculateForces(const ProcessInfo& r_process_info,
                             double OldLocalElasticContactForce[3],
                             double LocalElasticContactForce[3],
                             double LocalElasticExtraContactForce[3],
                             double LocalCoordSystem[3][3],
                             double LocalDeltDisp[3],
                             const double kn_el,
                             const double kt_el,
                             double& contact_sigma,
                             double& contact_tau,
                             double& failure_criterion_state,
                             double equiv_young,
                             double equiv_shear,
                             double indentation,
                             double calculation_area,
                             double& acumulated_damage,
                             SphericContinuumParticle* element1,
                             SphericContinuumParticle* element2,
                             int i_neighbour_count,
                             int time_steps,
                             bool& sliding,
                             int search_control,
                             DenseVector<int>& search_control_vector,
                             double &equiv_visco_damp_coeff_normal,
                             double &equiv_visco_damp_coeff_tangential,
                             double LocalRelVel[3],
                             double ViscoDampingLocalContactForce[3]) override;

        virtual void CalculateNormalForces(double LocalElasticContactForce[3],
            const double kn_el,
            double equiv_young,
            double indentation,
            double calculation_area,
            double& acumulated_damage,
            SphericContinuumParticle* element1,
            SphericContinuumParticle* element2,
            int i_neighbour_count,
            int time_steps,
            const ProcessInfo& r_process_info) override;

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

        void FindMaximumValueOfNormalAndTangentialDamageComponents();

        double mDamageNormal = 0.0;
        double mDamageTangential = 0.0;
        const double mDamageThresholdTolerance = 0.99;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEM_KDEM_soft_torque)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEM_KDEM_soft_torque)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} // namespace Kratos
#endif /* DEM_KDEM_WITH_DAMAGE_H_INCLUDED  defined */
