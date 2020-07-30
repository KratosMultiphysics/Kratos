#if !defined(DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_H_INCLUDED)
#define DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_with_damage_CL.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_with_damage_parallel_bond : public DEM_KDEM_with_damage {

        typedef DEM_KDEM_with_damage BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_with_damage_parallel_bond);

        DEM_KDEM_with_damage_parallel_bond() {}

        ~DEM_KDEM_with_damage_parallel_bond() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        void SetDebugPrintingOptionValue(Properties::Pointer pProp);

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double bonded_equiv_young,
                                             double equiv_poisson, double calculation_area,
                                             SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

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

        void CalculateViscoDamping(double LocalRelVel[3],
                                         double ViscoDampingLocalContactForce[3],
                                         double indentation,
                                         double equiv_visco_damp_coeff_normal,
                                         double equiv_visco_damp_coeff_tangential,
                                         bool& sliding,
                                         int failure_id) override;

        void CalculateNormalForces(double LocalElasticContactForce[3],
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

        void CalculateTangentialForces(double OldLocalElasticContactForce[3],
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

        double LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

        void AdjustEquivalentYoung(double& equiv_young, const SphericContinuumParticle* element, const SphericContinuumParticle* neighbor) override;

        virtual void AdjustTauStrengthAndUpdatedMaxTauStrength(double& tau_strength, double& updated_max_tau_strength, const double internal_friction,
                                                               double contact_sigma, SphericContinuumParticle* element1, SphericContinuumParticle* element2);

        double mUnbondedLocalElasticContactForce2 = 0.0;
        double mUnbondedNormalElasticConstant;
        double mUnbondedTangentialElasticConstant;
        double mViscoDampingLocalContactForce[3];
        double mBondedScalingFactor = 0.0;
        double mUnbondedScalingFactor = 0.0;
        bool mDebugPrintingOption;
        double mDamageEnergyCoeff;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseClassType)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseClassType)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} // namespace Kratos
#endif /* DEM_KDEM_WITH_DAMAGE_PARALLEL_BOND_H_INCLUDED defined */
