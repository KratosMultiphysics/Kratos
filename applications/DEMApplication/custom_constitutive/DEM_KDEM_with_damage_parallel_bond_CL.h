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

        void TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp) override;

        void Check(Properties::Pointer pProp) const override;

        void Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps) override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double bonded_equiv_young,
                                        double equiv_poisson, double calculation_area,
                                        SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) override;

        double GetYoungModulusForComputingRotationalMoments(const double& equiv_young) override;

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
                                double &equiv_visco_damp_coeff_normal,
                                double &equiv_visco_damp_coeff_tangential,
                                double LocalRelVel[3],
                                double ViscoDampingLocalContactForce[3]) override;

        void CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                              double& equiv_visco_damp_coeff_tangential,
                                              SphericContinuumParticle* element1,
                                              SphericContinuumParticle* element2,
                                              const double kn_el,
                                              const double kt_el) override;

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
            double ViscoDampingLocalContactForce[3],
            double LocalCoordSystem[3][3],
            double LocalDeltDisp[3],
            double LocalRelVel[3],
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
            const ProcessInfo& r_process_info) override;

        virtual void ComputeNormalUnbondedForce(double indentation);

        virtual double LocalMaxSearchDistance(const int i, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

        virtual void CalculateNormalAndTangentialDamageComponents(SphericContinuumParticle* element1, SphericContinuumParticle* element2);

        double mUnbondedLocalElasticContactForce2 = 0.0;
        double mUnbondedNormalElasticConstant = 0.0;
        double mUnbondedTangentialElasticConstant = 0.0;
        double mUnbondedViscoDampingLocalContactForce[3] = {0.0};
        double mBondedViscoDampingLocalContactForce[3] = {0.0};
        double mBondedScalingFactor = 0.0;
        double mUnbondedScalingFactor = 0.0;
        bool mDebugPrintingOption = false;
        double mDamageEnergyCoeff = 0.0;
        double mUnbondedEquivViscoDampCoeffTangential = 0.0;
        double mUnbondedEquivViscoDampCoeffNormal = 0.0;
        double mInitialIndentationForBondedPart = 0.0;
        double mAccumulatedBondedTangentialLocalDisplacement[2] = {0.0};

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
