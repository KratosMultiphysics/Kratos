/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: chengshun.shang1996@gmail.com
// Date: June 2022
/////////////////////////////////////////////////

#if !defined(DEM_PARALLEL_BOND_CL_H_INCLUDE)
#define DEM_PARALLEL_BOND_CL_H_INCLUDE

// Project includes
#include "DEM_continuum_constitutive_law.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_parallel_bond : public DEMContinuumConstitutiveLaw {

        typedef DEMContinuumConstitutiveLaw BaseClassType;

    public:

        // TODO: what is the function?
        KRATOS_CLASS_POINTER_DEFINITION(DEM_parallel_bond);

        DEM_parallel_bond() {}

        void TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp) override;
        std::string GetTypeOfLaw() override;
        void Check(Properties::Pointer pProp) const override;
        void Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps) override;

        ~DEM_parallel_bond() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        virtual void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;
        virtual double CalculateContactArea(double radius, double other_radius, Vector& v) override;
        void GetContactArea(const double radius, const double other_radius, const Vector& vector_of_initial_areas, const int neighbour_position, double& calculation_area) override;
        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                    double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) override;
        virtual void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation_particle);

        virtual void CalculateUnbondedViscoDampingForce(double LocalRelVel[3],
                                                double UnbondedViscoDampingLocalContactForce[3],
                                                SphericParticle* const element1,
                                                SphericParticle* const element2);

        // TODO: check whether it is necessary 
        double LocalMaxSearchDistance(const int i,
                                    SphericContinuumParticle* element1,
                                    SphericContinuumParticle* element2) override;
        //double GetContactSigmaMax();
        
        //TODO:CHECK
        virtual double GetYoungModulusForComputingRotationalMoments(const double& equiv_young);
        virtual double GetBondKn(double bond_equiv_young, double calculation_area, double distance);

        void CheckFailure(const int i_neighbour_count, 
                            SphericContinuumParticle* element1, 
                            SphericContinuumParticle* element2,
                            double& contact_sigma,
                            double& contact_tau, 
                            double LocalElasticContactForce[3],
                            double ViscoDampingLocalContactForce[3],
                            double ElasticLocalRotationalMoment[3],
                            double ViscoLocalRotationalMoment[3]) override;

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
                            double indentation_particle,
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

        virtual double ComputeNormalUnbondedForce(double unbonded_indentation);

        virtual void CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double indentation_particle,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps,
                const ProcessInfo& r_process_info,
                double& contact_sigma);
            
        void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                                double &equiv_visco_damp_coeff_tangential,
                                SphericContinuumParticle* element1,
                                SphericContinuumParticle* element2,
                                const double kn_el,
                                const double kt_el) override;

        void CalculateViscoDamping(double LocalRelVel[3],
                                double ViscoDampingLocalContactForce[3],
                                double indentation_particle,
                                double equiv_visco_damp_coeff_normal,
                                double equiv_visco_damp_coeff_tangential,
                                bool& sliding,
                                int failure_id,
                                int i_neighbour_count,
                                SphericContinuumParticle* element1,
                                SphericContinuumParticle* element2);

        virtual void CalculateTangentialForces(double OldLocalElasticContactForce[3],
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
                                            double indentation_particle,
                                            double calculation_area,
                                            double& failure_criterion_state,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2,
                                            int i_neighbour_count,
                                            bool& sliding,
                                            const ProcessInfo& r_process_info) override;

        void CalculateMoments(SphericContinuumParticle* element, 
                            SphericContinuumParticle* neighbor, 
                            double equiv_young, 
                            double distance, 
                            double calculation_area,
                            double LocalCoordSystem[3][3], 
                            double ElasticLocalRotationalMoment[3], 
                            double ViscoLocalRotationalMoment[3], 
                            double equiv_poisson, 
                            double indentation, 
                            double indentation_particle,
                            double normalLocalContactForce,
                            double GlobalContactForce[3],
                            double LocalCoordSystem_2[3],
                            const int i_neighbor_count) override;

        void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                SphericContinuumParticle* neighbor,
                                                double equiv_young,
                                                double distance,
                                                double calculation_area,
                                                double LocalCoordSystem[3][3],
                                                double ElasticLocalRotationalMoment[3],
                                                double ViscoLocalRotationalMoment[3],
                                                double equiv_poisson,
                                                double indentation) override; 

        void CalculateBondRotationalDamping(SphericContinuumParticle* element,
                                                SphericContinuumParticle* neighbor,
                                                double LocalCoordSystem[3][3],
                                                double ViscoLocalRotationalMoment[3]);
        
        void AddContributionOfShearStrainParallelToBond(double OldLocalElasticContactForce[3],
                                                    double LocalElasticExtraContactForce[3],
                                                    array_1d<double, 3>& OldElasticExtraContactForce,
                                                    double LocalCoordSystem[3][3],
                                                    const double kt_el,
                                                    const double calculation_area,
                                                    SphericContinuumParticle* element1,
                                                    SphericContinuumParticle* element2);

        double GetTangentialStiffness() override;

        void GetLocalUnbondedContactForce(double LocalUnbondedContactForce[3], double LocalElasticContactForce[3]) override;

        double rand_normal(const double mean, const double stddev);

        double mUnbondedLocalElasticContactForce2 = 0.0;
        double mUnbondedNormalElasticConstant = 0.0;
        double mUnbondedTangentialElasticConstant = 0.0;
        double mUnbondedViscoDampingLocalContactForce[3] = {0.0};
        double mBondedViscoDampingLocalContactForce[3] = {0.0};
        double mBondedLocalContactForce[3] = {0.0};
        double mBondedLocalElasticContactForce2 = 0.0;
        double mBondedScalingFactor[3] = {0.0};
        double mUnBondedScalingFactor[3] = {0.0};
        double mUnbondedEquivViscoDampCoeffTangential = 0.0;
        double mUnbondedEquivViscoDampCoeffNormal = 0.0;
        double mInitialIndentationForBondedPart = 0.0;
        double mAccumulatedBondedTangentialLocalDisplacement[2] = {0.0};
        //double mBondedLocalContactNormalTorque[3] = {0.0};
        //double mBondedLocalContactTangentTorque[3] = {0.0};
        double mKn;
        double mKt;
        bool mDebugPrintingOption;
        double mBondSigmaMax = 0.0;
        double mBondTauZero = 0.0;
        double mCalculation_area = 0.0;
        int mTransientDampingStepsCount = 0;

    protected:

    private:
        using DEMContinuumConstitutiveLaw::CalculateNormalForces;
        using DEMContinuumConstitutiveLaw::CalculateViscoDamping;
        using DEMContinuumConstitutiveLaw::CalculateTangentialForces;

    };   

} // namespace Kratos

#endif /*DEM_PARALLEL_BOND_CL_H_INCLUDE defined*/