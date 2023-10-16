/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu, chengshun.shang1996@gmail.com
// Date: June 2023
/////////////////////////////////////////////////

#if !defined(DEM_SMOOTH_JOINT_CL_H_INCLUDE)
#define DEM_SMOOTH_JOINT_CL_H_INCLUDE

// Project includes
#include "DEM_continuum_constitutive_law.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_smooth_joint : public DEMContinuumConstitutiveLaw {

        typedef DEMContinuumConstitutiveLaw BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_smooth_joint);

        DEM_smooth_joint() {}

        void Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps) override;
        void TransferParametersToProperties(const Parameters& parameters, Properties::Pointer pProp) override;
        void Check(Properties::Pointer pProp) const override;

        ~DEM_smooth_joint() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        std::string GetTypeOfLaw() override;

        virtual void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;
        virtual double CalculateContactArea(double radius, double other_radius, Vector& v) override;
        void GetContactArea(const double radius, const double other_radius, const Vector& vector_of_initial_areas, const int neighbour_position, double& calculation_area) override;
        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                    double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) override;
        //virtual void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation);

        double LocalMaxSearchDistance(const int i,
                                    SphericContinuumParticle* element1,
                                    SphericContinuumParticle* element2) override;
        
        //virtual double GetYoungModulusForComputingRotationalMoments(const double& equiv_young);

        virtual void CheckFailure(const int i_neighbour_count, 
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
                const ProcessInfo& r_process_info,
                double& contact_sigma);

        virtual void CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double LocalElasticExtraContactForce[3],
                                            double ViscoDampingLocalContactForce[3],
                                            double LocalCoordSystem[3][3],
                                            double LocalDeltDisp[3],
                                            double LocalRelVel[3],
                                            const double kt_el,
                                            const double equiv_shear,
                                            double& contact_tau,
                                            double indentation,
                                            double calculation_area,
                                            double& failure_criterion_state,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2,
                                            int i_neighbour_count,
                                            bool& sliding,
                                            const ProcessInfo& r_process_info,
                                            int time_steps);

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
                            double LocalElasticContactForce[3],
                            double normalLocalContactForce,
                            double GlobalElasticContactForces[3],
                            double LocalCoordSystem_2[3],
                            const int i_neighbor_count) override;     
        
        void AddContributionOfShearStrainParallelToBond(double OldLocalElasticContactForce[3],
                                                    double LocalElasticExtraContactForce[3],
                                                    array_1d<double, 3>& OldElasticExtraContactForce,
                                                    double LocalCoordSystem[3][3],
                                                    const double kt_el,
                                                    const double calculation_area,
                                                    SphericContinuumParticle* element1,
                                                    SphericContinuumParticle* element2);


        double mAccumulatedJointTangentialLocalDisplacement[2] = {0.0};
        double mKn;
        double mKt;
        double mLocalJointNormal[3] = {0.0};
        double mInitialDistanceJoint = 0.0;

    protected:

    private:
        using DEMContinuumConstitutiveLaw::CalculateNormalForces;
        using DEMContinuumConstitutiveLaw::CalculateTangentialForces;

    };   

} // namespace Kratos

#endif /*DEM_SMOOTH_JOINT_CL_H_INCLUDE defined*/