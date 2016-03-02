
#if !defined(DEM_KDEM_CL_H_INCLUDED)
#define  DEM_KDEM_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"


namespace Kratos {

    class DEM_KDEM : public DEMContinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM);

        DEM_KDEM() {
        }    

        void Initialize(const ProcessInfo& r_process_info);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_KDEM() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const;

        void CalculateContactArea(double radius, double other_radius, double& calculation_area);
        double CalculateContactArea(double radius, double other_radius, std::vector<double>& v);
        void GetContactArea(const double radius, const double other_radius, const std::vector<double> & vector_of_initial_areas, const int neighbour_position, double& calculation_area);
        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young, double equiv_poisson, double calculation_area);

        void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                double kn_el,
                double kt_el);

        void CalculateForces(const ProcessInfo& r_process_info,
                             double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                const double kn_el,
                double kt_el,
                double& contact_sigma,
                double& contact_tau,
                double& failure_criterion_state,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps,
                bool& sliding,
                int search_control,
                vector<int>& search_control_vector);


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


        void CalculateTangentialForces(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                double kt_el,
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
                vector<int>& search_control_vector);
        
        
        void CalculateViscoDamping(double LocalRelVel[3],
                                   double ViscoDampingLocalContactForce[3],
                                   double indentation,
                                   double equiv_visco_damp_coeff_normal,
                                   double equiv_visco_damp_coeff_tangential,
                                   bool sliding,
                                   int failure_id);
        

        virtual void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                              SphericContinuumParticle* neighbor,
                                              double equiv_young,
                                              double distance,
                                              double calculation_area,
                                              double LocalCoordSystem[3][3],
                                              array_1d<double, 3>& mContactMoment);
        
        void AddPoissonContribution(const double equiv_poisson, 
                                    double LocalCoordSystem[3][3], 
                                    double& normal_force, 
                                    double calculation_area, Matrix* mSymmStressTensor);

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
#endif /* DEM_KDEM_H_INCLUDED  defined */
