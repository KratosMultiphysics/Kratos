
#if !defined(DEM_DEMPACK_CL_H_INCLUDED)
#define  DEM_DEMPACK_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
//#include "DEM_discontinuum_constitutive_law.h"



namespace Kratos {

    class DEM_Dempack : public DEMContinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Dempack);

        DEM_Dempack() {
        }

        //DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

        double mHistoryMaxInd;
        double mHistoryMaxForce;
        double mHistoryDamage;
        double mHistoryDegradation;
        double mHistoryDisp;
        double mHistoryShearFlag;


        virtual void Initialize();

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_Dempack() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactArea(double radius, double other_radius, double &calculation_area);
        virtual void CalculateElasticConstants(double &kn_el, double &kt_el, double initial_dist, double equiv_young, double equiv_poisson, double calculation_area);

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                double kn_el,
                double kt_el);

        virtual void CalculateForces(const ProcessInfo& r_process_info,
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

        virtual void CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps);


        virtual void CalculateTangentialForces(double LocalElasticContactForce[3],
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
        
        virtual void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                      SphericContinuumParticle* neighbor,
                                                      double equiv_young,
                                                      double distance,
                                                      double calculation_area,
                                                      double LocalCoordSystem[3][3],
                                                      double ElasticLocalRotationalMoment[3],
                                                      double ViscoLocalRotationalMoment[3]);

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
#endif /* DEM_DEMPACK_H_INCLUDED  defined */
