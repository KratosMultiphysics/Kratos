
#if !defined(DEM_EXPONENTIALHC_CL_H_INCLUDED)
#define  DEM_EXPONENTIALHC_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"


namespace Kratos {

    class DEM_ExponentialHC : public DEMContinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_ExponentialHC);

        DEM_ExponentialHC() {
        }


        //    virtual void NonlinearNormalForceCalculation(double LocalElasticContactForce[3], double kn1, double kn2, double distance, double max_dist, double initial_dist);            

        double mHistoryMaxInd;
        double mHistoryMaxForce;
        double mHistoryDamage;
        double mHistoryDegradation;
        double mHistoryDisp;
        double mHistoryShearFlag;
        double mGamma1;
        double mGamma2;
        double mGamma3;
        double mMaxDef;

        void Initialize(const ProcessInfo& r_process_info);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_ExponentialHC() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const;

        void CalculateContactArea(double radius, double other_radius, double &calculation_area);
        void CalculateElasticConstants(double &kn_el, double &kt_el, double initial_dist, double equiv_young, double equiv_poisson, double calculation_area);

        void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                double kn_el,
                double kt_el);

        void CalculateForces(ProcessInfo& r_process_info,
                             double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                const double kn_el,
                double kt_el,
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
                double indentation,
                double calculation_area,
                double& failure_criterion_state,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                bool& sliding,
                int search_control,
                vector<int>& search_control_vector);


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

