
#if !defined(DEM_DEMPACK_2D_CL_H_INCLUDED)
#define  DEM_DEMPACK_2D_CL_H_INCLUDED

/* Project includes */
#include "DEM_Dempack_CL.h"
//#include "DEM_discontinuum_constitutive_law.h"



namespace Kratos {

    class DEM_Dempack2D : public DEM_Dempack {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Dempack2D);

        DEM_Dempack2D() {
        }

        //DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

        double mHistoryMaxInd;
        double mHistoryMaxForce;
        double mHistoryDamage;
        double mHistoryDegradation;
        double mHistoryDisp;
        double mHistoryShearFlag;


        void Initialize(const ProcessInfo& rCurrentProcessInfo);

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_Dempack2D() {
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

        void CalculateForces(ProcessInfo& rCurrentProcessInfo,
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
                int &mNeighbourFailureId_count,
                int &mIniNeighbourFailureId_mapping,
                double &mNeighbourDelta_count,
                int time_steps,
                bool& sliding,
                int search_control,
                vector<int>& search_control_vector,
                double mapping_new_cont);

        void CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int &mNeighbourFailureId_count,
                int &mIniNeighbourFailureId_mapping,
                double &mNeighbourDelta_count,
                int time_steps);


        void CalculateTangentialForces(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                double kt_el,
                double indentation,
                double calculation_area,
                double& failure_criterion_state,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int &mNeighbourFailureId_count,
                int &mIniNeighbourFailureId_mapping,
                bool& sliding,
                int search_control,
                vector<int>& search_control_vector,
                double mapping_new_cont);


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
#endif /* DEM_DEMPACK_2D_H_INCLUDED  defined */
