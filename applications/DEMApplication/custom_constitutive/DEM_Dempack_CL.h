
#if !defined(DEM_DEMPACK_CL_H_INCLUDED)
#define  DEM_DEMPACK_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
//#include "DEM_discontinuum_constitutive_law.h"



namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_Dempack : public DEMContinuumConstitutiveLaw {
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


        virtual void Initialize(SphericContinuumParticle* owner_sphere) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        ~DEM_Dempack() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        virtual void CalculateContactArea(double radius,
                                          double other_radius,
                                          double &calculation_area) override;

        virtual void GetContactArea(const double radius,
                                    const double other_radius,
                                    const Vector& vector_of_initial_areas,
                                    const int neighbour_position,
                                    double& calculation_area) override;

        virtual void CalculateElasticConstants(double &kn_el,
                                               double &kt_el,
                                               double initial_dist,
                                               double equiv_young,
                                               double equiv_poisson,
                                               double calculation_area,
                                               SphericContinuumParticle* element1,
                                               SphericContinuumParticle* element2) override;

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                const double kn_el,
                const double kt_el) override;

        virtual double LocalMaxSearchDistance(const int i,
                                              SphericContinuumParticle* element1,
                                              SphericContinuumParticle* element2) override;

        virtual void CalculateForces(const ProcessInfo& r_process_info,
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

        virtual void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                      SphericContinuumParticle* neighbor,
                                                      double equiv_young,
                                                      double distance,
                                                      double calculation_area,
                                                      double LocalCoordSystem[3][3],
                                                      double ElasticLocalRotationalMoment[3],
                                                      double ViscoLocalRotationalMoment[3],
                                                      double equiv_poisson,
                                                      double indentation) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override{
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override{
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_DEMPACK_H_INCLUDED  defined */
