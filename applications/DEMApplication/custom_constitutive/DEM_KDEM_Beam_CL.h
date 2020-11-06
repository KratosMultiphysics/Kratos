
#if !defined(DEM_KDEM_Beam_CL_H_INCLUDED)
#define  DEM_KDEM_Beam_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
#include "DEM_KDEM_CL.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_Beam : public DEM_KDEM {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_Beam);

        DEM_KDEM_Beam() {}

        ~DEM_KDEM_Beam() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;
        void Check(Properties::Pointer pProp) const override;

        virtual void CalculateElasticConstants(double& kn_el,
                                               double& kt_el_0,
                                               double& kt_el_1,
                                               double initial_dist,
                                               double equiv_young,
                                               double equiv_poisson,
                                               double calculation_area,
                                               SphericContinuumParticle* element1,
                                               SphericContinuumParticle* element2) override;

        virtual void CalculateViscoDampingCoeff(double& equiv_visco_damp_coeff_normal,
                                                double& equiv_visco_damp_coeff_tangential_0,
                                                double& equiv_visco_damp_coeff_tangential_1,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                const double kn_el,
                                                const double kt_el_0,
                                                const double kt_el_1);

        virtual void CalculateForces(const ProcessInfo& r_process_info,
                                     double OldLocalElasticContactForce[3],
                                     double LocalElasticContactForce[3],
                                     double LocalElasticExtraContactForce[3],
                                     double LocalCoordSystem[3][3],
                                     double LocalDeltDisp[3],
                                     const double kn_el,
                                     const double kt_el_0,
                                     const double kt_el_1,
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
                                     double& equiv_visco_damp_coeff_normal,
                                     double& equiv_visco_damp_coeff_tangential_0,
                                     double& equiv_visco_damp_coeff_tangential_1,
                                     double LocalRelVel[3],
                                     double ViscoDampingLocalContactForce[3]) override;

        virtual void CalculateNormalForces(double LocalElasticContactForce[3],
                                           const double kn_el,
                                           double indentation);

        virtual void CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                               double LocalElasticContactForce[3],
                                               double LocalDeltDisp[3],
                                               const double kt_el_0,
                                               const double kt_el_1);

        virtual void CalculateViscoDamping(double LocalRelVel[3],
                                           double ViscoDampingLocalContactForce[3],
                                           double equiv_visco_damp_coeff_normal,
                                           double equiv_visco_damp_coeff_tangential_0,
                                           double equiv_visco_damp_coeff_tangential_1);

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
#endif /* DEM_KDEM_Beam_H_INCLUDED  defined */
