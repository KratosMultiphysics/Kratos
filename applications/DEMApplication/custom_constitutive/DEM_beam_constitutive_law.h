
#if !defined(DEM_BEAM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_BEAM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/discrete_element.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "containers/array_1d.h"


namespace Kratos {

    class Properties; //forward declaration
    class SphericContinuumParticle; // forward declaration of spheric cont particle

    class KRATOS_API(DEM_APPLICATION) DEMBeamConstitutiveLaw : public Flags {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMBeamConstitutiveLaw);

        DEMBeamConstitutiveLaw();

        virtual void Initialize(SphericContinuumParticle* element1, SphericContinuumParticle* element2, Properties::Pointer pProps);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true);

        virtual void SetConstitutiveLawInPropertiesWithParameters(Properties::Pointer pProp, const Parameters& parameters, bool verbose = true);

        virtual void Check(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual ~DEMBeamConstitutiveLaw();

        virtual DEMBeamConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateElasticConstants(double& kn_el,
                                               double& kt_el_0,
                                               double& kt_el_1,
                                               double initial_dist,
                                               double equiv_young,
                                               double equiv_poisson,
                                               double calculation_area,
                                               SphericContinuumParticle* element1,
                                               SphericContinuumParticle* element2, double indentation);

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
                                     double ViscoDampingLocalContactForce[3]);

        virtual void CalculateNormalForces(double LocalElasticContactForce[3],
                                           const double kn_el,
                                           double indentation);

        virtual void CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                               double LocalElasticContactForce[3],
                                               double LocalDeltDisp[3],
                                               double LocalRelVel[3],
                                               const double kt_el_0,
                                               const double kt_el_1);

        virtual void CalculateViscoDamping(double LocalRelVel[3],
                                           double ViscoDampingLocalContactForce[3],
                                           double equiv_visco_damp_coeff_normal,
                                           double equiv_visco_damp_coeff_tangential_0,
                                           double equiv_visco_damp_coeff_tangential_1);
                                
        virtual void CalculateMoments(SphericContinuumParticle* element, 
                                      SphericContinuumParticle* neighbor, 
                                      double equiv_young, 
                                      double distance, 
                                      double calculation_area,
                                      double LocalCoordSystem[3][3], 
                                      double ElasticLocalRotationalMoment[3], 
                                      double ViscoLocalRotationalMoment[3], 
                                      double equiv_poisson, 
                                      double indentation,
                                      double normalLocalContactForce,
                                      double GlobalElasticContactForces[3],
                                      double LocalCoordSystem_2[3],
                                      const int i_neighbor_count);

        virtual void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                      SphericContinuumParticle* neighbor,
                                                      double equiv_young,
                                                      double distance,
                                                      double calculation_area,
                                                      double LocalCoordSystem[3][3],
                                                      double ElasticLocalRotationalMoment[3],
                                                      double ViscoLocalRotationalMoment[3],
                                                      double equiv_poisson,
                                                      double indentation);

        virtual bool CheckRequirementsOfStressTensor();

    protected:

        Properties::Pointer mpProperties;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    //This definition is done here to avoid recursive inclusion of header files
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMBeamConstitutiveLaw::Pointer, DEM_BEAM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_BEAM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

