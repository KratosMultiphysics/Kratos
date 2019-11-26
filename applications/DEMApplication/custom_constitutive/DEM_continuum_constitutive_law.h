
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
//#include "includes/properties.h"


#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"


namespace Kratos {

    class Properties; //forward declaration
    class SphericContinuumParticle; // forward declaration of spheric cont particle

    class KRATOS_API(DEM_APPLICATION) DEMContinuumConstitutiveLaw : public Flags {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMContinuumConstitutiveLaw);

        DEMContinuumConstitutiveLaw();

        DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw& rReferenceContinuumConstitutiveLaw);

        virtual void Initialize(SphericContinuumParticle* owner_sphere);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true);

        virtual void Check(Properties::Pointer pProp) const;

        virtual std::string GetTypeOfLaw();

        virtual ~DEMContinuumConstitutiveLaw();

        virtual DEMContinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool& sliding,
                int failure_id);

        virtual void CalculateContactArea(double radius,
                double other_radius,
                double& calculation_area) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateContactArea) should not be called.","")
        };

        virtual double CalculateContactArea(double radius,
                double other_radius,
                Vector& v) {
            return 0.0;
        }

        virtual void GetContactArea(const double radius,
                                    const double other_radius,
                                    const Vector& vector_of_initial_areas,
                                    const int neighbour_position,
                                    double& calculation_area){

            CalculateContactArea(radius, other_radius, calculation_area);

        }

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateElasticConstants) should not be called.","")
        };

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                                                double &equiv_visco_damp_coeff_tangential,
                                                SphericContinuumParticle* element1,
                                                SphericContinuumParticle* element2,
                                                const double kn_el,
                                                const double kt_el) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateViscoDampingCoeff) should not be called.","")
        };

        virtual void CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2) {

        }

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
                double ViscoDampingLocalContactForce[3]) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateForces) should not be called.","")
        };

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
            const ProcessInfo& r_process_info) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateNormalForces) should not be called.","")
        }

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
                                            const ProcessInfo& r_process_info) {
            KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::CalculateTangentialForces) should not be called.","")
        };

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

        virtual void AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, double calculation_area, BoundedMatrix<double, 3, 3>* mSymmStressTensor,
                                            SphericContinuumParticle* element1, SphericContinuumParticle* element2, const ProcessInfo& r_process_info, const int i_neighbor_count, const double indentation);

        virtual double LocalMaxSearchDistance(const int i,
                                              SphericContinuumParticle* element1,
                                              SphericContinuumParticle* element2);

        virtual double LocalPeriod(const int i,
                                   SphericContinuumParticle* element1,
                                   SphericContinuumParticle* element2);

        virtual bool CheckRequirementsOfStressTensor();

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
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMContinuumConstitutiveLaw::Pointer, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

