
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"
//#include "../custom_elements/spheric_continuum_particle.h"


namespace Kratos {

    /**
     * Base class of constitutive laws.
     */

    class SphericContinuumParticle; // forward declaration of spheric cont particle

    class /*__declspec( dllexport )*/ DEMContinuumConstitutiveLaw : public Flags {
    
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMContinuumConstitutiveLaw);

        DEMContinuumConstitutiveLaw();

        DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw& rReferenceContinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
        
        virtual std::string GetTypeOfLaw();

        virtual ~DEMContinuumConstitutiveLaw();

        virtual DEMContinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding);

        virtual void CalculateContactArea(double radius,
                double other_radius,
                double &calculation_area) {
        };

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area) {
        };

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                double kn_el,
                double kt_el) {
        };

        virtual void CalculateForces(ProcessInfo& rCurrentProcessInfo,
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
                vector<int>& search_control_vector) {
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
                int time_steps) {
            std::cout << "error this base class should be overloaded" << std::endl;
        }

        virtual void CalculateTangentialForces(double LocalElasticContactForce[3],
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
                vector<int>& search_control_vector) {
        };

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    KRATOS_DEFINE_VARIABLE(DEMContinuumConstitutiveLaw::Pointer, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

