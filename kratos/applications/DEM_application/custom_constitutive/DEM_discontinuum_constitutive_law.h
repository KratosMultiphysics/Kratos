
#if !defined(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

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

namespace Kratos {
    /**
     * Base class of constitutive laws.
     */

    class SphericParticle; // forward declaration of spheric cont particle

    class /*__declspec( dllexport )*/ DEMDiscontinuumConstitutiveLaw : public Flags {
    public:



        KRATOS_CLASS_POINTER_DEFINITION(DEMDiscontinuumConstitutiveLaw);

        DEMDiscontinuumConstitutiveLaw();

        DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        virtual ~DEMDiscontinuumConstitutiveLaw();

        virtual DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;
        
        
        
        virtual void CalculateContactArea(double mRadius, double other_radius, double &calculation_area);

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area);


        virtual void CalculateForces(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                double kn_el,
                double kt_el,
                double indentation,
                double& failure_criterion_state,
                bool& sliding,
                SphericParticle* element1,
                SphericParticle* element2,
                int &mNeighbourFailureId_count,
                double mapping_new_cont);


        virtual void CalculateNormalForceLinear(double LocalElasticContactForce[3], const double kn_el, const double indentation);
        virtual void CalculateTangentialForceLinear(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                const double kt_el,
                const double indentation,
                double& failure_criterion_state,
                bool& sliding,
                SphericParticle* element1,
                SphericParticle* element2,
                int &mNeighbourFailureId_count,
                double mapping_new_cont);
        virtual void CalculateNormalForceHertz(double LocalElasticContactForce[3], const double kn_el, const double indentation);
        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding,
                int mDampType);

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericParticle* element1,
                SphericParticle* element2,
                double kn_el,
                double kt_el);

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

    KRATOS_DEFINE_VARIABLE(DEMDiscontinuumConstitutiveLaw::Pointer, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

