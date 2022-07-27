#if !defined(DEM_CONTACT_H_INCLUDED)
#define  DEM_CONTACT_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "custom_elements/spheric_particle.h"


namespace Kratos {

    class DemContact {

        //class Properties; //forward declaration
        //class SphericContinuumParticle; // forward declaration of spheric cont particle
        //class SphericParticle;
        //class DEMFlags;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DemContact);

        DemContact();

        virtual ~DemContact();

        static void ComputeParticleContactMoments(double NormalLocalContactForce,
                                                double Force[3],
                                                double& RollingResistance,
                                                double LocalCoordSystem2[3],
                                                SphericParticle* p_element,
                                                SphericParticle* p_neighbour,
                                                double indentation,
                                                unsigned int i)
        {
            double arm_length = p_element->GetInteractionRadius() - indentation;

            const double other_young = p_neighbour->GetYoung();
            arm_length = p_element->GetInteractionRadius() - indentation * other_young / (other_young + p_element->GetYoung());

            array_1d<double, 3> arm_vector;
            arm_vector[0] = -LocalCoordSystem2[0] * arm_length;
            arm_vector[1] = -LocalCoordSystem2[1] * arm_length;
            arm_vector[2] = -LocalCoordSystem2[2] * arm_length;

            array_1d<double, 3> moment_of_this_neighbour;
            GeometryFunctions::CrossProduct(arm_vector, Force, moment_of_this_neighbour);
            noalias(p_element->mContactMoment) += moment_of_this_neighbour; //TODO: RIGHT OR NOT?

            // ROLLING FRICTION
            if (p_element->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
                Properties& properties_of_this_contact = p_element->GetProperties().GetSubProperties(p_neighbour->GetProperties().Id());
                const double min_radius = std::min(p_element->GetRadius(), p_neighbour->GetRadius());
                const double equiv_rolling_friction_coeff = properties_of_this_contact[ROLLING_FRICTION] * min_radius;

                if (equiv_rolling_friction_coeff) {
                    p_element->ComputeRollingResistance(RollingResistance, NormalLocalContactForce, equiv_rolling_friction_coeff, i);
                }
            }
            
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_CONTACT_H_INCLUDED  defined */