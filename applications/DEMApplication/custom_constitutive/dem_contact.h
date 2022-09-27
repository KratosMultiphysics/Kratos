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
            noalias(p_element->mContactMoment) += moment_of_this_neighbour; 
            
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_CONTACT_H_INCLUDED  defined */