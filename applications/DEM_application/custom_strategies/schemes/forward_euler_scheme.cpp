

#include "DEM_application.h"
#include "forward_euler_scheme.h"

namespace Kratos {

    void ForwardEulerScheme::UpdateTranslationalVariables(
            const ModelPart::NodeIterator& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
                vel[k] += delta_t * force_reduction_factor * force[k] / mass;
            } else {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions  
    }

    void ForwardEulerScheme::UpdateRotationalVariables(
            const ModelPart::NodeIterator& i,
            array_1d<double, 3 >& rotational_displacement,
            array_1d<double, 3 >& delta_rotation,
            array_1d<double, 3 >& angular_velocity,
            const array_1d<double, 3 >& torque,
            const double moment_reduction_factor,
            const double moment_of_inertia,
            const double delta_t,
            const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotational_displacement[k] += delta_rotation[k];
                angular_velocity[k] += delta_t * moment_reduction_factor * torque[k] / moment_of_inertia;
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotational_displacement[k] += delta_rotation[k];
            }
        }

    }    
} //namespace Kratos
