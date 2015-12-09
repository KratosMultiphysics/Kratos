

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
            array_1d<double, 3 >& rotated_angle,
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
                rotated_angle[k] += delta_rotation[k];
                angular_velocity[k] += delta_t * moment_reduction_factor * torque[k] / moment_of_inertia;
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }    
    
    void ForwardEulerScheme::UpdateRotationalVariablesOfClusters(
        const Node < 3 > & i,
        array_1d<double, 3 >& rotated_angle,
        array_1d<double, 3 >& delta_rotation,
        array_1d<double, 3 >& angular_velocity,
        const array_1d<double, 3 >& angular_acceleration,
        const double delta_t,
        const bool Fix_Ang_vel[3]) {

        for (int j = 0; j < 3; j++) {
            if (Fix_Ang_vel[j] == false) {
                delta_rotation[j] = angular_velocity[j] * delta_t;
                angular_velocity[j] += angular_acceleration[j] * delta_t;
                rotated_angle[j] += delta_rotation[j];
            } else {
                angular_velocity[j] = 0.0;
                delta_rotation[j] = 0.0;
            }
        }
    }
    
    void ForwardEulerScheme::CalculateLocalAngularAccelerationByEulerEquations(
                                const Node < 3 > & i,
                                const array_1d<double, 3 >& local_angular_velocity,
                                const array_1d<double, 3 >& moments_of_inertia,
                                const array_1d<double, 3 >& local_torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& local_angular_acceleration){
        for (int j = 0; j < 3; j++) {
            //Euler equations in Explicit (Forward Euler) scheme:
            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
            local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
            
        }
    }
} //namespace Kratos
