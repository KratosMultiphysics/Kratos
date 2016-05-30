//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
#include "includes/dem_variables.h"
#include "hybrid_bashforth_scheme.h"

namespace Kratos {

    void HybridBashforthScheme::AddSpheresVariables(ModelPart & r_model_part){
        DEMIntegrationScheme::AddSpheresVariables(r_model_part);
    }

    void HybridBashforthScheme::AddClustersVariables(ModelPart & r_model_part){
        DEMIntegrationScheme::AddClustersVariables(r_model_part);
    }

    void HybridBashforthScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 > & i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            const array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3])
    {   
        array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);

        if (StepFlag == 1){
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = 0.5 * delta_t * (3 * vel[k] - old_vel[k]);
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }

        else {
            noalias(old_vel) = vel;

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    vel[k] += delta_t * force_reduction_factor * force[k] / mass;
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }
    }

    void HybridBashforthScheme::UpdateRotationalVariables(
                int StepFlag,
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                angular_velocity[k] += delta_t * angular_acceleration[k];
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }

    void HybridBashforthScheme::CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque,
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration){

        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] / moment_of_inertia;
        }
    }


    void HybridBashforthScheme::CalculateLocalAngularAccelerationByEulerEquations(
                                const Node < 3 > & i,
                                const array_1d<double, 3 >& local_angular_velocity,
                                const array_1d<double, 3 >& moments_of_inertia,
                                const array_1d<double, 3 >& local_torque,
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& local_angular_acceleration){

        for (int j = 0; j < 3; j++) {
            //Euler equations in Explicit (Symplectic Euler) scheme:
            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
            local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
        }
    }
} //namespace Kratos
