//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

#include "terminal_velocity_scheme.h"


namespace Kratos {


/// Destructor.
TerminalVelocityScheme::~TerminalVelocityScheme(){}

void TerminalVelocityScheme::UpdateTranslationalVariables(
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
//    if (StepFlag == 1){
//        for (int k = 0; k < 3; k++) {
//            if (Fix_vel[k] == false) {
//                noalias(vel) = i.FastGetSolutionStepValue(VELOCITY);
//                delta_displ[k] = delta_t * vel[k];
//                displ[k] += delta_displ[k];
//                coor[k] = initial_coor[k] + displ[k];
//            }
//        } // dimensions
//    }

    if (StepFlag == 1){
        const array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);

        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
                noalias(vel) = 0.5 * (3 * i.FastGetSolutionStepValue(VELOCITY) - old_vel);
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions
    }
    else {
        const double drag_coefficient_inv = 1.0 / i.FastGetSolutionStepValue(DRAG_COEFFICIENT);
        const array_1d<double, 3 >& fluid_vel = i.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        array_1d<double, 3 > slip_vel = fluid_vel - i.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3 > contact_force =  force - i.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);

        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
    //            KRATOS_WATCH(i.Id())
    //            KRATOS_WATCH(i.FastGetSolutionStepValue(FLUID_DENSITY_PROJECTED));
    //            KRATOS_WATCH(i.FastGetSolutionStepValue(FLUID_VISCOSITY_PROJECTED));

    //            KRATOS_WATCH(coor)
    //            KRATOS_WATCH(vel)
    //            KRATOS_WATCH(contact_force)
                vel[k] = 0.5 * (- vel[k]  + 3 * (fluid_vel[k] + drag_coefficient_inv * contact_force[k])); // all forces are always in equilibrium with the drag force
                //delta_displ[k] = delta_t * vel[k];
//                displ[k] += delta_displ[k];
//                coor[k] = initial_coor[k] + displ[k];

    //            KRATOS_WATCH(drag_coefficient_inv)
    //            KRATOS_WATCH(coor)
    //            KRATOS_WATCH(vel)
            } else {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions
        array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);
        old_vel = vel;
    }
}

void TerminalVelocityScheme::UpdateRotationalVariables(
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
            delta_rotation[k] = angular_velocity[k] * delta_t;
            rotated_angle[k] += delta_rotation[k];
            angular_velocity[k] += delta_t * angular_acceleration[k];
        } else {
            delta_rotation[k] = angular_velocity[k] * delta_t;
            rotated_angle[k] += delta_rotation[k];
        }
    }
}

void TerminalVelocityScheme::CalculateLocalAngularAcceleration(
                            const Node < 3 > & i,
                            const double moment_of_inertia,
                            const array_1d<double, 3 >& torque,
                            const double moment_reduction_factor,
                            array_1d<double, 3 >& angular_acceleration){

    for (int j = 0; j < 3; j++) {
        angular_acceleration[j] = moment_reduction_factor * torque[j] / moment_of_inertia;
    }
}

void TerminalVelocityScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
