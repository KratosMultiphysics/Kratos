//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
// Project includes
#include "terminal_velocity_scheme.h"
#include "swimming_dem_application_variables.h"

namespace Kratos {

    void TerminalVelocityScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node& i,
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

        if (StepFlag == 1){
            const array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);
            const array_1d<double, 3 >& current_vel = i.FastGetSolutionStepValue(VELOCITY);

            //noalias(vel) = old_vel;
            //noalias(vel) = 0.5 * (3 * i.FastGetSolutionStepValue(VELOCITY) - old_vel);
            // trapezoidal
            noalias(vel) = 0.5 * (old_vel + current_vel);
            //true trapezoidal
            //noalias(vel) = current_vel;

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }
        else {
            // const double drag_coefficient_inv = 1.0 / i.FastGetSolutionStepValue(DRAG_COEFFICIENT);
            const array_1d<double, 3 >& fluid_vel = i.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3 > contact_force =  force - i.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);
            array_1d<double, 3 >& force_old = i.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD);

            // hard-coded
            const double rad = i.FastGetSolutionStepValue(RADIUS);
            double disp_volume = (4. * M_PI / 3.) * pow(rad, 3);
            const double g = 9.8;
            const double mu = 1e-3;
            // const double rho_f = 1e3;
            const double rho_f = i.FastGetSolutionStepValue(NODAL_DENSITY_PROJECTED);
            const double rho_p = (mass / disp_volume);  // relative density
            i.FastGetSolutionStepValue(DENSITY) = rho_p;
            const double drag_coeff = 6.0 * Globals::Pi * mu * rad;

            // std::cout << "rho_p = " << rho_p << ", drag_coef = " << drag_coeff << ", disp_volume " << disp_volume << std::endl;
            for (int k = 0; k < 3; k++){
                if (Fix_vel[k] == false){
//                    vel[k] = 0.5 * (- vel[k]  + 3 * (fluid_vel[k] + drag_coefficient_inv * contact_force[k])); // all forces are always in equilibrium with the drag force
                    // adams-bashforth
                    // vel[k] = fluid_vel[k] + 0.5 * drag_coefficient_inv * (3 * contact_force[k] - force_old[k]); // all forces are always in equilibrium with the drag force
                    // vel[k] = fluid_vel[k] + 0.5 * drag_coefficient_inv * (3 * force[k] - force_old[k]); // all forces are always in equilibrium with the drag force
                    // trapezoidal
//                    vel[k] = fluid_vel[k] + 0.5 * drag_coefficient_inv * (contact_force[k] + force_old[k]); // all forces are always in equilibrium with the drag force

                //    delta_displ[k] = delta_t * vel[k];
                //    displ[k] += delta_displ[k];
                //    coor[k] = initial_coor[k] + displ[k];

                // Force drag to be in balance with the forces hardcoded (assume different densities in z direction)
                if ((k == 2) && (rho_f >= 0.0)) {
                    vel[k] = fluid_vel[k] - (disp_volume / drag_coeff) * (rho_p - rho_f) * g;
                } else {
                    vel[k] = fluid_vel[k];
                }

                }
                else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions

            // std::cout << "particle at " << initial_coor << " will have v = " << vel << std::endl;
            // std::cout << "f_contact = " << contact_force << ", f_old = " << force_old << ", fluid_vel = " << fluid_vel << std::endl;

            array_1d<double, 3 >& old_vel = i.FastGetSolutionStepValue(VELOCITY_OLD);
            noalias(old_vel) = vel;
            noalias(force_old) = contact_force;
        }
    }

void TerminalVelocityScheme::UpdateRotationalVariables(
                int StepFlag,
                Node& i,
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
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) {

        double moment_of_inertia_inv = 1.0 / moment_of_inertia;
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] * moment_of_inertia_inv;
        }
    }

void TerminalVelocityScheme::CalculateLocalAngularAccelerationByEulerEquations(
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration) {

        for (int j = 0; j < 3; j++) {
            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
            local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
        }
    }
} //namespace Kratos
