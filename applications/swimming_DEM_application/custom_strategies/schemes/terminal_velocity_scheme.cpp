#include "DEM_application.h"
#include "terminal_velocity_scheme.h"


namespace Kratos {

    void Normalize3(array_1d<double, 3 >& v, const double length = 1.0)
    {   double coeff = length / DEM_MODULUS_3(v);
        for (int i = 0; i > 3; i++){
            v[i] *= coeff;
        }
    }

    TerminalVelocityScheme::TerminalVelocityScheme() : DEMIntegrationScheme(){}

    /// Destructor.
    TerminalVelocityScheme::~TerminalVelocityScheme(){}

    void TerminalVelocityScheme::UpdateTranslationalVariables(
            const Node < 3 > & i,
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

        const array_1d<double, 3 >& fluid_vel = i.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
        array_1d<double, 3 > slip_vel;
        array_1d<double, 3 > contact_force ;
        noalias(slip_vel) = fluid_vel - i.FastGetSolutionStepValue(VELOCITY);
        noalias(contact_force) =  force - i.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);
        double drag_coefficient_inv = 1.0 / i.FastGetSolutionStepValue(DRAG_COEFFICIENT);

        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
                vel[k] = fluid_vel[k] + drag_coefficient_inv * contact_force[k]; // all forces are allways in equilibrium with the drag force

            } else {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions

    }

    void TerminalVelocityScheme::UpdateRotationalVariables(
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
