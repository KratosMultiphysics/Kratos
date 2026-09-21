//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//
// Project includes
#include <cmath>

#include "terminal_velocity_scheme.h"
#include "swimming_dem_application_variables.h"

namespace Kratos {

    TerminalVelocityScheme::TerminalVelocityScheme(Parameters rParameters)
        : mDynamicViscosity(0.0), mGravity(ZeroVector(3)), mIsConfigured(false)
    {
        KRATOS_ERROR_IF_NOT(rParameters.Has("dynamic_viscosity"))
            << "TerminalVelocityScheme: the parameters must contain \"dynamic_viscosity\" (no hidden default is used). "
            << "Received:\n" << rParameters.PrettyPrintJsonString() << std::endl;
        KRATOS_ERROR_IF_NOT(rParameters.Has("gravity"))
            << "TerminalVelocityScheme: the parameters must contain the gravity vector \"gravity\" (no hidden default is used). "
            << "Received:\n" << rParameters.PrettyPrintJsonString() << std::endl;
        rParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        SetDynamicViscosity(rParameters["dynamic_viscosity"].GetDouble());
        const Vector g = rParameters["gravity"].GetVector();
        KRATOS_ERROR_IF(g.size() != 3) << "TerminalVelocityScheme: \"gravity\" must have 3 components." << std::endl;
        array_1d<double, 3> gravity;
        for (unsigned k = 0; k < 3; ++k) gravity[k] = g[k];
        SetGravity(gravity);
    }

    void TerminalVelocityScheme::SetDynamicViscosity(const double Viscosity)
    {
        KRATOS_ERROR_IF_NOT(Viscosity > 0.0 && std::isfinite(Viscosity))
            << "TerminalVelocityScheme: the dynamic viscosity must be finite and positive, got " << Viscosity << std::endl;
        mDynamicViscosity = Viscosity;
        mIsConfigured = true;
    }

    void TerminalVelocityScheme::SetGravity(const array_1d<double, 3>& rGravity)
    {
        for (unsigned k = 0; k < 3; ++k)
            KRATOS_ERROR_IF_NOT(std::isfinite(rGravity[k])) << "TerminalVelocityScheme: non-finite gravity component." << std::endl;
        noalias(mGravity) = rGravity;
        mIsConfigured = mIsConfigured && mDynamicViscosity > 0.0;
    }

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

            // trapezoidal
            noalias(vel) = 0.5 * (old_vel + current_vel);

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions
        }
        else {
            KRATOS_ERROR_IF_NOT(mIsConfigured)
                << "TerminalVelocityScheme: the fluid dynamic viscosity and the gravity were not set. "
                << "Construct the scheme with Parameters {\"dynamic_viscosity\": mu, \"gravity\": [gx, gy, gz]} "
                << "(SwimmingDEM: custom_dem.terminal_velocity_scheme_parameters in the ProjectParameters)." << std::endl;

            const array_1d<double, 3 >& fluid_vel = i.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            const array_1d<double, 3 > contact_force =  force - i.FastGetSolutionStepValue(HYDRODYNAMIC_FORCE);
            array_1d<double, 3 >& force_old = i.FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD);

            // Terminal (Stokes) settling velocity of an inertia-free sphere:
            //     v = u_f + (V_p / (6 pi mu a)) (rho_p - rho_f) g   with  V_p / (6 pi mu a) = 2 a^2 / (9 mu).
            // The fluid density is the value projected from the fluid mesh onto the particle
            // (NODAL_DENSITY_PROJECTED); the particle density follows from its mass and volume.
            // A negative projected density flags a particle without fluid data (outside the
            // mesh or before the first projection): it is then advected with the fluid only.
            const double rad = i.FastGetSolutionStepValue(RADIUS);
            const double disp_volume = (4. * Globals::Pi / 3.) * rad * rad * rad;
            const double rho_f = i.FastGetSolutionStepValue(NODAL_DENSITY_PROJECTED);
            const double rho_p = mass / disp_volume;
            // DENSITY is not necessarily a nodal variable of the particles: writing it
            // with the unchecked accessor corrupted the nodal data buffer (heap corruption
            // revealed when a particle was destroyed).  Store it only if it exists.
            if (i.SolutionStepsDataHas(DENSITY)) {
                i.FastGetSolutionStepValue(DENSITY) = rho_p;
            }
            const double mobility = disp_volume / (6.0 * Globals::Pi * mDynamicViscosity * rad); // 2 a^2 / (9 mu)
            const double buoyant_factor = (rho_f >= 0.0) ? mobility * (rho_p - rho_f) : 0.0;

            for (int k = 0; k < 3; k++){
                if (Fix_vel[k] == false){
                    vel[k] = fluid_vel[k] + buoyant_factor * mGravity[k];
                }
                else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions

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
