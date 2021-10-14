//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Project includes
#include "central_differences_scheme.h"

namespace Kratos {

    void CentralDifferencesScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning CentralDifferencesScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void CentralDifferencesScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning CentralDifferencesScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void CentralDifferencesScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
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

        const double& alpha = i.FastGetSolutionStepValue(RAYLEIGH_ALPHA);
        const double& beta = i.FastGetSolutionStepValue(RAYLEIGH_BETA);
        const double& theta = i.FastGetSolutionStepValue(THETA_FACTOR);
        const double& g_coefficient = i.FastGetSolutionStepValue(G_COEFFICIENT);
        array_1d<double,3>& displ_old = i.FastGetSolutionStepValue(DISPLACEMENT_OLD);
        const array_1d<double,3>& internal_force = i.FastGetSolutionStepValue(INTERNAL_FORCE);
        const array_1d<double,3>& internal_force_old = i.FastGetSolutionStepValue(INTERNAL_FORCE_OLD);
        const array_1d<double,3>& external_force = i.FastGetSolutionStepValue(EXTERNAL_FORCE);
        const array_1d<double,3>& external_force_old = i.FastGetSolutionStepValue(EXTERNAL_FORCE_OLD);

        // TODO. Ignasi
        const array_1d<double,3>& nodal_mass_array = i.FastGetSolutionStepValue(NODAL_MASS_ARRAY);

        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
                delta_displ[k] = ( (2.0*(1.0+g_coefficient*delta_t)-alpha*delta_t)*nodal_mass_array[k]*displ[k]
                                + (alpha*delta_t-(1.0+g_coefficient*delta_t))*nodal_mass_array[k]*displ_old[k]
                                - delta_t*(beta+theta*delta_t)*internal_force[k]
                                + delta_t*(beta-delta_t*(1.0-theta))*internal_force_old[k]
                                + delta_t*delta_t*(theta*external_force[k]+(1.0-theta)*external_force_old[k]) ) * (1.0 / (nodal_mass_array[k]*(1.0+g_coefficient*delta_t))) - displ[k];
                displ_old[k] = displ[k];
                displ[k] = displ_old[k] + delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
                vel[k] = delta_displ[k]/delta_t;
            } else {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions

        // double mass_inv = 1.0 / (mass*(1.0+g_coefficient*delta_t));
        // for (int k = 0; k < 3; k++) {
        //     if (Fix_vel[k] == false) {
        //         delta_displ[k] = ( (2.0*(1.0+g_coefficient*delta_t)-alpha*delta_t)*mass*displ[k]
        //                         + (alpha*delta_t-(1.0+g_coefficient*delta_t))*mass*displ_old[k]
        //                         - delta_t*(beta+theta*delta_t)*internal_force[k]
        //                         + delta_t*(beta-delta_t*(1.0-theta))*internal_force_old[k]
        //                         + delta_t*delta_t*(theta*external_force[k]+(1.0-theta)*external_force_old[k]) ) * mass_inv - displ[k];
        //         displ_old[k] = displ[k];
        //         displ[k] = displ_old[k] + delta_displ[k];
        //         coor[k] = initial_coor[k] + displ[k];
        //         vel[k] = delta_displ[k]/delta_t;
        //     } else {
        //         delta_displ[k] = delta_t * vel[k];
        //         displ[k] += delta_displ[k];
        //         coor[k] = initial_coor[k] + displ[k];
        //     }
        // } // dimensions
    }

    void CentralDifferencesScheme::CalculateNewRotationalVariablesOfSpheres(
                int StepFlag,
                Node < 3 >& i,
                const double moment_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        // TODO.Ignasi: these four parameters should be taken from the process info... 
        //              or they could be member variables of the scheme...
        const double& alpha = i.FastGetSolutionStepValue(RAYLEIGH_ALPHA);
        const double& beta = i.FastGetSolutionStepValue(RAYLEIGH_BETA);
        const double& theta = i.FastGetSolutionStepValue(THETA_FACTOR);
        const double& g_coefficient = i.FastGetSolutionStepValue(G_COEFFICIENT);
        
        array_1d<double,3>& rotated_angle_old = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE_OLD);
        const array_1d<double,3>& internal_torque = i.FastGetSolutionStepValue(PARTICLE_INTERNAL_MOMENT);
        const array_1d<double,3>& internal_torque_old = i.FastGetSolutionStepValue(PARTICLE_INTERNAL_MOMENT_OLD);
        const array_1d<double,3>& external_torque = i.FastGetSolutionStepValue(PARTICLE_EXTERNAL_MOMENT);
        const array_1d<double,3>& external_torque_old = i.FastGetSolutionStepValue(PARTICLE_EXTERNAL_MOMENT_OLD);

        // TODO. Ignasi
        const array_1d<double,3>& particle_moment_intertia_array = i.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA_ARRAY);

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                delta_rotation[k] = ( (2.0*(1.0+g_coefficient*delta_t)-alpha*delta_t)*particle_moment_intertia_array[k]*rotated_angle[k]
                                    + (alpha*delta_t-(1.0+g_coefficient*delta_t))*particle_moment_intertia_array[k]*rotated_angle_old[k]
                                    - delta_t*(beta+theta*delta_t)*internal_torque[k]
                                    + delta_t*(beta-delta_t*(1.0-theta))*internal_torque_old[k]
                                    + delta_t*delta_t*(theta*external_torque[k]+(1.0-theta)*external_torque_old[k]) ) * (1.0 / (particle_moment_intertia_array[k]*(1.0+g_coefficient*delta_t))) - rotated_angle[k];
                rotated_angle_old[k] = rotated_angle[k];
                rotated_angle[k] = rotated_angle_old[k] + delta_rotation[k];
                angular_velocity[k] = delta_rotation[k]/delta_t;
            } else {
                delta_rotation[k] = delta_t * angular_velocity[k];
                rotated_angle[k] += delta_rotation[k];
            }
        } // dimensions

        // double moment_of_inertia_inv = 1.0 / (moment_of_inertia*(1.0+g_coefficient*delta_t));
        // for (int k = 0; k < 3; k++) {
        //     if (Fix_Ang_vel[k] == false) {
        //         delta_rotation[k] = ( (2.0*(1.0+g_coefficient*delta_t)-alpha*delta_t)*moment_of_inertia*rotated_angle[k]
        //                             + (alpha*delta_t-(1.0+g_coefficient*delta_t))*moment_of_inertia*rotated_angle_old[k]
        //                             - delta_t*(beta+theta*delta_t)*internal_torque[k]
        //                             + delta_t*(beta-delta_t*(1.0-theta))*internal_torque_old[k]
        //                             + delta_t*delta_t*(theta*external_torque[k]+(1.0-theta)*external_torque_old[k]) ) * moment_of_inertia_inv - rotated_angle[k];
        //         rotated_angle_old[k] = rotated_angle[k];
        //         rotated_angle[k] = rotated_angle_old[k] + delta_rotation[k];
        //         angular_velocity[k] = delta_rotation[k]/delta_t;
        //     } else {
        //         delta_rotation[k] = delta_t * angular_velocity[k];
        //         rotated_angle[k] += delta_rotation[k];
        //     }
        // } // dimensions
    }

    void CentralDifferencesScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
                int StepFlag,
                Node < 3 >& i,
                const array_1d<double, 3 > moments_of_inertia,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        array_1d<double, 3 >& local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);

        array_1d<double, 3 > local_angular_acceleration, local_torque, angular_acceleration;

        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        CalculateLocalAngularAccelerationByEulerEquations(local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);
        GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_acceleration, angular_acceleration);

        // NOTE: this is from symplectic euler...
        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);

        double ang = DEM_INNER_PRODUCT_3(delta_rotation, delta_rotation);

        if (ang) {
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        } //if ang
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
    }

    void CentralDifferencesScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
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

    void CentralDifferencesScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
