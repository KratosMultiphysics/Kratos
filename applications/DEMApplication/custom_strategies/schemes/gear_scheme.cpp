// Project includes
#include "gear_scheme.h"

namespace Kratos {

    void GearScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning GearScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GearScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning GearScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GearScheme::UpdateTranslationalVariables(
            int StepFlag,
            Node < 3 >& i,
            array_1d<double, 3 >& coor,
            array_1d<double, 3 >& displ,
            array_1d<double, 3 >& delta_displ,
            array_1d<double, 3 >& vel,
            const array_1d<double, 3 >& initial_coor,
            array_1d<double, 3 >& force,
            const double force_reduction_factor,
            const double mass,
            const double delta_t,
            const bool Fix_vel[3]) {

        double mass_inv = 1.0 / mass;
        array_1d<double, 3> accel_est_pred;

        if(StepFlag == 1) //predict
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    delta_displ[k] = vel[k] * delta_t +
                                     0.5 * force[k] * mass_inv * delta_t * delta_t +
                                     0.166666*mDeltaAccel[k] * delta_t * delta_t * delta_t;

                    vel[k] += 0.5*force_reduction_factor * force[k] * mass_inv * delta_t +
                             0.5 * mDeltaAccel[k] * delta_t * delta_t;

                    mOldAcceleration[k] = force[k] * mass_inv + mDeltaAccel[k] * delta_t;

                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];

                    // if (k ==0 && i.Id()==1){
                    //     KRATOS_WATCH("StepFlag == 1")
                    //     KRATOS_WATCH(mDeltaAccel[k])
                    //     KRATOS_WATCH(delta_displ[k])
                    //     KRATOS_WATCH(vel[k])
                    //     KRATOS_WATCH(mOldAcceleration[k])
                    // }

                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 2) //correct
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    accel_est_pred[k] = force[k] * mass_inv - mOldAcceleration[k];

                    //delta_displ[k] += 1/12*mDeltaAccel[k] * delta_t * delta_t;

                    vel[k] += 0.5*force_reduction_factor * force[k] * mass_inv * delta_t +
                              5/12 * mDeltaAccel[k] * delta_t;
                    //vel[k] += 0.5*force_reduction_factor * force[k] * mass_inv * delta_t;

                    mDeltaAccel[k] += accel_est_pred[k]/delta_t;

                    //displ[k] += delta_displ[k];
                    //coor[k] = initial_coor[k] + displ[k];

                    // if (k ==0 && i.Id()==1){
                    //     KRATOS_WATCH("StepFlag == 2")
                    //     KRATOS_WATCH(accel_est_pred[k])
                    //     KRATOS_WATCH(mDeltaAccel[k])
                    //     KRATOS_WATCH(delta_displ[k])
                    //     KRATOS_WATCH(vel[k])
                    //     KRATOS_WATCH(mOldAcceleration[k])
                    // }
                }
            }
        } // dimensions
    }

    void GearScheme::CalculateNewRotationalVariablesOfSpheres(
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

        array_1d<double, 3 > angular_acceleration;
        CalculateLocalAngularAcceleration(moment_of_inertia, torque, moment_reduction_factor, angular_acceleration);

        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);
    }

    void GearScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
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

        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);

        if (StepFlag == 1) //PREDICT
        {
            double ang = DEM_INNER_PRODUCT_3(delta_rotation, delta_rotation);

            if (ang) {
                GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
            } //if ang
        }
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
    }

    void GearScheme::UpdateRotationalVariables(
                int StepFlag,
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        if (StepFlag == 1) //PREDICT
        {
             for (int k = 0; k < 3; k++) {
                 if (Fix_Ang_vel[k] == false) {
                     delta_rotation[k] = angular_velocity[k] * delta_t + 0.5 * delta_t * delta_t * angular_acceleration[k];
                     rotated_angle[k] += delta_rotation[k];
                     angular_velocity[k] += 0.5 * angular_acceleration[k] * delta_t;
                } else {
                     delta_rotation[k] = angular_velocity[k] * delta_t;
                     rotated_angle[k] += delta_rotation[k];
                }
            }
        }

        else if(StepFlag == 2) //CORRECT
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_Ang_vel[k] == false) {
                    angular_velocity[k] += 0.5 * angular_acceleration[k] * delta_t;
                }
            }
        }//CORRECT
    }

    void GearScheme::CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) {

        double moment_of_inertia_inv = 1.0 / moment_of_inertia;
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] * moment_of_inertia_inv;
        }
    }

    void GearScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
