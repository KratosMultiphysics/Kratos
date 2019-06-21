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
        array_1d<double, 3> delta_accel_init;
        array_1d<double, 3> delta_accel;
        array_1d<double, 3> temp_accel;


        if(StepFlag == 0) //initialize
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mOldAcceleration[k] = force[k] * mass_inv;

                    if (k ==0 && i.Id()==1){
                    KRATOS_WATCH("StepFlag == 0")
                    KRATOS_WATCH(mOldAcceleration[k])
                    }

                    mOldVelocity[k] = vel[k];
                    mOldDisp[k] = displ[k];

                    delta_displ[k] = vel[k] * delta_t + 0.5 * force[k] * mass_inv * delta_t * delta_t ;
                    vel[k] += 1.0 * force_reduction_factor * force[k] * mass_inv * delta_t ;

                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];

                }
            }
        }
        else if(StepFlag == 1) //predict
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    temp_accel[k] = force[k] * mass_inv;

                    delta_accel_init[k] = force[k] * mass_inv - mOldAcceleration[k]; //relative!! (A_a1-A_a2)/A_t

                    if (k ==0 && i.Id()==1){
                        KRATOS_WATCH("StepFlag == 1")
                        KRATOS_WATCH(temp_accel[k])
                        KRATOS_WATCH(delta_accel_init[k])
                    }


                    delta_displ[k] = mOldVelocity[k] * delta_t +
                                     0.5 * force[k] * mass_inv * delta_t * delta_t +
                                     0.166666*delta_accel_init[k] * delta_t * delta_t;

                    vel[k] = mOldVelocity[k] +
                             force_reduction_factor * force[k] * mass_inv * delta_t +
                             0.5 * delta_accel_init[k] * delta_t;

                    //accel pr = oldaccel + delta_accel_init*delta_t
                    //force[k] * mass_inv = mOldAcceleration[k] + delta_accel_init[k] *delta_t;

                    force[k] = (mOldAcceleration[k] + delta_accel_init[k]) * mass;

                    mOldAcceleration[k] = force[k] * mass_inv - delta_accel_init[k];

                    displ[k] = mOldDisp[k] + delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];

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

                    delta_accel[k] = force[k] * mass_inv - mOldAcceleration[k]; //absolute (A_a1-A_a2)

                    if (k ==0 && i.Id()==1){
                        KRATOS_WATCH("StepFlag == 2")
                        KRATOS_WATCH(delta_accel[k])
                    }

                    delta_displ[k] += 1/12*delta_accel[k] * delta_t * delta_t;

                    vel[k] += 5/12 * delta_accel[k] * delta_t;

                    // if (k ==1){
                    //     KRATOS_WATCH("step_ini")
                    //     KRATOS_WATCH(force[k])
                    // }

                    force[k] += delta_accel_init[k] * mass;

                    // if (k ==1){
                    //     KRATOS_WATCH(force[k])
                    // }

                    //delta_accel_init[k] = delta_accel_init[k] + delta_accel_init[k]/delta_t;

                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
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
