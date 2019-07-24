#include "translational_RK4_scheme.h"

namespace Kratos {

    void TranslationalRungeKuttaScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void TranslationalRungeKuttaScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void TranslationalRungeKuttaScheme::UpdateTranslationalVariables(
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

        double mass_inv = 1.0 / mass;
        KRATOS_WATCH(i)

        if(StepFlag == 0) //Init
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    mInitialDispl[k] = displ[k];

                    KRATOS_WATCH(StepFlag)

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k];

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k) * 0.5 ;
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];

                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 1) //Step1
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    KRATOS_WATCH(StepFlag)

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k) * 0.5;

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k) * 0.5 ;
                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 2) //Step2
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    KRATOS_WATCH(StepFlag)


                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k) * 0.5;

                    KRATOS_WATCH(mRungeKuttaL(StepFlag,k))
                    KRATOS_WATCH(mRungeKuttaK(StepFlag,k))

                    delta_displ[k] = mRungeKuttaK(StepFlag,k);
                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            }
        }
        else if(StepFlag == 3) //Step3
        {
            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {

                    KRATOS_WATCH(StepFlag)
                    KRATOS_WATCH(mRungeKuttaL)
                    KRATOS_WATCH(mRungeKuttaK)

                    mRungeKuttaL(StepFlag,k) = delta_t * force[k] * mass_inv;
                    KRATOS_WATCH("1")
                    mRungeKuttaK(StepFlag,k) = delta_t * vel[k] +  delta_t * mRungeKuttaL(StepFlag-1,k);
                    KRATOS_WATCH("2")

                    delta_displ[k] = 1.0/6.0 * (mRungeKuttaK(0,k) + 2.0 * mRungeKuttaK(1,k) + 2.0 * mRungeKuttaK(2,k) + mRungeKuttaK(3,k));
                    KRATOS_WATCH(delta_displ[k])

                    displ[k] = mInitialDispl[k] + delta_displ[k];
                    KRATOS_WATCH(displ[k])

                    coor[k] = initial_coor[k] + displ[k];
                    KRATOS_WATCH(coor[k])

                    vel[k] += 1.0/6.0 * (mRungeKuttaL(0,k) + 2.0 * mRungeKuttaL(1,k) + 2.0 * mRungeKuttaL(2,k) + mRungeKuttaL(3,k));
                    KRATOS_WATCH(vel[k])

                }
            }
        }
    }

    void TranslationalRungeKuttaScheme::CalculateNewRotationalVariablesOfSpheres(
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

    void TranslationalRungeKuttaScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
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

    void TranslationalRungeKuttaScheme::UpdateRotationalVariables(
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
        }
    }

    void TranslationalRungeKuttaScheme::CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) {

        double moment_of_inertia_inv = 1.0 / moment_of_inertia;
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] * moment_of_inertia_inv;
        }
    }

    void TranslationalRungeKuttaScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
}
