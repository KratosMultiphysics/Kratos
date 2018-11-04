// Project includes
#include "quaternion_integration_scheme.h"

namespace Kratos {

    void QuaternionIntegrationScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning QuaternionIntegrationScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void QuaternionIntegrationScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
//         if(verbose) KRATOS_INFO("DEM") << "Assigning QuaternionIntegrationScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void QuaternionIntegrationScheme::UpdateTranslationalVariables(
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
            const bool Fix_vel[3])
    {
        KRATOS_THROW_ERROR(std::runtime_error, "This scheme (QuaternionIntegrationScheme) should not calculate translational motion, so the function (QuaternionIntegrationScheme::UpdateTranslationalVariables) shouldn't be accessed", 0);
    }

    void QuaternionIntegrationScheme::CalculateNewRotationalVariablesOfSpheres(
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

        Quaternion<double  > Orientation = Quaternion<double>::Identity();

        array_1d<double, 3 >& half_local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_AUX_ANGULAR_VELOCITY);
        Quaternion<double  >& half_Orientation            = i.FastGetSolutionStepValue(AUX_ORIENTATION);

        array_1d<double, 3 > local_angular_acceleration, local_torque;

        array_1d<double, 3 > moments_of_inertia;
        moments_of_inertia[0] = moment_of_inertia;
        moments_of_inertia[1] = moment_of_inertia;
        moments_of_inertia[2] = moment_of_inertia;

        array_1d<double, 3 > global_torque = ZeroVector(3);

        for (int j = 0; j < 3; j++) {
            if (Fix_Ang_vel[j] == false) {
                global_torque[j] = torque[j];
            }
        }

        if (StepFlag != 1 && StepFlag != 2) {
            array_1d<double, 3 > local_angular_velocity, quarter_angular_velocity;

            CalculateLocalAngularAcceleration(moment_of_inertia, global_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(quarter_angular_velocity)    = angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            noalias(half_local_angular_velocity) = angular_velocity + 0.5  * local_angular_acceleration * delta_t;

            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, half_Orientation, rotation_aux);

            GeometryFunctions::QuaternionVectorGlobal2Local(half_Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(half_local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(local_angular_velocity) = angular_velocity + local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(half_Orientation, local_angular_velocity, angular_velocity);

            UpdateRotatedAngle(rotated_angle, delta_rotation, angular_velocity, delta_t);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
        }

        else if (StepFlag == 1) { //PREDICT
            array_1d<double, 3 > quarter_angular_velocity, local_angular_velocity;

            CalculateLocalAngularAcceleration(moment_of_inertia, global_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(quarter_angular_velocity)    = angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            noalias(half_local_angular_velocity) = angular_velocity + 0.5  * local_angular_acceleration * delta_t;

            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, half_Orientation, rotation_aux);
        }//if StepFlag == 1

        else if (StepFlag == 2) { //CORRECT
            array_1d<double, 3 > local_angular_velocity;

            GeometryFunctions::QuaternionVectorGlobal2Local(half_Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(half_local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(local_angular_velocity) = angular_velocity + local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(half_Orientation, local_angular_velocity, angular_velocity);

            UpdateRotatedAngle(rotated_angle, delta_rotation, angular_velocity, delta_t);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
        }//if StepFlag == 2
    }

    void QuaternionIntegrationScheme::CalculateNewRotationalVariablesOfRigidBodyElements(
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

        array_1d<double, 3 >& local_angular_velocity      = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);
        array_1d<double, 3 >& half_local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_AUX_ANGULAR_VELOCITY);
        Quaternion<double  >& half_Orientation            = i.FastGetSolutionStepValue(AUX_ORIENTATION);
        array_1d<double, 3 > local_angular_acceleration, local_torque;

        array_1d<double, 3 > global_torque = ZeroVector(3);

        for (int j = 0; j < 3; j++) {
            if (Fix_Ang_vel[j] == false) {
                global_torque[j] = torque[j];
            }
        }

        if (StepFlag != 1 && StepFlag != 2) {
            array_1d<double, 3 > quarter_local_angular_velocity, quarter_angular_velocity;

            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(quarter_local_angular_velocity) = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            noalias(half_local_angular_velocity)    = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, half_Orientation, rotation_aux);

            GeometryFunctions::QuaternionVectorGlobal2Local(half_Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(half_local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(local_angular_velocity) += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(half_Orientation, local_angular_velocity, angular_velocity);

            UpdateRotatedAngle(rotated_angle, delta_rotation, angular_velocity, delta_t);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
        }

        else if (StepFlag == 1) { //PREDICT
            array_1d<double, 3 > quarter_local_angular_velocity, quarter_angular_velocity;

            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(local_angular_velocity,moments_of_inertia,local_torque,moment_reduction_factor,local_angular_acceleration);

            noalias(quarter_local_angular_velocity) = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            noalias(half_local_angular_velocity)    = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, half_Orientation, rotation_aux);
        }//if StepFlag == 1

        else if (StepFlag == 2) { //CORRECT
            GeometryFunctions::QuaternionVectorGlobal2Local(half_Orientation, global_torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(half_local_angular_velocity, moments_of_inertia, local_torque, moment_reduction_factor, local_angular_acceleration);

            noalias(local_angular_velocity) += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(half_Orientation, local_angular_velocity, angular_velocity);

            UpdateRotatedAngle(rotated_angle, delta_rotation, angular_velocity, delta_t);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
        }//if StepFlag == 2
    }

    void QuaternionIntegrationScheme::UpdateRotatedAngle(
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t) {

        noalias(delta_rotation) = angular_velocity * delta_t;
        noalias(rotated_angle) += delta_rotation;
    }

    void QuaternionIntegrationScheme::CalculateLocalAngularAcceleration(
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration) {

        double moment_of_inertia_inv = 1.0 / moment_of_inertia;
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] * moment_of_inertia_inv;
        }
    }

    void QuaternionIntegrationScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
