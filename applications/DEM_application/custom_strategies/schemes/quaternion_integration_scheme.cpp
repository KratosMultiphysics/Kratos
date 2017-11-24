// Project includes
#include "quaternion_integration_scheme.h"

namespace Kratos {

    void QuaternionIntegrationScheme::SetIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        if(verbose) std::cout << "\nAssigning QuaternionIntegrationScheme to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void QuaternionIntegrationScheme::UpdateTranslationalVariables(
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
        KRATOS_THROW_ERROR(std::runtime_error, "This scheme (QuaternionIntegrationScheme) should not calculate translational motion, so the function (QuaternionIntegrationScheme::UpdateTranslationalVariables) shouldn't be accessed", 0);
    }

    void QuaternionIntegrationScheme::CalculateNewRotationalVariables(
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

        array_1d<double, 3 > local_angular_acceleration, local_torque, quarter_local_angular_velocity, quarter_angular_velocity, AuxAngularVelocity;
        array_1d<double, 3 > & local_angular_velocity  = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);
        Quaternion<double  > & AuxOrientation          = i.FastGetSolutionStepValue(AUX_ORIENTATION);
        array_1d<double, 3 > & LocalAuxAngularVelocity = i.FastGetSolutionStepValue(LOCAL_AUX_ANGULAR_VELOCITY);
        Quaternion<double  > Orientation               = Quaternion<double>(1.0, 0.0, 0.0, 0.0);
        
        array_1d<double, 3 > moments_of_inertia;
        moments_of_inertia[0] = moment_of_inertia;
        moments_of_inertia[1] = moment_of_inertia;
        moments_of_inertia[2] = moment_of_inertia;

        if (StepFlag != 1 && StepFlag != 2) {
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,local_angular_velocity,moments_of_inertia,local_torque,moment_reduction_factor,local_angular_acceleration);
            quarter_local_angular_velocity = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            LocalAuxAngularVelocity        = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, AuxOrientation, rotation_aux);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(AuxOrientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,LocalAuxAngularVelocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);
            local_angular_velocity += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }
        
        else if (StepFlag == 1) { //PREDICT
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,local_angular_velocity,moments_of_inertia,local_torque,moment_reduction_factor,local_angular_acceleration);
            quarter_local_angular_velocity = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            LocalAuxAngularVelocity        = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, AuxOrientation, rotation_aux);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }//if StepFlag == 1
                    
        else if (StepFlag == 2) { //CORRECT
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(AuxOrientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,LocalAuxAngularVelocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);
            local_angular_velocity += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }//if StepFlag == 2
    }
    
    void QuaternionIntegrationScheme::CalculateNewRotationalVariables(
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

        array_1d<double, 3 > local_angular_acceleration, local_torque, quarter_local_angular_velocity, quarter_angular_velocity, AuxAngularVelocity;
        array_1d<double, 3 > & local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);
        Quaternion<double  > & AuxOrientation = i.FastGetSolutionStepValue(AUX_ORIENTATION);
        array_1d<double, 3 > & LocalAuxAngularVelocity = i.FastGetSolutionStepValue(LOCAL_AUX_ANGULAR_VELOCITY);

        if (StepFlag != 1 && StepFlag != 2) {
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,local_angular_velocity,moments_of_inertia,local_torque,moment_reduction_factor,local_angular_acceleration);
            quarter_local_angular_velocity = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            LocalAuxAngularVelocity        = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, AuxOrientation, rotation_aux);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(AuxOrientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,LocalAuxAngularVelocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);
            local_angular_velocity += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }
        
        else if (StepFlag == 1) { //PREDICT
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,local_angular_velocity,moments_of_inertia,local_torque,moment_reduction_factor,local_angular_acceleration);
            quarter_local_angular_velocity = local_angular_velocity + 0.25 * local_angular_acceleration * delta_t;
            LocalAuxAngularVelocity        = local_angular_velocity + 0.5  * local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, quarter_local_angular_velocity, quarter_angular_velocity);
            array_1d<double, 3 > rotation_aux = 0.5 * quarter_angular_velocity * delta_t;
            GeometryFunctions::UpdateOrientation(Orientation, AuxOrientation, rotation_aux);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }//if StepFlag == 1
                    
        else if (StepFlag == 2) { //CORRECT
            //Angular velocity and torques are saved in the local framework:
            GeometryFunctions::QuaternionVectorGlobal2Local(AuxOrientation, torque, local_torque);
            CalculateLocalAngularAccelerationByEulerEquations(i,LocalAuxAngularVelocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);
            local_angular_velocity += local_angular_acceleration * delta_t;
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
            GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
            UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        }//if StepFlag == 2
    }

    void QuaternionIntegrationScheme::UpdateRotationalVariables(
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
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
                angular_velocity[k] += delta_t * angular_acceleration[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }

    void QuaternionIntegrationScheme::UpdateRotationalVariablesOfCluster(
                Node < 3 >& i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
        }
        
        array_1d<double, 3 > angular_velocity_aux;
        
        double LocalTensorInv[3][3];
        GeometryFunctions::ConstructInvLocalTensor(moments_of_inertia, LocalTensorInv);
        GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
        UpdateAngularVelocity(Orientation, LocalTensorInv, angular_momentum, angular_velocity_aux);
        for (int j = 0; j < 3; j++) {
            if (Fix_Ang_vel[j] == false){
                angular_velocity[j] = angular_velocity_aux[j];
            }
        }
    }
    
    void QuaternionIntegrationScheme::UpdateRotationalVariables(
                Node < 3 >& i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        
        delta_rotation = angular_velocity * delta_t;
        rotated_angle += delta_rotation;
    }
    
    void QuaternionIntegrationScheme::QuaternionCalculateMidAngularVelocities(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                const double dt,
                const array_1d<double, 3>& InitialAngularVel,
                array_1d<double, 3>& FinalAngularVel) {
        
        array_1d<double, 3 > aux = InitialAngularVel;
        DEM_MULTIPLY_BY_SCALAR_3(aux, dt);
        array_1d<double, 3 > TempDeltaRotation = aux;

        Quaternion<double> TempOrientation;
        double GlobalTensorInv[3][3];
            
        GeometryFunctions::UpdateOrientation(Orientation, TempOrientation, TempDeltaRotation);
        GeometryFunctions::QuaternionTensorLocal2Global(TempOrientation, LocalTensorInv, GlobalTensorInv);
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensorInv, angular_momentum, FinalAngularVel);
    }
    
    void QuaternionIntegrationScheme::UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity) {
        
        double GlobalTensorInv[3][3];
        
        GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensorInv, GlobalTensorInv);
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensorInv, angular_momentum, angular_velocity);
    }

    void QuaternionIntegrationScheme::CalculateLocalAngularAcceleration(
                Node < 3 >& i,
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
                Node < 3 >& i,
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

    void QuaternionIntegrationScheme::CalculateAngularVelocityRK(
                                    const Quaternion<double  >& Orientation,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& angular_momentum,
                                    array_1d<double, 3 >& angular_velocity,
                                    const double delta_t,
                                    const bool Fix_Ang_vel[3]) {
            
            double LocalTensorInv[3][3];
            
            GeometryFunctions::ConstructInvLocalTensor(moments_of_inertia, LocalTensorInv);
            
            array_1d<double, 3 > angular_velocity1 = angular_velocity;
            array_1d<double, 3 > angular_velocity2, angular_velocity3, angular_velocity4;

            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum, 0.5*delta_t, angular_velocity1, angular_velocity2);
            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum, 0.5*delta_t, angular_velocity2, angular_velocity3);
            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum,     delta_t, angular_velocity3, angular_velocity4);

            for (int j = 0; j < 3; j++) {
                if (Fix_Ang_vel[j] == false){
                    angular_velocity[j] = 0.16666666666666667 * (angular_velocity1[j] + 2*angular_velocity2[j] + 2*angular_velocity3[j] + angular_velocity4[j]);
                }
            }
    }
} //namespace Kratos
