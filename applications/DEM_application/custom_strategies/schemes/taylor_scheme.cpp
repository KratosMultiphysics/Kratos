// Project includes
#include "taylor_scheme.h"

namespace Kratos {
    
        void TaylorScheme::SetIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
            if(verbose) std::cout << "\nAssigning TaylorScheme to properties " << pProp->Id() << std::endl;
            pProp->SetValue(DEM_INTEGRATION_SCHEME_POINTER, this->CloneShared());
        }

        /*void TaylorScheme::AddSpheresVariables(ModelPart & r_model_part, bool TRotationOption){
            DEMIntegrationScheme::AddSpheresVariables(r_model_part, TRotationOption);}

        void TaylorScheme::AddClustersVariables(ModelPart & r_model_part, bool TRotationOption){
            DEMIntegrationScheme::AddClustersVariables(r_model_part, TRotationOption);}*/

        void TaylorScheme::UpdateTranslationalVariables(
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
                const bool Fix_vel[3]) {

            for (int k = 0; k < 3; k++) {
                if (Fix_vel[k] == false) {
                    delta_displ[k] = delta_t * (vel [k] + (0.5 * delta_t / mass) * force[k]);
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                    vel[k] += (delta_t / mass) * force[k];
                } else {
                    delta_displ[k] = delta_t * vel[k];
                    displ[k] += delta_displ[k];
                    coor[k] = initial_coor[k] + displ[k];
                }
            } // dimensions                                     

        }

        void TaylorScheme::UpdateRotationalVariables(
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
                    delta_rotation[k] = delta_t * (angular_velocity[k] + (0.5 * delta_t * angular_acceleration[k]));
                    rotated_angle[k] += delta_rotation[k];
                    angular_velocity[k] += delta_t * angular_acceleration[k];
                } else {
                    delta_rotation[k] = angular_velocity[k] * delta_t;
                    rotated_angle[k] += delta_rotation[k];
                }
            }
        }
        
        void TaylorScheme::UpdateRotationalVariablesOfCluster(
                const Node < 3 > & i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3])  {

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
    
    void TaylorScheme::UpdateRotationalVariables(
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

        for (int k = 0; k < 3; k++) {
            if (Fix_Ang_vel[k] == false) {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }
    
    void TaylorScheme::QuaternionCalculateMidAngularVelocities(
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
    
    void TaylorScheme::UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity) {
        
        double GlobalTensorInv[3][3];
        
        GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensorInv, GlobalTensorInv);
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensorInv, angular_momentum, angular_velocity);
    }

        void TaylorScheme::CalculateLocalAngularAcceleration(
                const Node < 3 > & i,
                const double moment_of_inertia,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& angular_acceleration)  {

            for (int j = 0; j < 3; j++) {
                angular_acceleration[j] = moment_reduction_factor * torque[j] / moment_of_inertia;
            }
        }

        void TaylorScheme::CalculateLocalAngularAccelerationByEulerEquations(
                const Node < 3 > & i,
                const array_1d<double, 3 >& local_angular_velocity,
                const array_1d<double, 3 >& moments_of_inertia,
                const array_1d<double, 3 >& local_torque,
                const double moment_reduction_factor,
                array_1d<double, 3 >& local_angular_acceleration)  {

            for (int j = 0; j < 3; j++) {
                local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
                local_angular_acceleration[j] = local_angular_acceleration[j] * moment_reduction_factor;
            }
        }
        
        void TaylorScheme::CalculateAngularVelocityRK(
                                    const Quaternion<double  >& Orientation,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& angular_momentum,
                                    array_1d<double, 3 > & angular_velocity,
                                    const double delta_t,
                                    const bool Fix_Ang_vel[3]) {
                                        
            double dt = delta_t;
            
            double LocalTensorInv[3][3];
            
            GeometryFunctions::ConstructInvLocalTensor(moments_of_inertia, LocalTensorInv);
            
            array_1d<double, 3 > angular_velocity1 = angular_velocity;
            array_1d<double, 3 > angular_velocity2, angular_velocity3, angular_velocity4;

            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum, 0.5*dt, angular_velocity1, angular_velocity2);
            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum, 0.5*dt, angular_velocity2, angular_velocity3);
            QuaternionCalculateMidAngularVelocities(Orientation, LocalTensorInv, angular_momentum,     dt, angular_velocity3, angular_velocity4);

            for (int j = 0; j < 3; j++) {
                if (Fix_Ang_vel[j] == false){
                    angular_velocity[j] = 0.16666666666666667 * (angular_velocity1[j] + 2*angular_velocity2[j] + 2*angular_velocity3[j] + angular_velocity4[j]);
                }
            }
    }

} // namespace Kratos.
