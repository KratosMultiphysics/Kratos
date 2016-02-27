// Project includes
#include "verlet_velocity_scheme.h"

namespace Kratos {

    void VerletVelocityScheme::AddSpheresVariables(ModelPart & r_model_part){
         DEMIntegrationScheme::AddSpheresVariables(r_model_part);}
    
    void VerletVelocityScheme::AddClustersVariables(ModelPart & r_model_part){
         DEMIntegrationScheme::AddClustersVariables(r_model_part);}

    void VerletVelocityScheme::UpdateTranslationalVariables(
            int StepFlag,
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

            double mass_inv = 1.0 / mass;
            if(StepFlag == 1) //PREDICT
            {
                for (int k = 0; k < 3; k++) 
                {
                    if (Fix_vel[k] == false) 
                    {
                      delta_displ[k] = vel[k] * delta_t + 0.5 * force[k] * mass_inv * delta_t * delta_t ;
                      displ[k] += delta_displ[k];
                      coor[k] = initial_coor[k] + displ[k];
                      vel[k] += 0.5 * force_reduction_factor * force[k] * mass_inv * delta_t ;
                    }
                    else 
                    {
                      delta_displ[k] = delta_t * vel[k];
                      displ[k] += delta_displ[k];
                      coor[k] = initial_coor[k] + displ[k];
                    }
                }  
           }

           else if(StepFlag == 2) //CORRECT
           {
             for (int k = 0; k < 3; k++) 
             {
                if (Fix_vel[k] == false) 
                {
                  vel[k] += 0.5 * force_reduction_factor * force[k] * mass_inv * delta_t ;
                }   
             }
           }


    }//VerletVelocityScheme
        
    void VerletVelocityScheme::UpdateRotationalVariables(
            int StepFlag,
            const Node < 3 > & i,
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
    
    void VerletVelocityScheme::CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration){
        
        double moment_of_inertia_inv = 1.0 / moment_of_inertia;
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] * moment_of_inertia_inv;
        }
    }
    
    
    void VerletVelocityScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
}
