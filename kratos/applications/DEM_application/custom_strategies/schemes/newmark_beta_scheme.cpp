

#include "DEM_application.h"
#include "newmark_beta_scheme.h"

namespace Kratos {
    
    void NewmarkBetaScheme::AddSpheresVariables(ModelPart & r_model_part){
        
        DEMIntegrationScheme::AddSpheresVariables(r_model_part);
        //r_model_part.AddNodalSolutionStepVariable(OLD_FORCE?);
    }
    
    void NewmarkBetaScheme::AddClustersVariables(ModelPart & r_model_part){
        
        DEMIntegrationScheme::AddClustersVariables(r_model_part);
        //r_model_part.AddNodalSolutionStepVariable(OLD_FORCE?);                      
    }        
    
    
    void NewmarkBetaScheme::CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag) {
        KRATOS_TRY
        ProcessInfo& r_process_info = model_part.GetProcessInfo();
        double delta_t = r_process_info[DELTA_TIME];
        double virtual_mass_coeff = r_process_info[NODAL_MASS_COEFF];
        bool if_virtual_mass_option = (bool) r_process_info[VIRTUAL_MASS_OPTION];
        double force_reduction_factor = 1.0;
        if (if_virtual_mass_option) {
            force_reduction_factor = 1.0 - virtual_mass_coeff;
            if (virtual_mass_coeff < 0.0) KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= ", virtual_mass_coeff)
        }
        vector<unsigned int> node_partition;
        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

    #pragma omp parallel for shared(delta_t)
        for (int k = 0; k < (int) number_of_threads; k++) {
            NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k + 1];

            for (ModelPart::NodeIterator i_iterator = i_begin; i_iterator != i_end; ++i_iterator) {
                Node < 3 > & i = *i_iterator;
                if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
                array_1d<double, 3 >& vel = i.FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3 >& displ = i.FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3 >& delta_displ = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double, 3 >& coor = i.Coordinates();
                array_1d<double, 3 >& initial_coor = i.GetInitialPosition();
                array_1d<double, 3 >& force = i.FastGetSolutionStepValue(TOTAL_FORCES);                

                double mass = i.FastGetSolutionStepValue(NODAL_MASS);

                bool Fix_vel[3] = {false, false, false};

                Fix_vel[0] = i.Is(DEMFlags::FIXED_VEL_X);
                Fix_vel[1] = i.Is(DEMFlags::FIXED_VEL_Y);
                Fix_vel[2] = i.Is(DEMFlags::FIXED_VEL_Z);

                UpdateTranslationalVariables(StepFlag, i, coor, displ, delta_displ, vel, initial_coor, force, force_reduction_factor, mass, delta_t, Fix_vel);
            } //nodes in the thread
        } //threads
        KRATOS_CATCH(" ")
    }

    void NewmarkBetaScheme::UpdateTranslationalVariables(
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

        const array_1d<double, 3 >& old_force = i.FastGetSolutionStepValue(TOTAL_FORCES);
        
        for (int k = 0; k < 3; k++) {
            if (Fix_vel[k] == false) {
                double mass_inv = 1.0 / mass;
                double old_acceleration = mass_inv * old_force[k];
                double acceleration = force_reduction_factor * mass_inv * force[k];
                delta_displ[k] = delta_t * vel[k] + delta_t * delta_t * ((0.5 - mBeta) * acceleration + mBeta * old_acceleration);
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
                vel[k] += 0.5 * delta_t * (acceleration + old_acceleration);
            } else {
                delta_displ[k] = delta_t * vel[k];
                displ[k] += delta_displ[k];
                coor[k] = initial_coor[k] + displ[k];
            }
        } // dimensions  
    }

    void NewmarkBetaScheme::UpdateRotationalVariables(
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
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
                angular_velocity[k] += delta_t * angular_acceleration[k];
            } else {
                delta_rotation[k] = angular_velocity[k] * delta_t;
                rotated_angle[k] += delta_rotation[k];
            }
        }
    }            
    
    void NewmarkBetaScheme::CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration){
        
        for (int j = 0; j < 3; j++) {
            angular_acceleration[j] = moment_reduction_factor * torque[j] / moment_of_inertia;           
        }
    }
    
    
    void NewmarkBetaScheme::CalculateLocalAngularAccelerationByEulerEquations(
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
