//        
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/cluster3D.h"
#include "dem_integration_scheme.h"
#include "DEM_application_variables.h"


namespace Kratos {

    DEMIntegrationScheme::DEMIntegrationScheme(){}
    DEMIntegrationScheme::~DEMIntegrationScheme(){}
    
    void DEMIntegrationScheme::SetIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        //if(verbose) std::cout << "\nAssigning DEMDiscontinuumConstitutiveLaw to properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    /*void DEMIntegrationScheme::AddSpheresVariables(ModelPart & r_model_part, bool TRotationOption){
        
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES);
        r_model_part.AddNodalSolutionStepVariable(NODAL_MASS);   
        
        if (TRotationOption){
            r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA); 
            r_model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY); 
            r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT); 
            r_model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE); 
            r_model_part.AddNodalSolutionStepVariable(DELTA_ROTATION);         
        }
    }
    
    void DEMIntegrationScheme::AddClustersVariables(ModelPart & r_model_part, bool TRotationOption){
        
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES);
        r_model_part.AddNodalSolutionStepVariable(NODAL_MASS);
        r_model_part.AddNodalSolutionStepVariable(ORIENTATION);
        r_model_part.AddNodalSolutionStepVariable(ORIENTATION_REAL); // JIG: SHOULD BE REMOVED IN THE FUTURE
        r_model_part.AddNodalSolutionStepVariable(ORIENTATION_IMAG); // JIG: SHOULD BE REMOVED IN THE FUTURE
        r_model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA);
        r_model_part.AddNodalSolutionStepVariable(ANGULAR_MOMENTUM); 
        
        if (TRotationOption){
            r_model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(LOCAL_AUX_ANGULAR_VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(LOCAL_ANGULAR_VELOCITY);
            r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT);
            r_model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE);
            r_model_part.AddNodalSolutionStepVariable(AUX_ORIENTATION);
            r_model_part.AddNodalSolutionStepVariable(DELTA_ROTATION);
        }
    }    */
    
    void DEMIntegrationScheme::CalculateTranslationalMotionOfNode(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        array_1d<double, 3 >& vel = i.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3 >& displ = i.FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 >& delta_displ = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        array_1d<double, 3 >& coor = i.Coordinates();
        array_1d<double, 3 >& initial_coor = i.GetInitialPosition();
        array_1d<double, 3 >& force = i.FastGetSolutionStepValue(TOTAL_FORCES);
        
        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(force, "NAN in Force in Integration Scheme");
        #endif  
        
        double mass = i.FastGetSolutionStepValue(NODAL_MASS);                   

        bool Fix_vel[3] = {false, false, false};

        Fix_vel[0] = i.Is(DEMFlags::FIXED_VEL_X);
        Fix_vel[1] = i.Is(DEMFlags::FIXED_VEL_Y);
        Fix_vel[2] = i.Is(DEMFlags::FIXED_VEL_Z);

        UpdateTranslationalVariables(StepFlag, i, coor, displ, delta_displ, vel, initial_coor, force, force_reduction_factor, mass, delta_t, Fix_vel);       
    }
    
    void DEMIntegrationScheme::CalculateRotationalMotionOfNode(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
    
        double moment_of_inertia = i.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
        array_1d<double, 3 >& angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3 >& torque = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
        array_1d<double, 3 >& rotated_angle = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        array_1d<double, 3 >& delta_rotation = i.FastGetSolutionStepValue(DELTA_ROTATION);
        
        #ifdef KRATOS_DEBUG
        DemDebugFunctions::CheckIfNan(torque, "NAN in Torque in Integration Scheme");
        #endif

        bool Fix_Ang_vel[3] = {false, false, false};
        Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
        Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
        Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

        array_1d<double, 3 > angular_acceleration;                    
        CalculateLocalAngularAcceleration(i, moment_of_inertia, torque, force_reduction_factor,angular_acceleration);

        UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);
    }
    
    void DEMIntegrationScheme::Move(Node<3> & i, const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;
        CalculateTranslationalMotionOfNode(i, delta_t, force_reduction_factor, StepFlag);  
        
        if(rotation_option) {
            CalculateRotationalMotionOfNode(i, delta_t, force_reduction_factor, StepFlag); 
        }                        
    }
    
    void DEMIntegrationScheme::MoveCluster(Cluster3D* cluster_element, Node<3> & i, const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
        CalculateTranslationalMotionOfNode(i, delta_t, force_reduction_factor, StepFlag);   
        if(rotation_option) {
            RotateClusterNode(i, delta_t, force_reduction_factor, StepFlag);                
            cluster_element->UpdatePositionOfSpheres();
        }  
        else {
            cluster_element->UpdateLinearDisplacementAndVelocityOfSpheres();
        }
    }    
        
    void DEMIntegrationScheme::UpdateTranslationalVariables(
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
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateTranslationalVariables) shouldn't be accessed, use derived class instead", 0);
    }    
    
    void DEMIntegrationScheme::UpdateLocalAngularVelocity(
                const Node < 3 > & i,
                array_1d<double, 3 >& partial_local_angular_velocity,
                array_1d<double, 3 >& local_angular_velocity,
                array_1d<double, 3 >& local_angular_acceleration,
                double dt,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateLocalAngularVelocity) shouldn't be accessed, use derived class instead", 0);            
    }
    
    void DEMIntegrationScheme::CalculateLocalAngularAcceleration(
                                const Node < 3 > & i,
                                const double moment_of_inertia,
                                const array_1d<double, 3 >& torque, 
                                const double moment_reduction_factor,
                                array_1d<double, 3 >& angular_acceleration){
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateLocalAngularAcceleration) shouldn't be accessed, use derived class instead", 0);            
    }
    
    void DEMIntegrationScheme::UpdateRotationalVariables(
                int StepFlag,
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                array_1d<double, 3 >& angular_acceleration,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);
    }

    void DEMIntegrationScheme::UpdateRotationalVariablesOfCluster(
                const Node < 3 > & i,
                const array_1d<double, 3 >& moments_of_inertia,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                Quaternion<double  >& Orientation,
                const array_1d<double, 3 >& angular_momentum,
                array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariablesOfCluster) shouldn't be accessed, use derived class instead", 0);
    }
    
    void DEMIntegrationScheme::UpdateRotationalVariables(
                const Node < 3 > & i,
                array_1d<double, 3 >& rotated_angle,
                array_1d<double, 3 >& delta_rotation,
                const array_1d<double, 3 >& angular_velocity,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);
    }
    
    void DEMIntegrationScheme::QuaternionCalculateMidAngularVelocities(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                const double dt,
                const array_1d<double, 3>& InitialAngularVel,
                array_1d<double, 3>& FinalAngularVel) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::QuaternionCalculateMidAngularVelocities) shouldn't be accessed, use derived class instead", 0);
    }
    
    void DEMIntegrationScheme::UpdateAngularVelocity(
                const Quaternion<double>& Orientation,
                const double LocalTensorInv[3][3],
                const array_1d<double, 3>& angular_momentum,
                array_1d<double, 3>& angular_velocity) {
        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateAngularVelocity) shouldn't be accessed, use derived class instead", 0);
    }    
    
    void DEMIntegrationScheme::CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag) {

        KRATOS_TRY

        ProcessInfo& r_process_info = model_part.GetProcessInfo();
        double delta_t = r_process_info[DELTA_TIME];
        vector<unsigned int> node_partition;
        double virtual_mass_coeff = r_process_info[NODAL_MASS_COEFF];
        bool virtual_mass_option = (bool) r_process_info[VIRTUAL_MASS_OPTION];
        double moment_reduction_factor = 1.0;
        if (virtual_mass_option) {
            moment_reduction_factor = 1.0 - virtual_mass_coeff;
            if ((moment_reduction_factor > 1.0) || (moment_reduction_factor < 0.0)) {
                KRATOS_THROW_ERROR(std::runtime_error, "The virtual mass coefficient is either larger than 1 or negative: virtual_mass_coeff= ", virtual_mass_coeff)
            }
        }      

        #pragma omp parallel for shared(delta_t)
        for (int k = 0; k < (int) pNodes.size(); k++) {
            ModelPart::NodeIterator i_iterator = pNodes.ptr_begin() + k;
            Node < 3 > & i = *i_iterator;
            if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
            CalculateRotationalMotionOfNode(i, delta_t, moment_reduction_factor, StepFlag);            
        }//for Node                  

        KRATOS_CATCH(" ")

    }//rotational_motion
    
    void DEMIntegrationScheme::CalculateLocalAngularAccelerationByEulerEquations(
                                    const Node < 3 > & i,
                                    const array_1d<double, 3 >& local_angular_velocity,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& local_torque, 
                                    const double moment_reduction_factor,
                                    array_1d<double, 3 >& local_angular_acceleration){
            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateLocalAngularAccelerationByEulerEquations) shouldn't be accessed, use derived class instead", 0);                        
    }  

    void DEMIntegrationScheme::CalculateAngularVelocityRK(
                                    const Quaternion<double  >& Orientation,
                                    const array_1d<double, 3 >& moments_of_inertia,
                                    const array_1d<double, 3 >& angular_momentum,
                                    array_1d<double, 3 >& angular_velocity,
                                    const double delta_t,
                                    const bool Fix_Ang_vel[3]) {
            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::CalculateAngularVelocityRK) shouldn't be accessed, use derived class instead", 0);                        
    }  
   
    void DEMIntegrationScheme::CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part, int StepFlag) { //must be done AFTER the translational motion!

        KRATOS_TRY

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ElementsArrayType::iterator ElementIterator;
        ProcessInfo& r_process_info = rcluster_model_part.GetProcessInfo();

        double delta_t = r_process_info[DELTA_TIME];
        double virtual_mass_coeff = r_process_info[NODAL_MASS_COEFF];
        bool virtual_mass_option = (bool) r_process_info[VIRTUAL_MASS_OPTION];
        double moment_reduction_factor = 1.0;
        if (virtual_mass_option) {
            moment_reduction_factor = 1.0 - virtual_mass_coeff;
            if ((moment_reduction_factor > 1.0) || (moment_reduction_factor < 0.0)) {
                KRATOS_THROW_ERROR(std::runtime_error, "The virtual mass coefficient is either larger than 1 or negative: virtual_mass_coeff= ", virtual_mass_coeff)
            }
        }

        vector<unsigned int> element_partition;
        ElementsArrayType& pElements = rcluster_model_part.GetCommunicator().LocalMesh().Elements();
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), element_partition);

        #pragma omp parallel for shared(delta_t)        
        for (int k = 0; k < (int) OpenMPUtils::GetNumThreads(); k++) {
            ElementIterator i_begin = pElements.ptr_begin() + element_partition[k];
            ElementIterator i_end = pElements.ptr_begin() + element_partition[k + 1];

            for (ElementsArrayType::iterator it = i_begin; it != i_end; ++it) {
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                Node<3> & i = cluster_element.GetGeometry()[0];

                RotateClusterNode(i, delta_t, moment_reduction_factor, StepFlag);                
                cluster_element.UpdatePositionOfSpheres();
            } //for Elements
        } //for number of threads
        KRATOS_CATCH(" ")
    }
    
    void DEMIntegrationScheme::RotateClusterNode(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {
        KRATOS_TRY
        array_1d<double, 3 > & moments_of_inertia = i.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        array_1d<double, 3 > & angular_momentum = i.FastGetSolutionStepValue(ANGULAR_MOMENTUM);
        array_1d<double, 3 > & angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3 > & local_angular_velocity = i.FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY);
        array_1d<double, 3 > & torque = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
        array_1d<double, 3 > & rotated_angle = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        Quaternion<double  > & Orientation = i.FastGetSolutionStepValue(ORIENTATION);
        array_1d<double, 3 > & delta_rotation = i.FastGetSolutionStepValue(DELTA_ROTATION);

        bool Fix_Ang_vel[3] = {false, false, false};

        Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
        Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
        Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);
                
        if (StepFlag != 1 && StepFlag != 2) {

            array_1d<double, 3 > angular_momentum_aux;
            angular_momentum_aux[0] = 0.0;
            angular_momentum_aux[1] = 0.0;
            angular_momentum_aux[2] = 0.0;

            if (Fix_Ang_vel[0] == true || Fix_Ang_vel[1] == true || Fix_Ang_vel[2] == true) {
                double LocalTensor[3][3];
                double GlobalTensor[3][3];
                GeometryFunctions::ConstructLocalTensor(moments_of_inertia, LocalTensor);
                GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensor, GlobalTensor);
                GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensor, angular_velocity, angular_momentum_aux);
            }


            for (int j = 0; j < 3; j++) {
                if (Fix_Ang_vel[j] == false){
                    angular_momentum[j] += moment_reduction_factor * torque[j] * delta_t;
                }
                else {
                    angular_momentum[j] = angular_momentum_aux[j];
                }
            }

            CalculateAngularVelocityRK(Orientation, moments_of_inertia, angular_momentum, angular_velocity, delta_t, Fix_Ang_vel);
            UpdateRotationalVariablesOfCluster(i, moments_of_inertia, rotated_angle, delta_rotation, Orientation, angular_momentum, angular_velocity, delta_t, Fix_Ang_vel);
            GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);

        }//if (StepFlag != 1 && StepFlag != 2)
                
        if (StepFlag == 1 || StepFlag == 2) {
            array_1d<double, 3 > local_angular_acceleration, local_torque, quarter_local_angular_velocity, quarter_angular_velocity, AuxAngularVelocity;
            Quaternion<double  > & AuxOrientation = i.FastGetSolutionStepValue(AUX_ORIENTATION);
            array_1d<double, 3 > & LocalAuxAngularVelocity = i.FastGetSolutionStepValue(LOCAL_AUX_ANGULAR_VELOCITY);

            if (StepFlag == 1) { //PREDICT
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
                    
            if (StepFlag == 2) { //CORRECT
                //Angular velocity and torques are saved in the local framework:
                GeometryFunctions::QuaternionVectorGlobal2Local(AuxOrientation, torque, local_torque);
                CalculateLocalAngularAccelerationByEulerEquations(i,LocalAuxAngularVelocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);
                        
                local_angular_velocity += local_angular_acceleration * delta_t;
                GeometryFunctions::QuaternionVectorLocal2Global(Orientation, local_angular_velocity, angular_velocity);
                        
                GeometryFunctions::QuaternionVectorLocal2Global(AuxOrientation, LocalAuxAngularVelocity, AuxAngularVelocity);
                UpdateRotationalVariables(i, rotated_angle, delta_rotation, AuxAngularVelocity, delta_t, Fix_Ang_vel);
                GeometryFunctions::UpdateOrientation(Orientation, delta_rotation);
            }//if StepFlag == 2
        }//if (StepFlag == 1 || StepFlag == 2)

        KRATOS_CATCH(" ")
    }
}
