//        
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//

#include "DEM_application.h"
#include "dem_integration_scheme.h"

namespace Kratos {

    DEMIntegrationScheme::DEMIntegrationScheme(){}
    DEMIntegrationScheme::~DEMIntegrationScheme(){}

    void DEMIntegrationScheme::AddSpheresVariables(ModelPart & r_model_part){
        
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES);
        r_model_part.AddNodalSolutionStepVariable(NODAL_MASS);   
        
        r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA); 
        r_model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY); 
        r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT); 
        r_model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE); 
        r_model_part.AddNodalSolutionStepVariable(DELTA_ROTATION);         
    }
    
    void DEMIntegrationScheme::AddClustersVariables(ModelPart & r_model_part){
        
        r_model_part.AddNodalSolutionStepVariable(VELOCITY);
        r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(DELTA_DISPLACEMENT);
        r_model_part.AddNodalSolutionStepVariable(TOTAL_FORCES);
        r_model_part.AddNodalSolutionStepVariable(NODAL_MASS);     
        
        r_model_part.AddNodalSolutionStepVariable(PRINCIPAL_MOMENTS_OF_INERTIA); 
        r_model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY); 
        r_model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT); 
        r_model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE); 
        r_model_part.AddNodalSolutionStepVariable(EULER_ANGLES); 
        r_model_part.AddNodalSolutionStepVariable(DELTA_ROTATION);   
    }

    void DEMIntegrationScheme::UpdateLinearDisplacementAndVelocityOfSpheres(ModelPart & rcluster_model_part) { //must be done AFTER the translational motion!

        KRATOS_TRY

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ElementsArrayType::iterator ElementIterator;

        vector<unsigned int> element_partition;
        ElementsArrayType& pElements = rcluster_model_part.GetCommunicator().LocalMesh().Elements();
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), element_partition);

#pragma omp parallel for
        for (int k = 0; k < (int) OpenMPUtils::GetNumThreads(); k++) {
            ElementIterator i_begin = pElements.ptr_begin() + element_partition[k];
            ElementIterator i_end = pElements.ptr_begin() + element_partition[k + 1];

            for (ElementsArrayType::iterator it = i_begin; it != i_end; ++it) {
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                cluster_element.UpdateLinearDisplacementAndVelocityOfSpheres();
            }
        }
        KRATOS_CATCH(" ")
    }
    
    void DEMIntegrationScheme::Calculate(ModelPart& model_part, int StepFlag) {
        KRATOS_TRY
        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        NodesArrayType& pLocalNodes = model_part.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pGhostNodes = model_part.GetCommunicator().GhostMesh().Nodes();
        CalculateTranslationalMotion(model_part, pLocalNodes, StepFlag);
        CalculateTranslationalMotion(model_part, pGhostNodes, StepFlag);

        if (!rCurrentProcessInfo[CONTAINS_CLUSTERS]) {
            if (rCurrentProcessInfo[ROTATION_OPTION] != 0) {
                CalculateRotationalMotion(model_part, pLocalNodes, StepFlag);
                CalculateRotationalMotion(model_part, pGhostNodes, StepFlag);
            }
        } else {

            if (rCurrentProcessInfo[ROTATION_OPTION] == 0) {
                UpdateLinearDisplacementAndVelocityOfSpheres(model_part);
            } else {
                CalculateRotationalMotionOfClusters(model_part, StepFlag);
            }
        }
        KRATOS_CATCH(" ")
    }
        
    void DEMIntegrationScheme::UpdateTranslationalVariables(
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

        KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateTranslationalVariables) shouldn't be accessed, use derived class instead", 0);
    }
    
    void DEMIntegrationScheme::CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag) {
        KRATOS_TRY
        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        double virtual_mass_coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];
        bool virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
        double force_reduction_factor = 1.0;
        if (virtual_mass_option) {
            force_reduction_factor = 1.0 - virtual_mass_coeff;
            if ((force_reduction_factor > 1.0) || (force_reduction_factor < 0.0)) {
                KRATOS_THROW_ERROR(std::runtime_error, "The virtual mass coefficient is either larger than 1 or negative: virtual_mass_coeff= ", virtual_mass_coeff)
            }
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
                //array_1d<double, 3 >& displ_ant = i.FastGetSolutionStepValue(DISPLACEMENT,2);
                array_1d<double, 3 >& delta_displ = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                /*
                if(i.Id()==1)
                {
                KRATOS_WATCH("  ")
                KRATOS_WATCH(delta_displ)
                KRATOS_WATCH(displ)
                KRATOS_WATCH(displ_ant)
                KRATOS_WATCH(displ-displ_ant)
                }
                */
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
    
    void DEMIntegrationScheme::CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes, int StepFlag) {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        vector<unsigned int> node_partition;
        bool if_trihedron_option = (bool) rCurrentProcessInfo[TRIHEDRON_OPTION];
        double virtual_mass_coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];
        bool virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
        double moment_reduction_factor = 1.0;
        if (virtual_mass_option) {
            moment_reduction_factor = 1.0 - virtual_mass_coeff;
            if ((moment_reduction_factor > 1.0) || (moment_reduction_factor < 0.0)) {
                KRATOS_THROW_ERROR(std::runtime_error, "The virtual mass coefficient is either larger than 1 or negative: virtual_mass_coeff= ", virtual_mass_coeff)
            }
        }      

        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

#pragma omp parallel for
        for (int k = 0; k < (int) number_of_threads; k++) {

            NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k + 1];

            for (ModelPart::NodeIterator i_iterator = i_begin; i_iterator != i_end; ++i_iterator) {
                Node < 3 > & i = *i_iterator;

                double moment_of_inertia = i.FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

                array_1d<double, 3 >& angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 >& torque = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
                array_1d<double, 3 >& rotated_angle = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                array_1d<double, 3 >& delta_rotation = i.FastGetSolutionStepValue(DELTA_ROTATION);
                double Orientation_real;
                array_1d<double, 3 > Orientation_imag;                    

                bool Fix_Ang_vel[3] = {false, false, false};

                Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
                Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
                Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

                array_1d<double, 3 > angular_acceleration;                    
                CalculateLocalAngularAcceleration(i, moment_of_inertia, torque, moment_reduction_factor,angular_acceleration);

                UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);

                if (if_trihedron_option) {
                    double theta[3] = {0.0};

                    theta[0] = rotated_angle[0] * 0.5;
                    theta[1] = rotated_angle[1] * 0.5;
                    theta[2] = rotated_angle[2] * 0.5;

                    double thetaMag = sqrt(theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2]);
                    if (thetaMag * thetaMag * thetaMag * thetaMag / 24.0 < DBL_EPSILON) { //Taylor: low angle                      
                        Orientation_real = 1 + thetaMag * thetaMag / 2;
                        Orientation_imag[0] = theta[0] * (1 - thetaMag * thetaMag / 6);
                        Orientation_imag[1] = theta[1] * (1 - thetaMag * thetaMag / 6);
                        Orientation_imag[2] = theta[2] * (1 - thetaMag * thetaMag / 6);
                    } else {
                        Orientation_real = cos(thetaMag);
                        Orientation_imag[0] = (theta[0] / thetaMag) * sin(thetaMag);
                        Orientation_imag[1] = (theta[1] / thetaMag) * sin(thetaMag);
                        Orientation_imag[2] = (theta[2] / thetaMag) * sin(thetaMag);
                    }

                    array_1d<double, 3>& EulerAngles = i.FastGetSolutionStepValue(EULER_ANGLES);

                    double test = Orientation_imag[0] * Orientation_imag[1] + Orientation_imag[2] * Orientation_real;

                    if (test > 0.49999999) { // singularity at north pole                     
                        EulerAngles[0] = 2 * atan2(Orientation_imag[0], Orientation_real);
                        EulerAngles[1] = KRATOS_M_PI;
                        EulerAngles[2] = 0.0;
                    } else if (test < -0.49999999) { // singularity at south pole                                       
                        EulerAngles[0] = -2 * atan2(Orientation_imag[0], Orientation_real);
                        EulerAngles[1] = -KRATOS_M_PI;
                        EulerAngles[2] = 0.0;
                    } else {
                        EulerAngles[0] = atan2(2 * Orientation_real * Orientation_imag[0] + 2 * Orientation_imag[1] * Orientation_imag[2], 1 - 2 * Orientation_imag[0] * Orientation_imag[0] - 2 * Orientation_imag[1] * Orientation_imag[1]);
                        EulerAngles[1] = asin(2 * Orientation_real * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[0]);
                        EulerAngles[2] = -atan2(2 * Orientation_real * Orientation_imag[2] + 2 * Orientation_imag[0] * Orientation_imag[1], 1 - 2 * Orientation_imag[1] * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[2]);
                    }
                }// Trihedron Option                  
            }//for Node                  
        }//for OMP            

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
    
    void DEMIntegrationScheme::CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part, int StepFlag) { //must be done AFTER the translational motion!

        KRATOS_TRY

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ElementsArrayType::iterator ElementIterator;
        ProcessInfo& rCurrentProcessInfo = rcluster_model_part.GetProcessInfo();

        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        double virtual_mass_coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];
        bool virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
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

#pragma omp parallel
        {
            double rotation_matrix[3][3];

#pragma omp for
            for (int k = 0; k < (int) OpenMPUtils::GetNumThreads(); k++) {
                ElementIterator i_begin = pElements.ptr_begin() + element_partition[k];
                ElementIterator i_end = pElements.ptr_begin() + element_partition[k + 1];

                for (ElementsArrayType::iterator it = i_begin; it != i_end; ++it) {
                    Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                    Node < 3 > & i = cluster_element.GetGeometry()[0];

                    array_1d<double, 3 > & moments_of_inertia = i.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
                    array_1d<double, 3 > & angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    array_1d<double, 3 > & torque = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
                    array_1d<double, 3 > & rotated_angle = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                    array_1d<double, 3 > & EulerAngles = i.FastGetSolutionStepValue(EULER_ANGLES);
                    array_1d<double, 3 > & delta_rotation = i.FastGetSolutionStepValue(DELTA_ROTATION);

                    GeometryFunctions::GetRotationMatrix(EulerAngles, rotation_matrix);

                    array_1d<double, 3 > local_angular_velocity, local_angular_acceleration, local_torque, angular_acceleration;

                    //Angular velocity and torques are saved in the local framework:
                    GeometryFunctions::VectorGlobal2Local(rotation_matrix, torque, local_torque);
                    GeometryFunctions::VectorGlobal2Local(rotation_matrix, angular_velocity, local_angular_velocity);

                    CalculateLocalAngularAccelerationByEulerEquations(i,local_angular_velocity,moments_of_inertia,local_torque, moment_reduction_factor,local_angular_acceleration);                        

                    //Angular acceleration is saved in the Global framework:
                    GeometryFunctions::VectorLocal2Global(rotation_matrix, local_angular_acceleration, angular_acceleration);

                    bool Fix_Ang_vel[3] = {false, false, false};

                    Fix_Ang_vel[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
                    Fix_Ang_vel[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
                    Fix_Ang_vel[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

                    UpdateRotationalVariables(StepFlag, i, rotated_angle, delta_rotation, angular_velocity, angular_acceleration, delta_t, Fix_Ang_vel);                                               

                    double ang = sqrt(delta_rotation[0] * delta_rotation[0] + delta_rotation[1] * delta_rotation[1] + delta_rotation[2] * delta_rotation[2]);
                    if (ang) {

                        array_1d<double, 3 > e1, e2;

                        e1[0] = rotation_matrix[0][0]; e1[1] = rotation_matrix[0][1]; e1[2] = rotation_matrix[0][2];
                        e2[0] = rotation_matrix[1][0]; e2[1] = rotation_matrix[1][1]; e2[2] = rotation_matrix[1][2];

                        array_1d<double, 3 > new_axes1, new_axes2, new_axes3, axis;

                        noalias(axis) = delta_rotation / ang;

                        GeometryFunctions::RotateRightHandedBasisAroundAxis(e1, e2, axis, ang, new_axes1, new_axes2, new_axes3);

                        rotation_matrix[0][0] = new_axes1[0]; rotation_matrix[0][1] = new_axes1[1]; rotation_matrix[0][2] = new_axes1[2]; 
                        rotation_matrix[1][0] = new_axes2[0]; rotation_matrix[1][1] = new_axes2[1]; rotation_matrix[1][2] = new_axes2[2]; 
                        rotation_matrix[2][0] = new_axes3[0]; rotation_matrix[2][1] = new_axes3[1]; rotation_matrix[2][2] = new_axes3[2];

                        GeometryFunctions::GetEulerAngles(rotation_matrix, EulerAngles);
                    } //if ang                                                                                                    
                    cluster_element.UpdatePositionOfSpheres(rotation_matrix);
                } //for Elements
            } //for number of threads
        } //End of parallel region
        KRATOS_CATCH(" ")
    }
}
