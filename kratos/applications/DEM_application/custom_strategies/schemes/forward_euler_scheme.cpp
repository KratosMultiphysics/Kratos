// External includes

// Project includes

#include "DEM_application.h"
#include "forward_euler_scheme.h"

namespace Kratos {

    void ForwardEulerScheme::Calculate(ModelPart& model_part) {

        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        NodesArrayType& pLocalNodes = model_part.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pGhostNodes = model_part.GetCommunicator().GhostMesh().Nodes();
        CalculateTranslationalMotion(model_part, pLocalNodes);
        CalculateTranslationalMotion(model_part, pGhostNodes);

        if (!rCurrentProcessInfo[CONTAINS_CLUSTERS]) {
            if (rCurrentProcessInfo[ROTATION_OPTION] != 0) {
                CalculateRotationalMotion(model_part, pLocalNodes);
                CalculateRotationalMotion(model_part, pGhostNodes);
            }
        } else {
            
            if (rCurrentProcessInfo[ROTATION_OPTION] == 0){
                UpdateLinearDisplacementAndVelocityOfSpheres(model_part);
            }
            else {
                CalculateRotationalMotionOfClusters(model_part);
            }
        }

        KRATOS_CATCH(" ")
    }

    void ForwardEulerScheme::CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes) {
        
        KRATOS_TRY 
        
        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        double virtual_mass_coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];
        bool if_virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
        vector<unsigned int> node_partition;

        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for shared(delta_t)
        for (int k = 0; k < (int) number_of_threads; k++) {
            NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k + 1];

            for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {
                if (i->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
                array_1d<double, 3 >& vel = i->FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3 >& displ = i->FastGetSolutionStepValue(DISPLACEMENT);
                array_1d<double, 3 >& delta_displ = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double, 3 >& coor = i->Coordinates();
                array_1d<double, 3 >& initial_coor = i->GetInitialPosition();
                array_1d<double, 3 >& force = i->FastGetSolutionStepValue(TOTAL_FORCES);

                double mass = i->FastGetSolutionStepValue(NODAL_MASS);

                double force_reduction_factor = 1.0;

                if (if_virtual_mass_option) {
                    force_reduction_factor = 1.0 - virtual_mass_coeff;                    
                    if (virtual_mass_coeff < 0.0) KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= ", virtual_mass_coeff)
                }
                
                bool Fix_vel[3] = {false, false, false};

                Fix_vel[0] = i->Is(DEMFlags::FIXED_VEL_X);
                Fix_vel[1] = i->Is(DEMFlags::FIXED_VEL_Y);
                Fix_vel[2] = i->Is(DEMFlags::FIXED_VEL_Z);
                
                for (int k = 0; k < 3; k++) {
                    if (Fix_vel[k] == false) {
                        delta_displ[k] = delta_t * vel[k];
                        displ[k] += delta_displ[k];
                        coor[k] = initial_coor[k] + displ[k];
                        vel[k] += delta_t * force_reduction_factor * force[k] / mass;                          
                    }
                    else {
                        delta_displ[k] = delta_t * vel[k];
                        displ[k] += delta_displ[k];
                        coor[k] = initial_coor[k] + displ[k];
                    }
                } // dimensions  
                
            } //nodes in the thread
        } //threads
        
        KRATOS_CATCH(" ")
    }

    void ForwardEulerScheme::CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes) {
        
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = model_part.GetProcessInfo();
        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        bool if_virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
        vector<unsigned int> node_partition;
        bool if_trihedron_option = (bool) rCurrentProcessInfo[TRIHEDRON_OPTION];
        double virtual_mass_coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];

        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        #pragma omp parallel for
        for (int k = 0; k < (int) number_of_threads; k++) {

            NodesArrayType::iterator i_begin = pNodes.ptr_begin() + node_partition[k];
            NodesArrayType::iterator i_end = pNodes.ptr_begin() + node_partition[k + 1];

            for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {

                double PMomentOfInertia = i->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

                array_1d<double, 3 >& AngularVel = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 >& RotaMoment = i->FastGetSolutionStepValue(PARTICLE_MOMENT);
                array_1d<double, 3 >& Rota_Displace = i->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                array_1d<double, 3 >& delta_rotation = i->FastGetSolutionStepValue(DELTA_ROTATION);
                double Orientation_real;
                array_1d<double, 3 > Orientation_imag;                
                
                double moment_reduction_factor = 1.0;

                if (if_virtual_mass_option) {
                    moment_reduction_factor = 1.0 - virtual_mass_coeff;                    
                    if (virtual_mass_coeff < 0.0) KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= ", virtual_mass_coeff)
                }
                   
                bool Fix_Ang_vel[3] = {false, false, false};

                Fix_Ang_vel[0] = i->Is(DEMFlags::FIXED_ANG_VEL_X);
                Fix_Ang_vel[1] = i->Is(DEMFlags::FIXED_ANG_VEL_Y);
                Fix_Ang_vel[2] = i->Is(DEMFlags::FIXED_ANG_VEL_Z);

                for (int k = 0; k < 3; k++) {
                    if (Fix_Ang_vel[k] == false) {
                        delta_rotation[k] = AngularVel[k] * delta_t;
                        Rota_Displace[k] += delta_rotation[k];
                        AngularVel[k] += delta_t * moment_reduction_factor * RotaMoment[k] / PMomentOfInertia;                        
                    }
                    else {
                        delta_rotation[k] = AngularVel[k] * delta_t;
                        Rota_Displace[k] += delta_rotation[k];
                    }
                }

                if (if_trihedron_option) {
                    double theta[3] = {0.0};

                    theta[0] = Rota_Displace[0] * 0.5;
                    theta[1] = Rota_Displace[1] * 0.5;
                    theta[2] = Rota_Displace[2] * 0.5;

                    double thetaMag = sqrt(theta[0] * theta[0] + theta[1] * theta[1] + theta[2] * theta[2]);
                    if (thetaMag * thetaMag * thetaMag * thetaMag / 24.0 < DBL_EPSILON) { //Taylor: low angle                      
                        Orientation_real = 1 + thetaMag * thetaMag / 2;
                        Orientation_imag[0] = theta[0] * (1 - thetaMag * thetaMag / 6);
                        Orientation_imag[1] = theta[1] * (1 - thetaMag * thetaMag / 6);
                        Orientation_imag[2] = theta[2] * (1 - thetaMag * thetaMag / 6);
                    }
                    else {
                        Orientation_real = cos(thetaMag);
                        Orientation_imag[0] = (theta[0] / thetaMag) * sin(thetaMag);
                        Orientation_imag[1] = (theta[1] / thetaMag) * sin(thetaMag);
                        Orientation_imag[2] = (theta[2] / thetaMag) * sin(thetaMag);
                    }

                    array_1d<double, 3>& EulerAngles = i->FastGetSolutionStepValue(EULER_ANGLES);

                    double test = Orientation_imag[0] * Orientation_imag[1] + Orientation_imag[2] * Orientation_real;

                    if (test > 0.49999999) { // singularity at north pole                     
                        EulerAngles[0] = 2 * atan2(Orientation_imag[0], Orientation_real);
                        EulerAngles[1] = KRATOS_M_PI;
                        EulerAngles[2] = 0.0;
                    }
                    else if (test < -0.49999999) { // singularity at south pole                                       
                        EulerAngles[0] = -2 * atan2(Orientation_imag[0], Orientation_real);
                        EulerAngles[1] = -KRATOS_M_PI;
                        EulerAngles[2] = 0.0;
                    }
                    else {
                        EulerAngles[0] =  atan2(2 * Orientation_real * Orientation_imag[0] + 2 * Orientation_imag[1] * Orientation_imag[2], 1 - 2 * Orientation_imag[0] * Orientation_imag[0] - 2 * Orientation_imag[1] * Orientation_imag[1]);
                        EulerAngles[1] =  asin( 2 * Orientation_real * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[0]);
                        EulerAngles[2] = -atan2(2 * Orientation_real * Orientation_imag[2] + 2 * Orientation_imag[0] * Orientation_imag[1], 1 - 2 * Orientation_imag[1] * Orientation_imag[1] - 2 * Orientation_imag[2] * Orientation_imag[2]);
                    }
                }// Trihedron Option                  
            }//for Node                  
        }//for OMP            

        KRATOS_CATCH(" ")

    }//rotational_motion                              

    void ForwardEulerScheme::CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part) { //must be done AFTER the translational motion!

        KRATOS_TRY

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ElementsArrayType::iterator ElementIterator;
        ProcessInfo& rCurrentProcessInfo = rcluster_model_part.GetProcessInfo();

        double delta_t = rCurrentProcessInfo[DELTA_TIME];
        bool if_virtual_mass_option = (bool) rCurrentProcessInfo[VIRTUAL_MASS_OPTION];
        double coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];

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

                array_1d<double, 3 > & PMomentsOfInertia = i.FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
                array_1d<double, 3 > & AngularVel = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3 > & RotaMoment = i.FastGetSolutionStepValue(PARTICLE_MOMENT);
                array_1d<double, 3 > & Rota_Displace = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
                array_1d<double, 3 > & EulerAngles = i.FastGetSolutionStepValue(EULER_ANGLES);
                array_1d<double, 3 > & delta_rotation = i.FastGetSolutionStepValue(DELTA_ROTATION);

                GeometryFunctions::GetRotationMatrix(EulerAngles, rotation_matrix);

                array_1d<double, 3 > LocalAngularVel;
                array_1d<double, 3 > LocalRotaAcc;
                array_1d<double, 3 > LocalRotaMoment;
                array_1d<double, 3 > GlobalRotaAcc;

                //Angular velocity and torques are saved in the local framework:
                GeometryFunctions::VectorGlobal2Local(rotation_matrix, RotaMoment, LocalRotaMoment);
                GeometryFunctions::VectorGlobal2Local(rotation_matrix, AngularVel, LocalAngularVel);

                for (int j = 0; j < 3; j++) {
                    //Euler equations in Explicit (Forward Euler) scheme:
                    LocalRotaAcc[j] = (LocalRotaMoment[j] - (LocalAngularVel[(j + 1) % 3] * PMomentsOfInertia[(j + 2) % 3] * LocalAngularVel[(j + 2) % 3] - LocalAngularVel[(j + 2) % 3] * PMomentsOfInertia[(j + 1) % 3] * LocalAngularVel[(j + 1) % 3])) / PMomentsOfInertia[j];
                    if (if_virtual_mass_option) {
                        LocalRotaAcc[j] = LocalRotaAcc[j] * (1 - coeff);
                    }
                }

                //Angular acceleration is saved in the Global framework:
                GeometryFunctions::VectorLocal2Global(rotation_matrix, LocalRotaAcc, GlobalRotaAcc);

                bool If_Fix_Rotation[3] = {false, false, false};

                If_Fix_Rotation[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
                If_Fix_Rotation[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
                If_Fix_Rotation[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

                for (int j = 0; j < 3; j++) {
                    if (If_Fix_Rotation[j] == false) {
                        delta_rotation[j] = AngularVel[j] * delta_t;
                        AngularVel[j] += GlobalRotaAcc[j] * delta_t;                         
                        Rota_Displace[j] += delta_rotation[j];
                    } else {
                        AngularVel[j] = 0.0;
                        delta_rotation[j] = 0.0;
                    }
                }

                double ang = sqrt(delta_rotation[0] * delta_rotation[0] + delta_rotation[1] * delta_rotation[1] + delta_rotation[2] * delta_rotation[2]);
                if (ang) {

                    array_1d<double, 3 > e1;
                    array_1d<double, 3 > e2;

                    e1[0] = rotation_matrix[0][0];
                    e1[1] = rotation_matrix[0][1];
                    e1[2] = rotation_matrix[0][2];

                    e2[0] = rotation_matrix[1][0];
                    e2[1] = rotation_matrix[1][1];
                    e2[2] = rotation_matrix[1][2];

                    array_1d<double, 3 > new_axes1;
                    array_1d<double, 3 > new_axes2;
                    array_1d<double, 3 > new_axes3;
                    array_1d<double, 3 > axis;

                    axis[0] = delta_rotation[0] / ang;
                    axis[1] = delta_rotation[1] / ang;
                    axis[2] = delta_rotation[2] / ang;

                    double cang = cos(ang);
                    double sang = sin(ang);

                    new_axes1[0] = axis[0] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[0] * cang + (-axis[2] * e1[1] + axis[1] * e1[2]) * sang;
                    new_axes1[1] = axis[1] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[1] * cang + (axis[2] * e1[0] - axis[0] * e1[2]) * sang;
                    new_axes1[2] = axis[2] * (axis[0] * e1[0] + axis[1] * e1[1] + axis[2] * e1[2]) * (1 - cang) + e1[2] * cang + (-axis[1] * e1[0] + axis[0] * e1[1]) * sang;

                    new_axes2[0] = axis[0] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[0] * cang + (-axis[2] * e2[1] + axis[1] * e2[2]) * sang;
                    new_axes2[1] = axis[1] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[1] * cang + (axis[2] * e2[0] - axis[0] * e2[2]) * sang;
                    new_axes2[2] = axis[2] * (axis[0] * e2[0] + axis[1] * e2[1] + axis[2] * e2[2]) * (1 - cang) + e2[2] * cang + (-axis[1] * e2[0] + axis[0] * e2[1]) * sang;

                    GeometryFunctions::CrossProduct(new_axes1, new_axes2, new_axes3);

                    rotation_matrix[0][0] = new_axes1[0];
                    rotation_matrix[0][1] = new_axes1[1];
                    rotation_matrix[0][2] = new_axes1[2];

                    rotation_matrix[1][0] = new_axes2[0];
                    rotation_matrix[1][1] = new_axes2[1];
                    rotation_matrix[1][2] = new_axes2[2];

                    rotation_matrix[2][0] = new_axes3[0];
                    rotation_matrix[2][1] = new_axes3[1];
                    rotation_matrix[2][2] = new_axes3[2];

                    GeometryFunctions::GetEulerAngles(rotation_matrix, EulerAngles);
                } //if ang                                                                                                    

                //GeometryFunctions::GetRotationMatrix(EulerAngles, rotation_matrix); //we get the new rotation matrix after having updated the Euler angles
                cluster_element.UpdatePositionOfSpheres(rotation_matrix, delta_t);
            } //for Elements
        } //for number of threads
        } //End of parallel region

        KRATOS_CATCH(" ")
    }//rotational_motion  
    
    
    void ForwardEulerScheme::UpdateLinearDisplacementAndVelocityOfSpheres(ModelPart& rcluster_model_part) { //must be done AFTER the translational motion!

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


}
