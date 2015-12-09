//        
// Author: Miguel AngelCeligueta, maceli@cimne.upc.edu
//


#if !defined(KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED)
#define KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED


// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

namespace Kratos {

    class DEMIntegrationScheme {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        KRATOS_CLASS_POINTER_DEFINITION(DEMIntegrationScheme);

        DEMIntegrationScheme() {
        }

        virtual ~DEMIntegrationScheme() {
        }

        virtual void UpdateLinearDisplacementAndVelocityOfSpheres(ModelPart & rcluster_model_part) { //must be done AFTER the translational motion!

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
        
        virtual void Calculate(ModelPart& model_part) {
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

                if (rCurrentProcessInfo[ROTATION_OPTION] == 0) {
                    UpdateLinearDisplacementAndVelocityOfSpheres(model_part);
                } else {
                    CalculateRotationalMotionOfClusters(model_part);
                }
            }

            KRATOS_CATCH(" ")
        }

        virtual void UpdateTranslationalVariables(
                const ModelPart::NodeIterator& i,
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

        virtual void CalculateTranslationalMotion(ModelPart& model_part, NodesArrayType& pNodes) {

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

                    UpdateTranslationalVariables(i, coor, displ, delta_displ, vel, initial_coor, force, force_reduction_factor, mass, delta_t, Fix_vel);

                } //nodes in the thread
            } //threads


            KRATOS_CATCH(" ")
        }

        virtual void UpdateRotationalVariables(
                const ModelPart::NodeIterator& i,
                array_1d<double, 3 >& rotational_displacement,
                array_1d<double, 3 >& delta_rotation,
                array_1d<double, 3 >& angular_velocity,
                const array_1d<double, 3 >& torque,
                const double moment_reduction_factor,
                const double moment_of_inertia,
                const double delta_t,
                const bool Fix_Ang_vel[3]) {

            KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMIntegrationScheme::UpdateRotationalVariables) shouldn't be accessed, use derived class instead", 0);

        }

        virtual void CalculateRotationalMotion(ModelPart& model_part, NodesArrayType& pNodes) {

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

                    double moment_of_inertia = i->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

                    array_1d<double, 3 >& angular_velocity = i->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                    array_1d<double, 3 >& torque = i->FastGetSolutionStepValue(PARTICLE_MOMENT);
                    array_1d<double, 3 >& rotational_displacement = i->FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
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

                    UpdateRotationalVariables(i, rotational_displacement, delta_rotation, angular_velocity, torque, moment_reduction_factor, moment_of_inertia, delta_t, Fix_Ang_vel);

                    if (if_trihedron_option) {
                        double theta[3] = {0.0};

                        theta[0] = rotational_displacement[0] * 0.5;
                        theta[1] = rotational_displacement[1] * 0.5;
                        theta[2] = rotational_displacement[2] * 0.5;

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

                        array_1d<double, 3>& EulerAngles = i->FastGetSolutionStepValue(EULER_ANGLES);

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

        virtual void CalculateRotationalMotionOfClusters(ModelPart& rcluster_model_part) { //must be done AFTER the translational motion!

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

                        for (int j = 0; j < 3; j++) {
                            //Euler equations in Explicit (Forward Euler) scheme:
                            local_angular_acceleration[j] = (local_torque[j] - (local_angular_velocity[(j + 1) % 3] * moments_of_inertia[(j + 2) % 3] * local_angular_velocity[(j + 2) % 3] - local_angular_velocity[(j + 2) % 3] * moments_of_inertia[(j + 1) % 3] * local_angular_velocity[(j + 1) % 3])) / moments_of_inertia[j];
                            if (if_virtual_mass_option) {
                                local_angular_acceleration[j] = local_angular_acceleration[j] * (1 - coeff);
                            }
                        }

                        //Angular acceleration is saved in the Global framework:
                        GeometryFunctions::VectorLocal2Global(rotation_matrix, local_angular_acceleration, angular_acceleration);

                        bool If_Fix_Rotation[3] = {false, false, false};

                        If_Fix_Rotation[0] = i.Is(DEMFlags::FIXED_ANG_VEL_X);
                        If_Fix_Rotation[1] = i.Is(DEMFlags::FIXED_ANG_VEL_Y);
                        If_Fix_Rotation[2] = i.Is(DEMFlags::FIXED_ANG_VEL_Z);

                        for (int j = 0; j < 3; j++) {
                            if (If_Fix_Rotation[j] == false) {
                                delta_rotation[j] = angular_velocity[j] * delta_t;
                                angular_velocity[j] += angular_acceleration[j] * delta_t;
                                rotated_angle[j] += delta_rotation[j];
                            } else {
                                angular_velocity[j] = 0.0;
                                delta_rotation[j] = 0.0;
                            }
                        }

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

                        cluster_element.UpdatePositionOfSpheres(rotation_matrix, delta_t);
                    } //for Elements
                } //for number of threads
            } //End of parallel region
            
            KRATOS_CATCH(" ")

            }

            virtual std::string Info() const {
                std::stringstream buffer;
                buffer << "DEMIntegrationScheme";
                return buffer.str();
            }

            /// Print information about this object.

            virtual void PrintInfo(std::ostream & rOStream) const {
                rOStream << "DEMIntegrationScheme";
            }

            /// Print object's data.

            virtual void PrintData(std::ostream & rOStream) const {
            }

            protected:

            private:

            DEMIntegrationScheme& operator=(DEMIntegrationScheme const& rOther) {
                return *this;
            }

            /// Copy constructor.

            DEMIntegrationScheme(DEMIntegrationScheme const& rOther) {
                *this = rOther;
            }
        }; // Class DEMIntegrationScheme

        /// input stream function

        inline std::istream& operator>>(std::istream& rIStream, DEMIntegrationScheme& rThis) {
            return rIStream;
        }

        /// output stream function

        inline std::ostream& operator<<(std::ostream& rOStream, const DEMIntegrationScheme& rThis) {
            rThis.PrintInfo(rOStream);
            rOStream << std::endl;
            rThis.PrintData(rOStream);
            return rOStream;
        }

    } // namespace Kratos.

#endif // KRATOS_DEM_INTEGRATION_SCHEME_H_INCLUDED  defined 

