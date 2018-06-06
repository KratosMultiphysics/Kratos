#ifndef POST_UTILITIES_H
#define POST_UTILITIES_H

#include "utilities/timer.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/Particle_Contact_Element.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utilities/openmp_utils.h"

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace Kratos {

    class PostUtilities {

    public:

        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::NodesContainerType NodesContainerType;

        KRATOS_CLASS_POINTER_DEFINITION(PostUtilities);

        /// Default constructor.

        PostUtilities() {};

        /// Destructor.

        virtual ~PostUtilities() {};

        void AddModelPartToModelPart(ModelPart& rCompleteModelPart, ModelPart& rModelPartToAdd)
        {
            ////WATCH OUT! This function respects the existing Id's!
            KRATOS_TRY;

            //preallocate the memory needed
            int tot_nodes = rCompleteModelPart.Nodes().size() + rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().size();
            int tot_elements = rCompleteModelPart.Elements().size() + rModelPartToAdd.GetCommunicator().LocalMesh().Elements().size();
            rCompleteModelPart.Nodes().reserve(tot_nodes);
            rCompleteModelPart.Elements().reserve(tot_elements);
            for (ModelPart::NodesContainerType::ptr_iterator node_it = rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_begin(); node_it != rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_end(); node_it++)
            {
                rCompleteModelPart.Nodes().push_back(*node_it);
            }

            for (ModelPart::ElementsContainerType::ptr_iterator elem_it = rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_begin(); elem_it != rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_end(); elem_it++)
            {
                rCompleteModelPart.Elements().push_back(*elem_it);
            }

            KRATOS_CATCH("");
        }

        void AddSpheresNotBelongingToClustersToMixModelPart(ModelPart& rCompleteModelPart, ModelPart& rModelPartToAdd)
        {
            ////WATCH OUT! This function respects the existing Id's!
            KRATOS_TRY;

            //preallocate the memory needed
            int tot_size = rCompleteModelPart.Nodes().size();
            for (ModelPart::NodesContainerType::ptr_iterator node_it = rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_begin(); node_it != rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_end(); node_it++)
            {
                ModelPart::NodeIterator i_iterator = node_it;
                Node < 3 > & i = *i_iterator;
                if (i.IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {tot_size += 1;}
            }
            rCompleteModelPart.Nodes().reserve(tot_size);
            rCompleteModelPart.Elements().reserve(tot_size);
            for (ModelPart::NodesContainerType::ptr_iterator node_it = rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_begin(); node_it != rModelPartToAdd.GetCommunicator().LocalMesh().Nodes().ptr_end(); node_it++)
            {
                ModelPart::NodeIterator i_iterator = node_it;
                Node < 3 > & i = *i_iterator;
                if (i.IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {rCompleteModelPart.Nodes().push_back(*node_it);}
            }

            for (ModelPart::ElementsContainerType::ptr_iterator elem_it = rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_begin(); elem_it != rModelPartToAdd.GetCommunicator().LocalMesh().Elements().ptr_end(); elem_it++)
            {
                Node < 3 >& i = (*elem_it)->GetGeometry()[0];
                if (i.IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {rCompleteModelPart.Elements().push_back(*elem_it);}
            }

            KRATOS_CATCH("");
        }

        array_1d<double,3> VelocityTrap(ModelPart& rModelPart, const array_1d<double,3>& low_point, const array_1d<double,3>& high_point) {

            ElementsArrayType& pElements = rModelPart.GetCommunicator().LocalMesh().Elements();

            OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());

            double velocity_X = 0.0, velocity_Y = 0.0, velocity_Z = 0.0;

            int number_of_elements = 0;

            #pragma omp parallel for reduction(+: velocity_X, velocity_Y, velocity_Z, number_of_elements)

            for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {

                ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
                ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

                for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {

                    array_1d<double,3> coor = (it)->GetGeometry()[0].Coordinates();

                    if (coor[0] >= low_point[0] && coor[0] <= high_point[0] &&
                        coor[1] >= low_point[1] && coor[1] <= high_point[1] &&
                        coor[2] >= low_point[2] && coor[2] <= high_point[2]) {

                        velocity_X += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X);
                        velocity_Y += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y);
                        velocity_Z += (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z);

                        number_of_elements++;
                    }

                } //elements for

                for (int i = 0; i < 3; ++i) {

                    if (high_point[i] < low_point[i]) {
                        KRATOS_THROW_ERROR(std::logic_error, "Check the limits of the Velocity Trap Box. Maximum coordinates smaller than minimum coordinates.", "");
                    }
                }

            } //parallel for

            if (number_of_elements) {
                velocity_X /= number_of_elements;
                velocity_Y /= number_of_elements;
                velocity_Z /= number_of_elements;
            }

            array_1d<double,3> velocity;

            velocity[0] = velocity_X;
            velocity[1] = velocity_Y;
            velocity[2] = velocity_Z;

            return velocity;

        }//VelocityTrap


        void IntegrationOfForces(ModelPart::NodesContainerType& mesh_nodes , array_1d<double, 3>& total_forces,
                                 array_1d<double, 3>& rotation_center, array_1d<double, 3>& total_moment) {

            for (ModelPart::NodesContainerType::ptr_iterator node_pointer_it = mesh_nodes.ptr_begin();
                                                             node_pointer_it != mesh_nodes.ptr_end(); ++node_pointer_it) {

                const array_1d<double, 3>& contact_forces_summed_at_structure_point = (*node_pointer_it)->FastGetSolutionStepValue(CONTACT_FORCES);
                noalias(total_forces) += contact_forces_summed_at_structure_point;

                array_1d<double, 3> vector_from_structure_center_to_structure_point;
                noalias(vector_from_structure_center_to_structure_point) = (*node_pointer_it)->Coordinates() - rotation_center;

                array_1d<double, 3> moment_to_add;
                GeometryFunctions::CrossProduct(vector_from_structure_center_to_structure_point, contact_forces_summed_at_structure_point, moment_to_add);

                noalias(total_moment) += moment_to_add;
            }
        }

        void IntegrationOfElasticForces(ModelPart::NodesContainerType& mesh_nodes, array_1d<double, 3>& total_forces) {

            for (ModelPart::NodesContainerType::ptr_iterator node_pointer_it = mesh_nodes.ptr_begin(); node_pointer_it != mesh_nodes.ptr_end(); ++node_pointer_it) {

                const array_1d<double, 3> elastic_forces_added_up_at_node = (*node_pointer_it)->FastGetSolutionStepValue(ELASTIC_FORCES);
                noalias(total_forces) += elastic_forces_added_up_at_node;
            }
        }

        array_1d<double, 3> ComputePoisson(ModelPart& rModelPart) {

            ElementsArrayType& pElements = rModelPart.GetCommunicator().LocalMesh().Elements();
            double total_poisson_value   = 0.0;
            unsigned int number_of_spheres_to_evaluate_poisson = 0;
            array_1d<double, 3> return_data = ZeroVector(3);

            // TODO: Add OpenMP code
            for (unsigned int k = 0; k < pElements.size(); k++) {

                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                Element* raw_p_element = &(*it);
                SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(raw_p_element);
                double& particle_poisson_value = p_sphere->GetGeometry()[0].FastGetSolutionStepValue(POISSON_VALUE);
                particle_poisson_value         = 0.0;
                double epsilon_XY              = 0.0;
                double epsilon_Z               = 0.0;
                unsigned int number_of_neighbors_per_sphere_to_evaluate_poisson = 0;
                array_1d<double, 3> other_to_me_vector;
                array_1d<double, 3> initial_other_to_me_vector;
                unsigned int number_of_neighbors = p_sphere->mNeighbourElements.size();

                for (unsigned int i = 0; i < number_of_neighbors; i++) {

                    if (p_sphere->mNeighbourElements[i] == NULL) continue;
                    noalias(other_to_me_vector)         = p_sphere->GetGeometry()[0].Coordinates() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                    noalias(initial_other_to_me_vector) = p_sphere->GetGeometry()[0].GetInitialPosition() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].GetInitialPosition();
                    double initial_distance_XY = sqrt(initial_other_to_me_vector[0] * initial_other_to_me_vector[0] + initial_other_to_me_vector[1] * initial_other_to_me_vector[1]);
                    double initial_distance_Z  = initial_other_to_me_vector[2];

                    if (initial_distance_XY && initial_distance_Z) {
                        epsilon_XY = -1 + sqrt(other_to_me_vector[0] * other_to_me_vector[0] + other_to_me_vector[1] * other_to_me_vector[1]) / initial_distance_XY;
                        epsilon_Z  = -1 + fabs(other_to_me_vector[2] / initial_distance_Z);
                    } else continue;

                    if (epsilon_Z) { // Should it be added here 'if p_sphere->Id() < p_sphere->mNeighbourElements[i]->Id()'?
                        if (((-epsilon_XY / epsilon_Z) > 0.5) || ((-epsilon_XY / epsilon_Z) < 0.0)) continue; // TODO: Check this
                        particle_poisson_value -= epsilon_XY / epsilon_Z;
                        number_of_neighbors_per_sphere_to_evaluate_poisson++;
                    } else continue;
                }

                if (number_of_neighbors_per_sphere_to_evaluate_poisson) {
                    particle_poisson_value /= number_of_neighbors_per_sphere_to_evaluate_poisson;
                    number_of_spheres_to_evaluate_poisson++;
                    total_poisson_value += particle_poisson_value;
                }
            }

            if (number_of_spheres_to_evaluate_poisson) total_poisson_value /= number_of_spheres_to_evaluate_poisson;

            return_data[0] = total_poisson_value;
            return return_data;

        } //ComputePoisson


        array_1d<double, 3> ComputePoisson2D(ModelPart& rModelPart) { // TODO: Adjust this function to the new changes made in the 3D version

            ElementsArrayType& pElements = rModelPart.GetCommunicator().LocalMesh().Elements();
            double total_poisson_value     = 0.0;
            unsigned int number_of_bonds_to_evaluate_poisson = 0;
            array_1d<double, 3> return_data = ZeroVector(3);
            double total_epsilon_y_value = 0.0;

            // TODO: Add OpenMP code
            for (unsigned int k = 0; k < pElements.size(); k++) {

                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                Element* raw_p_element = &(*it);
                SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(raw_p_element);
                double& particle_poisson_value = p_sphere->GetGeometry()[0].FastGetSolutionStepValue(POISSON_VALUE);
                particle_poisson_value         = 0.0;
                double epsilon_X               = 0.0;
                double epsilon_Y               = 0.0;
                unsigned int number_of_neighbors_to_evaluate_poisson = 0;
                array_1d<double, 3> other_to_me_vector;
                array_1d<double, 3> initial_other_to_me_vector;
                double average_sphere_epsilon_y_value = 0.0;

                unsigned int number_of_neighbors = p_sphere->mNeighbourElements.size();

                for (unsigned int i = 0; i < number_of_neighbors; i++)
                {
                    if (p_sphere->mNeighbourElements[i] == NULL) continue;
                    noalias(other_to_me_vector)         = p_sphere->GetGeometry()[0].Coordinates() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                    noalias(initial_other_to_me_vector) = p_sphere->GetGeometry()[0].GetInitialPosition() - p_sphere->mNeighbourElements[i]->GetGeometry()[0].GetInitialPosition();
                    double initial_distance_X = initial_other_to_me_vector[0];
                    double initial_distance_Y = initial_other_to_me_vector[1];
                    if (initial_distance_X && initial_distance_Y) {
                        epsilon_X = -1 + fabs(other_to_me_vector[0] / initial_distance_X);
                        epsilon_Y = -1 + fabs(other_to_me_vector[1] / initial_distance_Y);
                    }
                    if (epsilon_Y) {
                        particle_poisson_value -= epsilon_X / epsilon_Y;
                        number_of_neighbors_to_evaluate_poisson++;
                            total_poisson_value -= epsilon_X / epsilon_Y;
                            number_of_bonds_to_evaluate_poisson++;
                    }
                    average_sphere_epsilon_y_value += epsilon_Y;
                }
                if (number_of_neighbors_to_evaluate_poisson) particle_poisson_value /= number_of_neighbors_to_evaluate_poisson;

                total_epsilon_y_value += average_sphere_epsilon_y_value / number_of_neighbors;
            }

            if (number_of_bonds_to_evaluate_poisson) total_poisson_value /= number_of_bonds_to_evaluate_poisson;

            total_epsilon_y_value /= pElements.size();

            return_data[0] = total_poisson_value;
            return_data[1] = total_epsilon_y_value;

            return return_data;

        } //ComputePoisson2D


        void ComputeEulerAngles(ModelPart& rSpheresModelPart, ModelPart& rClusterModelPart) {

            ProcessInfo& r_process_info = rSpheresModelPart.GetProcessInfo();
            bool if_trihedron_option = (bool) r_process_info[TRIHEDRON_OPTION];

            typedef ModelPart::NodesContainerType NodesArrayType;
            NodesArrayType& pSpheresNodes = rSpheresModelPart.GetCommunicator().LocalMesh().Nodes();
            NodesArrayType& pClusterNodes = rClusterModelPart.GetCommunicator().LocalMesh().Nodes();

            #pragma omp parallel for
            for (int k = 0; k < (int) pSpheresNodes.size(); k++) {

                ModelPart::NodeIterator i_iterator = pSpheresNodes.ptr_begin() + k;
                Node < 3 > & i = *i_iterator;

                array_1d<double, 3 >& rotated_angle = i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);

                if (if_trihedron_option && i.IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                    array_1d<double, 3 >& EulerAngles = i.FastGetSolutionStepValue(EULER_ANGLES);
                    GeometryFunctions::EulerAnglesFromRotationAngle(EulerAngles, rotated_angle);
                } // if_trihedron_option && Not BELONGS_TO_A_CLUSTER
            }//for Node

            #pragma omp parallel for
            for (int k = 0; k < (int) pClusterNodes.size(); k++) {

                ModelPart::NodeIterator i_iterator = pClusterNodes.ptr_begin() + k;
                Node < 3 > & i = *i_iterator;

                Quaternion<double>& Orientation = i.FastGetSolutionStepValue(ORIENTATION);
                array_1d<double, 3 >& EulerAngles = i.FastGetSolutionStepValue(EULER_ANGLES);
                Orientation.ToEulerAngles(EulerAngles);

                double& OrientationReal = i.FastGetSolutionStepValue(ORIENTATION_REAL);
                OrientationReal = Orientation.w();

                array_1d<double, 3 >& OrientationImag = i.FastGetSolutionStepValue(ORIENTATION_IMAG);
                OrientationImag[0] = Orientation.x();
                OrientationImag[1] = Orientation.y();
                OrientationImag[2] = Orientation.z();
            }//for Node
        } //ComputeEulerAngles


        double QuasiStaticAdimensionalNumber(ModelPart& rParticlesModelPart, ModelPart& rContactModelPart, ProcessInfo& r_process_info) {

            double adimensional_value = 0.0;

            ElementsArrayType& pParticleElements = rParticlesModelPart.GetCommunicator().LocalMesh().Elements();

            OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pParticleElements.size(), this->GetElementPartition());

            array_1d<double,3> particle_forces;

            const array_1d<double,3>& gravity = r_process_info[GRAVITY];

            double total_force = 0.0;

            //#pragma omp parallel for
            #pragma omp parallel for reduction(+:total_force)
            for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {
                ElementsArrayType::iterator it_begin = pParticleElements.ptr_begin() + this->GetElementPartition()[k];
                ElementsArrayType::iterator it_end   = pParticleElements.ptr_begin() + this->GetElementPartition()[k + 1];

                for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {

                    Element::GeometryType& geom = it->GetGeometry();

                    if (geom[0].IsNot(DEMFlags::FIXED_VEL_X) && geom[0].IsNot(DEMFlags::FIXED_VEL_Y) && geom[0].IsNot(DEMFlags::FIXED_VEL_Z))
                    {
                        particle_forces  = geom[0].FastGetSolutionStepValue(TOTAL_FORCES);
                        double mass = geom[0].FastGetSolutionStepValue(NODAL_MASS);


                        particle_forces[0] += mass * gravity[0];
                        particle_forces[1] += mass * gravity[1];
                        particle_forces[2] += mass * gravity[2];

                        double module = 0.0;
                        GeometryFunctions::module(particle_forces, module);

                        total_force += module;

                    } //if

                }//balls

            }//paralel

            ElementsArrayType& pContactElements        = rContactModelPart.GetCommunicator().LocalMesh().Elements();

            OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pContactElements.size(), this->GetElementPartition());

            array_1d<double,3> contact_forces;
            double total_elastic_force = 0.0;

            #pragma omp parallel for reduction(+:total_elastic_force)
            for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {
                ElementsArrayType::iterator it_begin = pContactElements.ptr_begin() + this->GetElementPartition()[k];
                ElementsArrayType::iterator it_end   = pContactElements.ptr_begin() + this->GetElementPartition()[k + 1];

                for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){

                    Element::GeometryType& geom = it->GetGeometry();

                    if (geom[0].IsNot(DEMFlags::FIXED_VEL_X) && geom[0].IsNot(DEMFlags::FIXED_VEL_Y) && geom[0].IsNot(DEMFlags::FIXED_VEL_Z) &&
                        geom[1].IsNot(DEMFlags::FIXED_VEL_X) && geom[1].IsNot(DEMFlags::FIXED_VEL_Y) && geom[1].IsNot(DEMFlags::FIXED_VEL_Z)) {

                        contact_forces  = it->GetValue(LOCAL_CONTACT_FORCE);
                        double module = 0.0;
                        GeometryFunctions::module(contact_forces, module);
                        total_elastic_force += module;
                    }
                }
            }

            if (total_elastic_force != 0.0)
            {
                adimensional_value =  total_force/total_elastic_force;
            }
            else
            {
                KRATOS_THROW_ERROR(std::runtime_error,"There are no elastic forces= ", total_elastic_force)
            }

            return adimensional_value;

        }//QuasiStaticAdimensionalNumber


        DenseVector<unsigned int>& GetElementPartition() {return (mElementPartition);};

    protected:

        DenseVector<unsigned int> mElementPartition;

    }; // Class PostUtilities

} // namespace Kratos.

#endif // POST_UTILITIES_H
