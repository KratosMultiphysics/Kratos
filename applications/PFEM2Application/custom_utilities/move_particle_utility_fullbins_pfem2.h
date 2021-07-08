/*
==============================================================================
KratosIncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2011-09-21 12:30:32 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_MOVE_PARTICLE_UTILITY_FULLBINS_PFEM2_INCLUDED)
#define KRATOS_MOVE_PARTICLE_UTILITY_FLUID_FULLBINS_PFEM2_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"

///
#include "includes/dof.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/deprecated_variables.h"
#include "includes/global_pointer_variables.h"
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "processes/node_erase_process.h"
///

#include "utilities/geometry_utilities.h"

#include "includes/model_part.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/cell.h"
#include "utilities/binbased_fast_point_locator.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/point.h"

#include "pfem_2_application_variables.h"
#include "pfem_particle_fluidonly.h"

//#include "utilities/enrich_2d_2dofs.h"
#include "utilities/enrichment_utilities.h"
#include "utilities/openmp_utils.h"

#include "time.h"

//#include "processes/process.h"

namespace Kratos
{
//this class is to be modified by the user to customize the interpolation process
template <unsigned int TDim>
class MoveParticleUtilityFullBinsPFEM2
{
public:
    typedef SpatialContainersConfigure<TDim> Configure;
    typedef SpatialContainersConfigure<TDim> ConfigureType;
    typedef typename ConfigureType::EntityType EntityType;
    typedef typename Configure::PointType PointType;
    /// The definition of the node
    typedef Node<3> NodeType;
    /// The definition of the geometry
    typedef Geometry<NodeType> GeometryType;
    //typedef PointType::CoordinatesArrayType           CoordinatesArrayType;
    typedef typename Configure::ContainerType ContainerType;
    //typedef Configure::PointerType                    PointerType;
    typedef typename Configure::IteratorType IteratorType;
    typedef typename Configure::ResultContainerType ResultContainerType;
    //typedef Configure::ResultPointerType              ResultPointerType;
    typedef typename Configure::ResultIteratorType ResultIteratorType;
    typedef PointerVector<PFEM_Particle_Fluid, PFEM_Particle_Fluid *, std::vector<PFEM_Particle_Fluid *>> ParticlePointerVector;
    //typedef Configure::ContactPairType                ContactPairType;
    //typedef Configure::ContainerContactType           ContainerContactType;
    //typedef Configure::IteratorContactType            IteratorContactType;
    //typedef Configure::PointerContactType             PointerContactType;
    //typedef Configure::PointerTypeIterator            PointerTypeIterator;

    KRATOS_CLASS_POINTER_DEFINITION(MoveParticleUtilityFullBinsPFEM2);

    //template<unsigned int TDim>
    MoveParticleUtilityFullBinsPFEM2(ModelPart &model_part, int maximum_number_of_particles)
        : mr_model_part(model_part), mmaximum_number_of_particles(maximum_number_of_particles)
    {
        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Initializing Utility" << std::endl;

        Check();

        //tools to move the domain, in case we are using a moving domain approach.
        mintialized_transfer_tool = false;
        mcalculation_domain_complete_displacement = ZeroVector(3);
        mcalculation_domain_added_displacement = ZeroVector(3);

        //storing water and air density and their inverses, just in case it is needed for the streamline integration
        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR];
        mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER];


        //loop in elements to change their ID to their position in the array. Easier to get information later.
        //DO NOT PARALELIZE THIS! IT MUST BE SERIAL!!!!!!!!!!!!!!!!!!!!!!
        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        for (unsigned int ii = 0; ii < mr_model_part.Elements().size(); ii++)
        {
            ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
            ielem->SetId(ii + 1);
        }
        mlast_elem_id = (mr_model_part.ElementsEnd() - 1)->Id();
        int node_id = 0;
        // we look for the smallest edge. could be used as a weighting function when going lagrangian->eulerian instead of traditional shape functions(method currently used)
        ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator pnode = inodebegin + ii;
                array_1d<double, 3> position_node;
                double distance = 0.0;
                position_node = pnode->Coordinates();
                GlobalPointersVector<Node<3>> &rneigh = pnode->GetValue(NEIGHBOUR_NODES);
                //we loop all the nodes to check all the edges
                const double number_of_neighbours = double(rneigh.size());
                for (GlobalPointersVector<Node<3>>::iterator inode = rneigh.begin(); inode != rneigh.end(); inode++)
                {
                    array_1d<double, 3> position_difference;
                    position_difference = inode->Coordinates() - position_node;
                    double current_distance = sqrt(pow(position_difference[0], 2) + pow(position_difference[1], 2) + pow(position_difference[2], 2));
                    //if (current_distance>distance)
                    //  distance=current_distance;
                    distance += current_distance / number_of_neighbours;
                }
                //and we save the largest edge.
                pnode->FastGetSolutionStepValue(MEAN_SIZE) = distance;
            }
        }
        ModelPart::NodesContainerType::iterator lastnode = (mr_model_part.NodesEnd() - 1);
        mlast_node_id = lastnode->GetId();

        //we also calculate the element mean size in the same way, for the courant number
        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

        //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
            {
                ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;

                double elem_size;
                array_1d<double, 3> Edge(3, 0.0);
                Edge = ielem->GetGeometry()[1].Coordinates() - ielem->GetGeometry()[0].Coordinates();
                elem_size = Edge[0] * Edge[0];
                for (unsigned int d = 1; d < TDim; d++)
                    elem_size += Edge[d] * Edge[d];

                for (unsigned int i = 2; i < (TDim + 1); i++)
                    for (unsigned int j = 0; j < i; j++)
                    {
                        Edge = ielem->GetGeometry()[i].Coordinates() - ielem->GetGeometry()[j].Coordinates();
                        double Length = Edge[0] * Edge[0];
                        for (unsigned int d = 1; d < TDim; d++)
                            Length += Edge[d] * Edge[d];
                        if (Length < elem_size)
                            elem_size = Length;
                    }
                elem_size = sqrt(elem_size);
                ielem->SetValue(MEAN_SIZE, elem_size);
            }
        }

        //matrix containing the position of the 4/15/45 particles that we will seed at the beggining

        //             BoundedMatrix<double, 2*(1+TDim), 3 > pos;
        //             BoundedMatrix<double, 2*(1+TDim), (1+TDim) > N;

        BoundedMatrix<double, 5 * (1 + TDim), 3> pos;
        BoundedMatrix<double, 5 * (1 + TDim), (1 + TDim)> N;

        //             BoundedMatrix<double, 5*(1+TDim)*3, 3 > pos;
        //             BoundedMatrix<double, 5*(1+TDim)*3, (1+TDim) > N;

        int particle_id = 0;
        mnelems = mr_model_part.Elements().size();

        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "About to resize vectors" << std::endl;

        //setting the right size to the vector containing the particles assigned to each element
        //particles vector. this vector contains ALL the particles in the simulation.
        mparticles_vector.resize(mnelems * mmaximum_number_of_particles);

        //and this vector contains the current number of particles that are in each element (currently zero)
        // mnumber_of_particles_in_elems.resize(mnelems);
        // mnumber_of_particles_in_elems = ZeroVector(mnelems);

        //when moving the particles, an auxiliary vector is necessary (to store the previous number) - perhaps delete
        // mnumber_of_particles_in_elems_aux.resize(mnelems); - perhaps delete

        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "About to create particles" << std::endl;
        //now we seed: LOOP IN ELEMENTS

        for (unsigned int ii = 0; ii < mr_model_part.Elements().size(); ii++)
        {
            ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
            int &number_of_particles = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);
            number_of_particles = 0;
            //int & number_of_water_particles = ielem->GetValue(NUMBER_OF_WATER_PARTICLES);

            Geometry<Node<3>> &geom = ielem->GetGeometry();
            //mareas_vector[i_int]=CalculateArea(geom); UNUSED SO COMMENTED
            ComputeGaussPointPositions_initial(geom, pos, N); //we also have the standard (4), and 45
                                                              //                 ComputeGaussPointPositions_45(geom, pos, N); //we also have the standard (4), and 45
            //now we seed the particles in the current element
            for (unsigned int j = 0; j < pos.size1(); j++)
            {
                ++particle_id;

                PFEM_Particle_Fluid &pparticle = mparticles_vector[particle_id - 1];
                pparticle.X() = pos(j, 0);
                pparticle.Y() = pos(j, 1);
                pparticle.Z() = pos(j, 2);

                pparticle.GetEraseFlag() = false;

                array_1d<float, 3> &vel = pparticle.GetVelocity();
                float &distance = pparticle.GetDistance();
                noalias(vel) = ZeroVector(3);
                distance = 0.0;

                for (unsigned int k = 0; k < (TDim + 1); k++)
                {
                    noalias(vel) += (N(j, k) * geom[k].FastGetSolutionStepValue(VELOCITY));
                    distance += N(j, k) * geom[k].FastGetSolutionStepValue(DISTANCE);
                }

                if (distance <= 0.0)
                    distance = -1.0;
                else
                    distance = 1.0;

                number_of_particles++;
            }
        }

        bool nonzero_mesh_velocity = false;
        //seeing if we have to use the mesh_velocity or not
        for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
             inode != mr_model_part.NodesEnd(); inode++)
        {
            const array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(MESH_VELOCITY);
            for (unsigned int i = 0; i != 3; i++)
            {
                if (fabs(velocity[i]) > 1.0e-9)
                    nonzero_mesh_velocity = true;
            }
            if (nonzero_mesh_velocity == true)
                break;
        }

        if (nonzero_mesh_velocity == true)
            muse_mesh_velocity_to_convect = true; // if there is mesh velocity, then we have to take it into account when moving the particles
        else
            muse_mesh_velocity_to_convect = false; //otherwise, we can avoid reading the values since we know it is zero everywhere (to save time!)

        //for now, I don't consider the mesh velocity, therefore I put here a throw error
        if (muse_mesh_velocity_to_convect == true)
            KRATOS_THROW_ERROR(std::logic_error, "muse_mesh_velocity_to_convect = true", "");

        m_nparticles = particle_id; //we save the last particle created as the total number of particles we have. For the moment this is true.
        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Number of Particles Created : " << m_nparticles << std::endl;
        mparticle_printing_tool_initialized = false;
    }

    ~MoveParticleUtilityFullBinsPFEM2()
    {
    }

    void MountBin(typename BinBasedFastPointLocator<TDim>::Pointer pSearchStructure)
    {
        KRATOS_TRY

        mpSearchStructure=pSearchStructure;

        KRATOS_CATCH("")
    }

    void CalculateVelOverElemSize()
    {
        KRATOS_TRY

        //ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();

        const double nodal_weight = 1.0 / (1.0 + double(TDim));

        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        vector<unsigned int> element_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif

        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);


        if (muse_mesh_velocity_to_convect == false)
        {
            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
                {
                    ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
                    Geometry<Node<3>> &geom = ielem->GetGeometry();
                    array_1d<double, 3> vector_mean_velocity = ZeroVector(3);

                    for (unsigned int i = 0; i != (TDim + 1); i++)
                        vector_mean_velocity += geom[i].FastGetSolutionStepValue(VELOCITY);
                    vector_mean_velocity *= nodal_weight;


                    const double mean_velocity = sqrt(pow(vector_mean_velocity[0], 2) + pow(vector_mean_velocity[1], 2) + pow(vector_mean_velocity[2], 2));
                    ielem->SetValue(VELOCITY_OVER_ELEM_SIZE, mean_velocity / (ielem->GetValue(MEAN_SIZE)));
                }
            }
        }
        else
        {
            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
                {
                    ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
                    Geometry<Node<3>> &geom = ielem->GetGeometry();

                    array_1d<double, 3> vector_mean_velocity = ZeroVector(3);

                    for (unsigned int i = 0; i != (TDim + 1); i++)
                        vector_mean_velocity += geom[i].FastGetSolutionStepValue(VELOCITY) - geom[i].FastGetSolutionStepValue(MESH_VELOCITY);
                    vector_mean_velocity *= nodal_weight;

                    const double mean_velocity = sqrt(pow(vector_mean_velocity[0], 2) + pow(vector_mean_velocity[1], 2) + pow(vector_mean_velocity[2], 2));
                    ielem->SetValue(VELOCITY_OVER_ELEM_SIZE, mean_velocity / (ielem->GetValue(MEAN_SIZE)));
                }
            }
        }

        KRATOS_CATCH("")
    }

    //name self explained
    void ResetBoundaryConditions(bool fully_reset_nodes)
    {
        KRATOS_TRY

        if (fully_reset_nodes)
        {
            ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin + ii;

                    if (inode->IsFixed(VELOCITY_X))
                    {
                        inode->FastGetSolutionStepValue(VELOCITY_X) = inode->GetSolutionStepValue(VELOCITY_X, 1);
                    }
                    if (inode->IsFixed(VELOCITY_Y))
                    {
                        inode->FastGetSolutionStepValue(VELOCITY_Y) = inode->GetSolutionStepValue(VELOCITY_Y, 1);
                    }
                    if (TDim == 3)
                    {
                        if (inode->IsFixed(VELOCITY_Z))
                        {
                            inode->FastGetSolutionStepValue(VELOCITY_Z) = inode->GetSolutionStepValue(VELOCITY_Z, 1);
                        }
                    }
                    if (inode->IsFixed(PRESSURE))
                        inode->FastGetSolutionStepValue(PRESSURE) = inode->GetSolutionStepValue(PRESSURE, 1);
                    //                          inode->GetSolutionStepValue(PRESSURE,1)=inode->FastGetSolutionStepValue(PRESSURE);
                }
            }
        }
        else //for fractional step only!
        {
            //For now only fully_reset_nodes is implemented. Therefore I put a throw error here.
            KRATOS_THROW_ERROR(std::logic_error, "Only fully_reset_nodes is implemented", "");
            ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin + ii;

                    const array_1d<double, 3> original_velocity = inode->FastGetSolutionStepValue(VELOCITY);

                    if (inode->IsFixed(VELOCITY_X) || inode->IsFixed(VELOCITY_Y) || inode->IsFixed(VELOCITY_Z))
                    {

                        const array_1d<double, 3> &normal = inode->FastGetSolutionStepValue(NORMAL);
                        const double normal_scalar_sq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
                        const array_1d<double, 3> normal_adimensionalized = normal / sqrt(normal_scalar_sq);
                        array_1d<double, 3> &velocity = inode->FastGetSolutionStepValue(VELOCITY);

                        array_1d<double, 3> normal_velocity;
                        for (unsigned int j = 0; j != 3; j++)
                            normal_velocity[j] = fabs(normal_adimensionalized[j]) * original_velocity[j];

                        if (inode->IsFixed(VELOCITY_X))
                        {
                            velocity[0] = original_velocity[0] - normal_velocity[0];
                        }
                        if (inode->IsFixed(VELOCITY_Y))
                        {
                            velocity[1] = original_velocity[1] - normal_velocity[1];
                        }
                        if (TDim == 3)
                            if (inode->IsFixed(VELOCITY_Z))
                            {
                                velocity[2] = original_velocity[2] - normal_velocity[2];
                            }
                    }

                    if (inode->IsFixed(PRESSURE))
                        inode->FastGetSolutionStepValue(PRESSURE) = inode->GetSolutionStepValue(PRESSURE, 1);
                }
            }
        }
        KRATOS_CATCH("")
    }

    //setting the normal component of the velocity to zero
    void ResetBoundaryConditionsSlip()
    {
        KRATOS_TRY

        {
            ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin + ii;

                    if (inode->Is(SLIP))
                    {

                        array_1d<double, 3> &velocity = inode->FastGetSolutionStepValue(VELOCITY);
                        const array_1d<double, 3> &normal = inode->FastGetSolutionStepValue(NORMAL);
                        const double normal_scalar_sq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
                        const array_1d<double, 3> normal_adimensionalized = normal / sqrt(normal_scalar_sq);
                        //calculating the normal component of the velocity
                        array_1d<double, 3> normal_velocity;
                        for (unsigned int j = 0; j != 3; j++)
                            normal_velocity[j] = normal_adimensionalized[j] * velocity[j];

                        const double dot_prod = normal_velocity[0] * velocity[0] + normal_velocity[1] * velocity[1] + normal_velocity[2] * velocity[2];
                        //if the dot product of velocity *  normal velocity is lower than zero, then they have opposite signs and we must invert the direction:
                        if (dot_prod < 0.0)
                            normal_velocity *= -1.0;

                        velocity -= normal_velocity; //substracting the normal component
                    }
                    else if (inode->IsFixed(VELOCITY_X) && inode->IsFixed(VELOCITY_Y))
                    {
                        inode->FastGetSolutionStepValue(VELOCITY) = inode->GetSolutionStepValue(VELOCITY, 1);
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }

    //name self explained
    void ResetBoundaryConditionsOnProjectedVelocity(bool fully_reset_nodes)
    {
        KRATOS_TRY

        if (fully_reset_nodes)
        {
            ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin + ii;

                    if (inode->IsFixed(VELOCITY_X))
                    {
                        inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) = inode->GetSolutionStepValue(VELOCITY_X, 1);
                    }
                    if (inode->IsFixed(VELOCITY_Y))
                    {
                        inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) = inode->GetSolutionStepValue(VELOCITY_Y, 1);
                    }
                    if (TDim == 3)
                    {
                        if (inode->IsFixed(VELOCITY_Z))
                        {
                            inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) = inode->GetSolutionStepValue(VELOCITY_Z, 1);
                        }
                    }
                    if (inode->IsFixed(PRESSURE))
                        inode->FastGetSolutionStepValue(PRESSURE) = inode->GetSolutionStepValue(PRESSURE, 1);
                    //                          inode->GetSolutionStepValue(PRESSURE,1)=inode->FastGetSolutionStepValue(PRESSURE);
                }
            }
        }
        else //for fractional step only!
        {
            //For now only fully_reset_nodes is implemented. Therefore I put a throw error here.
            KRATOS_THROW_ERROR(std::logic_error, "Only fully_reset_nodes is implemented", "");
            ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
            vector<unsigned int> node_partition;
            #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
            #else
            int number_of_threads = 1;
            #endif
            OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

            #pragma omp parallel for
            for (int kkk = 0; kkk < number_of_threads; kkk++)
            {
                for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
                {
                    ModelPart::NodesContainerType::iterator inode = inodebegin + ii;

                    const array_1d<double, 3> original_velocity = inode->FastGetSolutionStepValue(VELOCITY);

                    if (inode->IsFixed(VELOCITY_X) || inode->IsFixed(VELOCITY_Y) || inode->IsFixed(VELOCITY_Z))
                    {

                        const array_1d<double, 3> &normal = inode->FastGetSolutionStepValue(NORMAL);
                        const double normal_scalar_sq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
                        const array_1d<double, 3> normal_adimensionalized = normal / sqrt(normal_scalar_sq);
                        array_1d<double, 3> &velocity = inode->FastGetSolutionStepValue(VELOCITY);

                        array_1d<double, 3> normal_velocity;
                        for (unsigned int j = 0; j != 3; j++)
                            normal_velocity[j] = fabs(normal_adimensionalized[j]) * original_velocity[j];

                        if (inode->IsFixed(VELOCITY_X))
                        {
                            velocity[0] = original_velocity[0] - normal_velocity[0];
                        }
                        if (inode->IsFixed(VELOCITY_Y))
                        {
                            velocity[1] = original_velocity[1] - normal_velocity[1];
                        }
                        if (TDim == 3)
                            if (inode->IsFixed(VELOCITY_Z))
                            {
                                velocity[2] = original_velocity[2] - normal_velocity[2];
                            }
                    }

                    if (inode->IsFixed(PRESSURE))
                        inode->FastGetSolutionStepValue(PRESSURE) = inode->GetSolutionStepValue(PRESSURE, 1);
                }
            }
        }
        KRATOS_CATCH("")
    }

    void CalculateDeltaVelocity()
    {
        KRATOS_TRY
        ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                inode->FastGetSolutionStepValue(DELTA_VELOCITY) = inode->FastGetSolutionStepValue(VELOCITY) - inode->FastGetSolutionStepValue(PROJECTED_VELOCITY);
            }
        }

        KRATOS_CATCH("")
    }

    void CalculateDeltaVelocity_Nitsche()
    {
//         KRATOS_TRY
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
//         vector<unsigned int> node_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 inode->FastGetSolutionStepValue(DELTA_VELOCITY) = inode->FastGetSolutionStepValue(VELOCITY) - inode->FastGetSolutionStepValue(PROJECTED_VELOCITY);
//                 inode->FastGetSolutionStepValue(DELTA_VELOCITY_NITSCHE) = inode->FastGetSolutionStepValue(VELOCITY_NITSCHE) - inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE);
//             }
//         }

//         KRATOS_CATCH("")
    }

    void CopyVectorVarToPreviousTimeStep(const Variable<array_1d<double, 3>> &OriginVariable,
                                         ModelPart::NodesContainerType &rNodes)
    {
        KRATOS_TRY
        ModelPart::NodesContainerType::iterator inodebegin = rNodes.begin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                noalias(inode->GetSolutionStepValue(OriginVariable, 1)) = inode->FastGetSolutionStepValue(OriginVariable);
            }
        }
        KRATOS_CATCH("")
    }

    void CopyScalarVarToPreviousTimeStep(const Variable<double> &OriginVariable,
                                         ModelPart::NodesContainerType &rNodes)
    {
        KRATOS_TRY
        ModelPart::NodesContainerType::iterator inodebegin = rNodes.begin();
        vector<unsigned int> node_partition;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif
        OpenMPUtils::CreatePartition(number_of_threads, rNodes.size(), node_partition);

#pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                inode->GetSolutionStepValue(OriginVariable, 1) = inode->FastGetSolutionStepValue(OriginVariable);
            }
        }
        KRATOS_CATCH("")
    }

    //to move all the particles across the streamlines. heavy task!
    void MoveParticles(const bool discriminate_streamlines) //,const bool pressure_gradient_integrate)
    {
        KRATOS_TRY
        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        const array_1d<double, 3> gravity = CurrentProcessInfo[GRAVITY];
        Vector N(TDim + 1);
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        vector<unsigned int> element_partition;
        vector<unsigned int> node_partition;
        vector<unsigned int> particle_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
        OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();


        //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
            {
                ModelPart::ElementsContainerType::iterator old_element = ielembegin + ii;
                int &number_of_particles = old_element->GetValue(NUMBER_OF_FLUID_PARTICLES);
                number_of_particles = 0;
                //we reset the local vectors for a faster access;
            }
        }

        bool nonzero_mesh_velocity = false;
        muse_mesh_velocity_to_convect = false;
        //seeing if we have to use the mesh_velocity or not
        // for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin(); inode != mr_model_part.NodesEnd(); inode++)
        // {
        //     const array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(MESH_VELOCITY);
        //     for (unsigned int i = 0; i != 3; i++)
        //     {
        //         if (fabs(velocity[i]) > 1.0e-9)
        //             nonzero_mesh_velocity = true;
        //     }
        //     if (nonzero_mesh_velocity == true)
        //         break;
        // }

        //for now, I don't consider the mesh velocity, therefore I put here a throw error
        // if (nonzero_mesh_velocity == true)
        // {
        //     KRATOS_THROW_ERROR(std::logic_error, "there is mesh velocity", "");
        //     muse_mesh_velocity_to_convect = true; // if there is mesh velocity, then we have to take it into account when moving the particles
        // }
        // else
        // {
        //     muse_mesh_velocity_to_convect = false; //otherwise, we can avoid reading the values since we know it is zero everywhere (to save time!)
        // }
        std::cout << "Convecting Particles" << std::endl;
        //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)

        const bool local_use_mesh_velocity_to_convect = muse_mesh_velocity_to_convect;
        #pragma omp parallel for firstprivate(results,N)
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            const array_1d<double, 3> mesh_displacement = mcalculation_domain_added_displacement; //if it is a standard problem, displacements are zero and therefore nothing is added.
            for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
            {
                PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
                bool &erase_flag = pparticle.GetEraseFlag();
                if (erase_flag == false)
                {
                    array_1d<double, 3>& position = pparticle.Coordinates();
                    Element::Pointer pelement;
                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                    bool is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
                    if (is_found)
                    {
                        MoveParticle(position, pparticle, pelement, result_begin, max_results, mesh_displacement, discriminate_streamlines, N, delta_t, gravity);
                        //              if (number_of_particles_in_current_elem<mmaximum_number_of_particles && erase_flag==false)
                        if (erase_flag == false)
                        {
                            int &number_of_particles_in_current_elem = pelement->GetValue(NUMBER_OF_FLUID_PARTICLES);
                            #pragma omp atomic
                            number_of_particles_in_current_elem++;
                        }
                        else
                        {
                            pparticle.GetEraseFlag() = true; //so we just delete it!
                        }
                    }
                    else
                    {
                        erase_flag = true;
                        KRATOS_WATCH(pparticle);
                        KRATOS_THROW_ERROR(std::logic_error, "particle not found although it was erase_flag false","");
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void MoveParticles_Nitsche(const bool discriminate_streamlines) //,const bool pressure_gradient_integrate)
    {

//         KRATOS_TRY

//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         double delta_t = CurrentProcessInfo[DELTA_TIME];
//         const array_1d<double, 3> gravity = CurrentProcessInfo[GRAVITY];

//         array_1d<double, TDim + 1> N;
//         const unsigned int max_results = 10000;

//         vector<unsigned int> element_partition;
//         vector<unsigned int> node_partition;
//         vector<unsigned int> particle_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
//         ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

// //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::ElementsContainerType::iterator old_element = ielembegin + ii;
//                 int &number_of_particles = old_element->GetValue(NUMBER_OF_FLUID_PARTICLES);
//                 number_of_particles = 0;

//                 //we reset the local vectors for a faster access;
//             }
//         }

//         bool nonzero_mesh_velocity = false;
//         //seeing if we have to use the mesh_velocity or not
//         for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
//              inode != mr_model_part.NodesEnd(); inode++)
//         {
//             const array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(MESH_VELOCITY);
//             for (unsigned int i = 0; i != 3; i++)
//             {
//                 if (fabs(velocity[i]) > 1.0e-9)
//                     nonzero_mesh_velocity = true;
//             }
//             if (nonzero_mesh_velocity == true)
//                 break;
//         }

//         if (nonzero_mesh_velocity == true)
//             muse_mesh_velocity_to_convect = true; // if there is mesh velocity, then we have to take it into account when moving the particles
//         else
//             muse_mesh_velocity_to_convect = false; //otherwise, we can avoid reading the values since we know it is zero everywhere (to save time!)

//         std::cout << "convecting particles" << std::endl;
//         //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)

//         const bool local_use_mesh_velocity_to_convect = muse_mesh_velocity_to_convect;
//         //
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {

//             array_1d<double, 3> position;
//             const array_1d<double, 3> mesh_displacement = mcalculation_domain_added_displacement; //if it is a standard problem, displacements are zero and therefore nothing is added.

//             ResultContainerType results(max_results);
//             if ((results.size()) != max_results)
//                 results.resize(max_results);
//             array_1d<double, TDim + 1> N;

//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {
//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 N = ZeroVector(TDim + 1);
//                 position = pparticle.Coordinates();
//                 ResultIteratorType result_begin = results.begin();
//                 bool &erase_flag = pparticle.GetEraseFlag();
//                 if (erase_flag == false)
//                 {
//                     Element::Pointer pelement;
//                     bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                     if (isfound == true)
//                     {
//                         erase_flag = false;
//                     }
//                     else
//                     {
//                         erase_flag = true;
//                     }
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();
//                     if (erase_flag == false)
//                     {
//                         result_begin = results.begin();
//                         MoveParticle_Nitsche(pparticle, pelement, result_begin, max_results, mesh_displacement, discriminate_streamlines, N, delta_t, gravity);
//                         int &number_of_particles_in_current_elem = pelement->GetValue(NUMBER_OF_FLUID_PARTICLES);
//                         //                   if (number_of_particles_in_current_elem<mmaximum_number_of_particles && erase_flag==false)
//                         if (erase_flag == false)
//                         {
// #pragma omp atomic
//                             number_of_particles_in_current_elem++;
//                         }
//                         else
//                             pparticle.GetEraseFlag() = true; //so we just delete it!
//                     }
//                 }
//             }
//         }

//         KRATOS_CATCH("")
    }

    void MoveParticlesRK4(const bool discriminate_streamlines) //,const bool pressure_gradient_integrate)
    {

        // KRATOS_TRY

        // ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        // double delta_t = CurrentProcessInfo[DELTA_TIME];
        // const array_1d<double, 3> gravity = CurrentProcessInfo[GRAVITY];
        // const unsigned int max_results = 10000;

        // vector<unsigned int> element_partition;
        // vector<unsigned int> node_partition;
        // vector<unsigned int> particle_partition;
        // #ifdef _OPENMP
        // int number_of_threads = omp_get_max_threads();
        // #else
        // int number_of_threads = 1;
        // #endif
        // OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
        // OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
        // OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
        // ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        // ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

        // //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
        // #pragma omp parallel for
        // for (int kkk = 0; kkk < number_of_threads; kkk++)
        // {
        //     for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
        //     {
        //         ModelPart::ElementsContainerType::iterator old_element = ielembegin + ii;
        //         int &number_of_particles = old_element->GetValue(NUMBER_OF_FLUID_PARTICLES);
        //         number_of_particles = 0;

        //         //we reset the local vectors for a faster access;
        //     }
        // }

        // bool nonzero_mesh_velocity = false;
        // //seeing if we have to use the mesh_velocity or not
        // for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
        //      inode != mr_model_part.NodesEnd(); inode++)
        // {
        //     const array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(MESH_VELOCITY);
        //     for (unsigned int i = 0; i != 3; i++)
        //     {
        //         if (fabs(velocity[i]) > 1.0e-9)
        //             nonzero_mesh_velocity = true;
        //     }
        //     if (nonzero_mesh_velocity == true)
        //         break;
        // }

        // //for now, I don't consider the mesh velocity, therefore I put here a throw error
        // if (nonzero_mesh_velocity == true)
        // {
        //     KRATOS_THROW_ERROR(std::logic_error, "there is mesh velocity", "");
        //     muse_mesh_velocity_to_convect = true; // if there is mesh velocity, then we have to take it into account when moving the particles
        // }
        // else
        // {
        //     muse_mesh_velocity_to_convect = false; //otherwise, we can avoid reading the values since we know it is zero everywhere (to save time!)
        // }
        // std::cout << "Convecting Particles-RK4" << std::endl;
        // //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)

        // const bool local_use_mesh_velocity_to_convect = muse_mesh_velocity_to_convect;
        // #pragma omp parallel for
        // for (int kkk = 0; kkk < number_of_threads; kkk++)
        // {
        //     array_1d<double, 3> position;
        //     const array_1d<double, 3> mesh_displacement = mcalculation_domain_added_displacement; //if it is a standard problem, displacements are zero and therefore nothing is added.
        //     ResultContainerType results(max_results);
        //     array_1d<double, TDim + 1> N;
        //     for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
        //     {
        //         PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
        //         bool &erase_flag = pparticle.GetEraseFlag();
        //         if (erase_flag == false)
        //         {
        //             if ((results.size()) != max_results)
        //                 results.resize(max_results);

        //             N = ZeroVector(TDim + 1);
        //             position = pparticle.Coordinates();
        //             ResultIteratorType result_begin = results.begin();
        //             // KRATOS_WATCH(ielembegin);
        //             std::cout<<*ielembegin<<std::endl;
        //             KRATOS_THROW_ERROR(std::invalid_argument, "ROTATIONS ONLY IMPLEMENTED AROUND Z AXIS! (xy plane) ", "");
        //             Element::Pointer pelement(*(ielembegin + pparticle.GetElementId() - 1).base());
        //             Geometry<Node<3>> &geom = pelement->GetGeometry();
        //             bool isfound = CalculatePosition(geom, position[0], position[1], position[2], N);
        //             if (isfound == false)
        //             {
        //                 erase_flag = true;
        //                 KRATOS_WATCH(pparticle);
        //                 KRATOS_THROW_ERROR(std::logic_error, "particle not found although it was erase_flag false","");
        //             }
        //             else
        //             {
        //                 MoveParticleRK4(pparticle, pelement, result_begin, max_results, mesh_displacement, discriminate_streamlines, N, delta_t, gravity);
        //                 int &number_of_particles_in_current_elem = pelement->GetValue(NUMBER_OF_FLUID_PARTICLES);
        //                 //              if (number_of_particles_in_current_elem<mmaximum_number_of_particles && erase_flag==false)
        //                 if (erase_flag == false)
        //                 {
        //                     pparticle.GetElementId() = pelement->Id();
        //                     #pragma omp atomic
        //                     number_of_particles_in_current_elem++;
        //                 }
        //                 else
        //                 {
        //                     pparticle.GetElementId() = 0;
        //                     pparticle.GetEraseFlag() = true; //so we just delete it!
        //                 }
        //             }
        //         }
        //     }
        // }

        // KRATOS_CATCH("")
    }

    //this part has to be worked on
    void MoveParticlesRK4_Nitsche(const bool discriminate_streamlines) //,const bool pressure_gradient_integrate)
    {

    //     KRATOS_TRY

    //     ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
    //     double delta_t = CurrentProcessInfo[DELTA_TIME];
    //     const array_1d<double, 3> gravity = CurrentProcessInfo[GRAVITY];

    //     array_1d<double, TDim + 1> N;
    //     const unsigned int max_results = 10000;

    //     vector<unsigned int> element_partition;
    //     vector<unsigned int> node_partition;
    //     vector<unsigned int> particle_partition;
    //     #ifdef _OPENMP
    //     int number_of_threads = omp_get_max_threads();
    //     #else
    //     int number_of_threads = 1;
    //     #endif
    //     OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
    //     OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
    //     OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
    //     ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
    //     ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();

    //     //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
    //     #pragma omp parallel for
    //     for (int kkk = 0; kkk < number_of_threads; kkk++)
    //     {
    //         for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
    //         {
    //             ModelPart::ElementsContainerType::iterator old_element = ielembegin + ii;

    //             int &number_of_particles = old_element->GetValue(NUMBER_OF_FLUID_PARTICLES);
    //             number_of_particles = 0;
    //             //we reset the local vectors for a faster access;
    //         }
    //     }

    //     bool nonzero_mesh_velocity = false;
    //     //seeing if we have to use the mesh_velocity or not
    //     for (ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
    //          inode != mr_model_part.NodesEnd(); inode++)
    //     {
    //         const array_1d<double, 3> velocity = inode->FastGetSolutionStepValue(MESH_VELOCITY);
    //         for (unsigned int i = 0; i != 3; i++)
    //         {
    //             if (fabs(velocity[i]) > 1.0e-9)
    //                 nonzero_mesh_velocity = true;
    //         }
    //         if (nonzero_mesh_velocity == true)
    //             break;
    //     }

    //     //for now, I don't consider the mesh velocity, therefore I put here a throw error
    //     if (nonzero_mesh_velocity == true)
    //     {
    //         KRATOS_THROW_ERROR(std::logic_error, "there is mesh velocity", "");
    //         muse_mesh_velocity_to_convect = true; // if there is mesh velocity, then we have to take it into account when moving the particles
    //     }
    //     else
    //     {
    //         muse_mesh_velocity_to_convect = false; //otherwise, we can avoid reading the values since we know it is zero everywhere (to save time!)
    //     }
    //     std::cout << "Convecting Particles-RK4" << std::endl;
    //     //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)

    //     const bool local_use_mesh_velocity_to_convect = muse_mesh_velocity_to_convect;
    //     #pragma omp parallel for
    //     for (int kkk = 0; kkk < number_of_threads; kkk++)
    //     {

    //         array_1d<double, 3> position;
    //         const array_1d<double, 3> mesh_displacement = mcalculation_domain_added_displacement; //if it is a standard problem, displacements are zero and therefore nothing is added.
    //         ResultContainerType results(max_results);
    //         array_1d<double, TDim + 1> N;
    //         for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
    //         {
    //             PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
    //             bool &erase_flag = pparticle.GetEraseFlag();
    //             if (erase_flag == false)
    //             {
    //                 if ((results.size()) != max_results)
    //                     results.resize(max_results);
    //                 N = ZeroVector(TDim + 1);
    //                 position = pparticle.Coordinates();
    //                 ResultIteratorType result_begin = results.begin();
    //                 Element::Pointer pelement(*(ielembegin + pparticle.GetElementId() - 1).base());
    //                 Geometry<Node<3>> &geom = pelement->GetGeometry();
    //                 bool isfound = CalculatePosition(geom, position[0], position[1], position[2], N);
    //                 if (isfound == false)
    //                 {
    //                     erase_flag = true;
    //                     KRATOS_WATCH(pparticle);
    //                     KRATOS_THROW_ERROR(std::logic_error, "particle not found although it was erase_flag false");
                        
    //                 }
    //                 else
    //                 {
                       
    //                 Geometry<Node<3>> &geom = pelement->GetGeometry();
    //                 if (erase_flag == false)
    //                 {
    //                     result_begin = results.begin();
    //                     MoveParticleRK4_Nitsche(pparticle, pelement, result_begin, max_results, mesh_displacement, discriminate_streamlines, N, delta_t, gravity);
    //                     int &number_of_particles_in_current_elem = pelement->GetValue(NUMBER_OF_FLUID_PARTICLES);
    //                     //                   if (number_of_particles_in_current_elem<mmaximum_number_of_particles && erase_flag==false)
    //                     if (erase_flag == false)
    //                     {
    //                         #pragma omp atomic
    //                         number_of_particles_in_current_elem++;
    //                     }
    //                     else
    //                         pparticle.GetEraseFlag() = true; //so we just delete it!
    //                 }
    //             }
    //         }
    //     }

    //     KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerian() //explicit
    {
        //Parallel. Works!
        KRATOS_TRY
        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        const double threshold = 0.0 / (double(TDim) + 1.0);
        // auto t1 = std::chrono::high_resolution_clock::now();
        Vector N(TDim + 1);
        const int max_results = 10000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Projecting info to the mesh - first-order" << std::endl;

        vector<unsigned int> node_partition;
        vector<unsigned int> particle_partition;
        vector<unsigned int> element_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
        OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
        ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();

        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 0.0;
                inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = ZeroVector(3);
                inode->FastGetSolutionStepValue(YP) = 0.0;
            }
        }

        #pragma omp parallel for firstprivate(results,N)
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
            {
                PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
                bool &erase_flag = pparticle.GetEraseFlag();
                if (erase_flag == false)
                {
                    const array_1d<double, 3> &position = pparticle.Coordinates();
                    Element::Pointer pelement;
                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                    // bool is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
                    bool is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
                    if (is_found == true)
                    {
                        array_1d<double, 3 * (TDim + 1)> nodes_positions;
                        array_1d<double, 3 * (TDim + 1)> nodes_addedvel = ZeroVector(3 * (TDim + 1));
                        array_1d<double, (TDim + 1)> nodes_added_distance = ZeroVector((TDim + 1));
                        array_1d<double, (TDim + 1)> nodes_addedweights = ZeroVector((TDim + 1));
                        Geometry<Node<3>> &geom = pelement->GetGeometry();
                        const array_1d<float, 3>& velocity = pparticle.GetVelocity();
                        const float particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
                        for (int j = 0; j != (TDim + 1); j++) //going through the 3/4 nodes of the element
                        {
                            //double sq_dist = 0;
                            //these lines for a weighting function based on the distance (or square distance) from the node insteadof the shape functions
                            //for (int k=0 ; k!=(TDim); k++) sq_dist += ((position[k] - nodes_positions[j*3+k])*(position[k] - nodes_positions[j*3+k]));
                            //double weight = (1.0 - (sqrt(sq_dist)*weighting_inverse_divisor[j] ) );

                            double weight = N(j);
                            //weight=N(j)*N(j)*N(j);
                            if (weight < threshold)
                                weight = 1e-10;
                            if (weight < 0.0)
                            {
                                KRATOS_WATCH(weight)
                            } //;weight=0.0;KRATOS_WATCH(velocity);KRATOS_WATCH(N);KRATOS_WATCH(number_of_particles_in_elem);}//{KRATOS_WATCH(weight); KRATOS_WATCH(geom[j].Id()); KRATOS_WATCH(position);}
                            else
                            {
                                nodes_addedweights[j] += weight;
                                nodes_added_distance[j] += weight * particle_distance;
                                for (int k = 0; k != (TDim); k++) //x,y,(z)
                                {
                                    nodes_addedvel[j * 3 + k] += weight * double(velocity[k]);
                                }
                            }
                        }

                        for (int i = 0; i != (TDim + 1); ++i)
                        {
                            geom[i].SetLock();
                            geom[i].FastGetSolutionStepValue(PROJECTED_DISTANCE) += nodes_added_distance[i];
                            geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_X) += nodes_addedvel[3 * i + 0];
                            geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) += nodes_addedvel[3 * i + 1];
                            geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) += nodes_addedvel[3 * i + 2]; //we are updating info to the previous time step!!

                            geom[i].FastGetSolutionStepValue(YP) += nodes_addedweights[i];
                            geom[i].UnSetLock();
                        }
                    }
                }
            }
        }

        #pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                double sum_weights = inode->FastGetSolutionStepValue(YP);
                if (sum_weights > 0.00001)
                {
                    //inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT)=(inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT))/sum_weights; //resetting the temperature
                    double &dist = inode->FastGetSolutionStepValue(PROJECTED_DISTANCE);
                    dist /= sum_weights;                                                                                                       //resetting the density
                    inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = (inode->FastGetSolutionStepValue(PROJECTED_VELOCITY)) / sum_weights; //resetting the velocity
                }

                else //this should never happen because other ways to recover the information have been executed before, but leaving it just in case..
                {
                    inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 3.0; //resetting the temperature
                    //inode->FastGetSolutionStepValue(DISTANCE)=inode->GetSolutionStepValue(DISTANCE,1); //resetting the temperature
                    inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = inode->GetSolutionStepValue(VELOCITY, 1);
                }
                ///finally, if there was an inlet that had a fixed position for the distance function, that has to remain unchanged:
                if (inode->IsFixed(DISTANCE))
                    inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = inode->GetSolutionStepValue(DISTANCE, 1);
            }
        }

        KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerianDistance_Nitsche() //explicit
    {
//         KRATOS_TRY

//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         //const double delta_t =CurrentProcessInfo[DELTA_TIME];
//         const double threshold = 0.0 / (double(TDim) + 1.0);

//         KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Projecting info to mesh" << std::endl;

//         //          const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
//         //KRATOS_WATCH(offset)                                                                  //(flag managed only by MoveParticles

//         //we must project data from the particles (lagrangian)  into the eulerian mesh
//         //ValuesVectorType eulerian_nodes_old_temperature;
//         //int nnodes = mr_model_part.Nodes().size();
//         //array_1d<double,(n_nodes)> eulerian_nodes_sumweights;

//         //we save data from previous time step of the eulerian mesh in case we must reuse it later cos no particle was found around the nodes
//         //though we could've use a bigger buffer, to be changed later!
//         //after having saved data, we reset them to zero, this way it's easier to add the contribution of the surrounding particles.
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
//         vector<unsigned int> node_partition;
//         vector<unsigned int> particle_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 0.0;
//                 inode->FastGetSolutionStepValue(YP_DISTANCE) = 0.0;
//             }
//         }

//         //adding contribution, loop on elements, since each element has stored the particles found inside of it
//         vector<unsigned int> element_partition;
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

//         ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
//         const unsigned int max_results = 10000;
//         ResultContainerType results(max_results);
//         array_1d<double, 3> position;
//         array_1d<double, TDim + 1> N;

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {

//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 position = pparticle.Coordinates(); //initial coordinates
//                 const float particle_distance = pparticle.GetDistance();
//                 array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                 array_1d<double, 3> last_useful_vel;
//                 double sum_Ns_without_other_phase_nodes;
//                 ResultIteratorType result_begin = results.begin();
//                 double only_integral = 0.0;
//                 if ((results.size()) != max_results)
//                     results.resize(max_results);
//                 //we have to find on which particle we are (using bins).
//                 bool &erase_flag = pparticle.GetEraseFlag();
//                 Element::Pointer pelement(*ielembegin.base());
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);

//                 array_1d<double, 3 * (TDim + 1)> nodes_positions;
//                 array_1d<double, 3 * (TDim + 1)> nodes_addedvel = ZeroVector(3 * (TDim + 1));
//                 array_1d<double, 3 * (TDim + 1)> nodes_addedvel_nitsche = ZeroVector(3 * (TDim + 1));

//                 array_1d<double, (TDim + 1)> nodes_added_distance = ZeroVector((TDim + 1));
//                 array_1d<double, (TDim + 1)> nodes_addedweights = ZeroVector((TDim + 1));
//                 array_1d<double, (TDim + 1)> nodes_addedweights_distance = ZeroVector((TDim + 1));
//                 array_1d<double, (TDim + 1)> nodes_addedweights_nitsche = ZeroVector((TDim + 1));
//                 //array_1d<double,(TDim+1)> weighting_inverse_divisor;

//                 if (isfound == true && pparticle.GetEraseFlag() == false)
//                 {
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();
//                     bool interface_element = CheckIfInterfaceElement(pelement);

//                     for (int i = 0; i != (TDim + 1); ++i)
//                     {
//                         nodes_positions[i * 3 + 0] = geom[i].X();
//                         nodes_positions[i * 3 + 1] = geom[i].Y();
//                         nodes_positions[i * 3 + 2] = geom[i].Z();
//                         //weighting_inverse_divisor[i]=1.0/((geom[i].FastGetSolutionStepValue(MEAN_SIZE))*1.01);
//                     }

//                     array_1d<double, 3> &position = pparticle.Coordinates();

//                     const array_1d<float, 3> &velocity = pparticle.GetVelocity();

//                     const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air

//                     array_1d<double, TDim + 1> N;
//                     bool is_found = CalculatePosition(nodes_positions, position[0], position[1], position[2], N);
//                     if (is_found == false) //something went wrong. if it was close enough to the edge we simply send it inside the element.
//                     {
//                         KRATOS_WATCH(N);
//                         for (int j = 0; j != (TDim + 1); j++)
//                             if (N[j] < 0.0 && N[j] > -1e-5)
//                                 N[j] = 1e-10;
//                     }

//                     for (int j = 0; j != (TDim + 1); j++) //going through the 3/4 nodes of the element
//                     {
//                         //double sq_dist = 0;
//                         //these lines for a weighting function based on the distance (or square distance) from the node insteadof the shape functions
//                         //for (int k=0 ; k!=(TDim); k++) sq_dist += ((position[k] - nodes_positions[j*3+k])*(position[k] - nodes_positions[j*3+k]));
//                         //double weight = (1.0 - (sqrt(sq_dist)*weighting_inverse_divisor[j] ) );

//                         double weight = N(j);
//                         //weight=N(j)*N(j)*N(j);
//                         if (weight < threshold)
//                             weight = 1e-10;
//                         if (weight < 0.0)
//                         {
//                             KRATOS_WATCH(weight)
//                         } //;weight=0.0;KRATOS_WATCH(velocity);KRATOS_WATCH(N);KRATOS_WATCH(number_of_particles_in_elem);}//{KRATOS_WATCH(weight); KRATOS_WATCH(geom[j].Id()); KRATOS_WATCH(position);}
//                         else
//                         {
//                             //nodes_addedtemp[j] += weight * particle_temp;

//                             nodes_added_distance[j] += weight * particle_distance;
//                             nodes_addedweights_distance[j] += weight;

//                             //nodes_added_oxygen[j] += weight*particle_oxygen;

//                         } //
//                     }

// #pragma omp critical
//                     {
//                         for (int i = 0; i != (TDim + 1); ++i)
//                         {
//                             geom[i].SetLock();
//                             geom[i].FastGetSolutionStepValue(PROJECTED_DISTANCE) += nodes_added_distance[i];
//                             geom[i].FastGetSolutionStepValue(YP_DISTANCE) += nodes_addedweights_distance[i];
//                             geom[i].UnSetLock();
//                         }
//                     }
//                 }
//             }
//         }

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 double sum_weights_distance = inode->FastGetSolutionStepValue(YP_DISTANCE);
//                 if (sum_weights_distance > 0.00001)
//                 {
//                     //inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT)=(inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT))/sum_weights; //resetting the temperature
//                     double &dist = inode->FastGetSolutionStepValue(PROJECTED_DISTANCE);
//                     dist /= sum_weights_distance; //resetting the density
//                 }
//                 else //this should never happen because other ways to recover the information have been executed before, but leaving it just in case..
//                 {
//                     inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 3.0; //resetting the temperature
//                 }
//                 ///finally, if there was an inlet that had a fixed position for the distance function, that has to remain unchanged:
//                 if (inode->IsFixed(DISTANCE))
//                     inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = inode->GetSolutionStepValue(DISTANCE, 1);
//             }
//         }

//         KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerianVelocity_Nitsche() //explicit
    {
//         KRATOS_TRY

//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         //const double delta_t =CurrentProcessInfo[DELTA_TIME];
//         const double threshold = 0.0 / (double(TDim) + 1.0);

//         KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Projecting info to mesh" << std::endl;

//         //          const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
//         //KRATOS_WATCH(offset)                                                                  //(flag managed only by MoveParticles

//         //we must project data from the particles (lagrangian)  into the eulerian mesh
//         //ValuesVectorType eulerian_nodes_old_temperature;
//         //int nnodes = mr_model_part.Nodes().size();
//         //array_1d<double,(n_nodes)> eulerian_nodes_sumweights;

//         //we save data from previous time step of the eulerian mesh in case we must reuse it later cos no particle was found around the nodes
//         //though we could've use a bigger buffer, to be changed later!
//         //after having saved data, we reset them to zero, this way it's easier to add the contribution of the surrounding particles.
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
//         vector<unsigned int> node_partition;
//         vector<unsigned int> particle_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(YP) = 0.0;
//                 inode->FastGetSolutionStepValue(YP_NITSCHE) = 0.0;
//             }
//         }

//         //adding contribution, loop on elements, since each element has stored the particles found inside of it
//         vector<unsigned int> element_partition;
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

//         ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
//         const unsigned int max_results = 10000;
//         ResultContainerType results(max_results);
//         array_1d<double, 3> position;
//         array_1d<double, TDim + 1> N;

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {

//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 position = pparticle.Coordinates(); //initial coordinates
//                 const float particle_distance = pparticle.GetDistance();
//                 array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                 array_1d<double, 3> last_useful_vel;
//                 double sum_Ns_without_other_phase_nodes;
//                 ResultIteratorType result_begin = results.begin();
//                 double only_integral = 0.0;
//                 if ((results.size()) != max_results)
//                     results.resize(max_results);
//                 //we have to find on which particle we are (using bins).
//                 bool &erase_flag = pparticle.GetEraseFlag();
//                 Element::Pointer pelement(*ielembegin.base());
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);

//                 array_1d<double, 3 * (TDim + 1)> nodes_positions;
//                 array_1d<double, 3 * (TDim + 1)> nodes_addedvel = ZeroVector(3 * (TDim + 1));
//                 array_1d<double, 3 * (TDim + 1)> nodes_addedvel_nitsche = ZeroVector(3 * (TDim + 1));

//                 array_1d<double, (TDim + 1)> nodes_addedweights = ZeroVector((TDim + 1));
//                 array_1d<double, (TDim + 1)> nodes_addedweights_nitsche = ZeroVector((TDim + 1));
//                 //array_1d<double,(TDim+1)> weighting_inverse_divisor;

//                 if (isfound == true && pparticle.GetEraseFlag() == false)
//                 {
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();
//                     bool interface_element = CheckIfInterfaceElement(pelement);

//                     for (int i = 0; i != (TDim + 1); ++i)
//                     {
//                         nodes_positions[i * 3 + 0] = geom[i].X();
//                         nodes_positions[i * 3 + 1] = geom[i].Y();
//                         nodes_positions[i * 3 + 2] = geom[i].Z();
//                         //weighting_inverse_divisor[i]=1.0/((geom[i].FastGetSolutionStepValue(MEAN_SIZE))*1.01);
//                     }

//                     array_1d<double, 3> &position = pparticle.Coordinates();

//                     const array_1d<float, 3> &velocity = pparticle.GetVelocity();

//                     const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air

//                     array_1d<double, TDim + 1> N;
//                     bool is_found = CalculatePosition(nodes_positions, position[0], position[1], position[2], N);
//                     if (is_found == false) //something went wrong. if it was close enough to the edge we simply send it inside the element.
//                     {
//                         KRATOS_WATCH(N);
//                         for (int j = 0; j != (TDim + 1); j++)
//                             if (N[j] < 0.0 && N[j] > -1e-5)
//                                 N[j] = 1e-10;
//                     }

//                     for (int j = 0; j != (TDim + 1); j++) //going through the 3/4 nodes of the element
//                     {
//                         //double sq_dist = 0;
//                         //these lines for a weighting function based on the distance (or square distance) from the node insteadof the shape functions
//                         //for (int k=0 ; k!=(TDim); k++) sq_dist += ((position[k] - nodes_positions[j*3+k])*(position[k] - nodes_positions[j*3+k]));
//                         //double weight = (1.0 - (sqrt(sq_dist)*weighting_inverse_divisor[j] ) );

//                         double weight = N(j);
//                         //weight=N(j)*N(j)*N(j);
//                         if (weight < threshold)
//                             weight = 1e-10;
//                         if (weight < 0.0)
//                         {
//                             KRATOS_WATCH(weight)
//                         } //;weight=0.0;KRATOS_WATCH(velocity);KRATOS_WATCH(N);KRATOS_WATCH(number_of_particles_in_elem);}//{KRATOS_WATCH(weight); KRATOS_WATCH(geom[j].Id()); KRATOS_WATCH(position);}
//                         else
//                         {
//                             //nodes_addedtemp[j] += weight * particle_temp;

//                             //nodes_added_oxygen[j] += weight*particle_oxygen;

//                             if (interface_element) //the particle is inside an interface element
//                             {
//                                 if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                 {
//                                     nodes_addedweights_nitsche[j] += weight;
//                                     for (int k = 0; k != (TDim); k++) //x,y,(z)
//                                     {
//                                         nodes_addedvel_nitsche[j * 3 + k] += weight * double(velocity[k]);
//                                     }
//                                 }
//                                 else
//                                 {
//                                     nodes_addedweights[j] += weight;
//                                     for (int k = 0; k != (TDim); k++) //x,y,(z)
//                                     {
//                                         nodes_addedvel[j * 3 + k] += weight * double(velocity[k]);
//                                     }
//                                 }
//                             }
//                             else
//                             {
//                                 nodes_addedweights[j] += weight;
//                                 for (int k = 0; k != (TDim); k++) //x,y,(z)
//                                 {
//                                     nodes_addedvel[j * 3 + k] += weight * double(velocity[k]);
//                                 }
//                             }

//                         } //
//                     }

//                     #pragma omp critical
//                     {
//                         for (int i = 0; i != (TDim + 1); ++i)
//                         {
//                             geom[i].SetLock();
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_X) += nodes_addedvel[3 * i + 0];
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) += nodes_addedvel[3 * i + 1];
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) += nodes_addedvel[3 * i + 2]; //we are updating info to the previous time step!!
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X) += nodes_addedvel_nitsche[3 * i + 0];
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y) += nodes_addedvel_nitsche[3 * i + 1];
//                             geom[i].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z) += nodes_addedvel_nitsche[3 * i + 2]; //we are updating info to the previous time step!!

//                             geom[i].FastGetSolutionStepValue(YP) += nodes_addedweights[i];
//                             geom[i].FastGetSolutionStepValue(YP_NITSCHE) += nodes_addedweights_nitsche[i];
//                             geom[i].UnSetLock();
//                         }
//                     }
//                 }
//             }
//         }

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 double sum_weights = inode->FastGetSolutionStepValue(YP);
//                 double sum_weights_nitsche = inode->FastGetSolutionStepValue(YP_NITSCHE);
//                 if (sum_weights > 0.00001)
//                 {
//                     //inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT)=(inode->FastGetSolutionStepValue(TEMPERATURE_OLD_IT))/sum_weights; //resetting the temperature
//                     inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = (inode->FastGetSolutionStepValue(PROJECTED_VELOCITY)) / sum_weights; //resetting the velocity
//                 }
//                 else //this should never happen because other ways to recover the information have been executed before, but leaving it just in case..
//                 {
//                     KRATOS_WATCH(*inode);
//                     KRATOS_WATCH(sum_weights);
//                     inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 3.0; //resetting the temperature
//                                                                                //                       KRATOS_THROW_ERROR(std::logic_error, "sum_weights<0.00001", "");
//                                                                                //                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY)=inode->GetSolutionStepValue(VELOCITY,1);
//                 }
//                 if (sum_weights_nitsche > 0.00001)
//                 {
//                     inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE) = (inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE)) / sum_weights_nitsche; //resetting the velocity
//                 }
//             }
//         }

//         KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerianSecondOrder() //explicit
    {
        // KRATOS_TRY
        // ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        // const double time = CurrentProcessInfo[TIME];
        // const double delta_t = CurrentProcessInfo[DELTA_TIME];
        // const double threshold = 0.0 / (double(TDim) + 1.0);
        // std::cout << "Projecting info to mesh with Least Squares-Second Order" << std::endl;
        // vector<unsigned int> node_partition;
        // vector<unsigned int> particle_partition;
        // vector<unsigned int> element_partition;
        // #ifdef _OPENMP
        // int number_of_threads = omp_get_max_threads();
        // #else
        // int number_of_threads = 1;
        // #endif
        // // number_of_threads = 1; //Forced to be serial. There is too much bottleneck
        // OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
        // OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
        // OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
        // ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
        // ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        // #pragma omp parallel for
        // for (int kkk = 0; kkk < number_of_threads; kkk++)
        // {
        //     for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
        //     {
        //         ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
        //         inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 0.0;
        //         inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) = 0.0;
        //         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = ZeroVector(3);
        //         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY) = ZeroVector(3);
        //         inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE) = 0.0;
        //     }
        // }

        // //loop on particles and calculate the lumped mass matrix values for the nodes
        // #pragma omp parallel for
        // for (int kkk = 0; kkk < number_of_threads; kkk++)
        // {
        //     for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
        //     {
        //         PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
        //         bool &erase_flag = pparticle.GetEraseFlag();
        //         if (erase_flag == false)
        //         {
        //             array_1d<double, TDim + 1> N;
        //             const float particle_distance = pparticle.GetDistance();
        //             array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        //             Element::Pointer pelement(*(ielembegin + pparticle.GetElementId() - 1).base());
        //             Geometry<Node<3>> &geom = pelement->GetGeometry();
        //             array_1d<double, 3> &position = pparticle.Coordinates();
        //             bool isfound = CalculatePosition(geom, position[0], position[1], position[2], N);
        //             // bool isfound = FindNodeOnMesh( position, N, pelement, result_begin, max_results);
        //             //array_1d<double,(TDim+1)> weighting_inverse_divisor;
        //             if (isfound == true)
        //             {
        //                 for (int j = 0; j != (TDim + 1); j++) //going through the 3/4 nodes of the element
        //                 {
        //                     geom[j].SetLock();
        //                     double weight = N(j);
        //                     if (weight < threshold)
        //                         weight = 1e-10;
        //                     if (weight < 0.0)
        //                     {
        //                         KRATOS_WATCH(weight)
        //                     }
        //                     else
        //                     {
        //                         for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
        //                         {
        //                             geom[j].FastGetSolutionStepValue(LUMPED_MASS_VALUE) += N[j] * N[k];
        //                         }
        //                     }
        //                     geom[j].UnSetLock();
        //                 }
        //             }
        //         }
        //     }
        // }


        // //iteration for the calculation of the projected velocity and distance is starting
        // bool continue_velocity_iteration = true;
        // bool continue_distance_iteration = true;
        // double eps_velocity = 0.0;
        // double eps_distance = 0.0;
        // int whileloopcounter = 0;
        // double projecteddeltadistanceratio = 0.0;
        // double projecteddeltavelocityratio = 0.0;
        // while (continue_velocity_iteration == true and continue_distance_iteration == true)
        // {
        //     whileloopcounter++;
        //     #pragma omp parallel for
        //     for (int kkk = 0; kkk < number_of_threads; kkk++)
        //     {
        //         for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
        //         {
        //             ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
        //             inode->FastGetSolutionStepValue(R_NODE_VELOCITY) = ZeroVector(3);
        //             inode->FastGetSolutionStepValue(R_NODE_DISTANCE) = 0.0;
        //         }
        //     }
        //     //Right hand side is changing in every iteration, therefore they have to be set to zero at every iteration
        //     #pragma omp parallel for
        //     for (int kkk = 0; kkk < number_of_threads; kkk++)
        //     {
                
        //         for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
        //         {
        //             PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
        //             bool &erase_flag = pparticle.GetEraseFlag();
        //             if (erase_flag == false)
        //             {
        //                 array_1d<double, TDim + 1> N;
        //                 const float particle_distance = pparticle.GetDistance();
        //                 array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        //                 Element::Pointer pelement(*(ielembegin + pparticle.GetElementId() - 1).base());
        //                 Geometry<Node<3>> &geom = pelement->GetGeometry();
        //                 array_1d<double, 3> &position = pparticle.Coordinates();
        //                 bool isfound = CalculatePosition(geom, position[0], position[1], position[2], N);
        //                 //  bool isfound = FindNodeOnMesh( position, N, pelement, result_begin, max_results);
        //                 //array_1d<double,(TDim+1)> weighting_inverse_divisor;
        //                 if (isfound == true && pparticle.GetEraseFlag() == false)
        //                 {
        //                     Geometry<Node<3>> &geom = pelement->GetGeometry();
        //                     array_1d<double, 3> &position = pparticle.Coordinates();
        //                     const array_1d<float, 3> &velocity = pparticle.GetVelocity();
        //                     const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
        //                     for (int j = 0; j != (TDim + 1); j++)                     //going through the 3/4 nodes of the element
        //                     {
        //                         geom[j].SetLock();
        //                         double weight = N(j);
        //                         if (weight < threshold)
        //                             weight = 1e-10;
        //                         if (weight < 0.0)
        //                         {
        //                             KRATOS_WATCH(weight)
        //                         }
        //                         else
        //                         {
        //                             geom[j].FastGetSolutionStepValue(R_NODE_DISTANCE) += N[j] * particle_distance;
        //                             geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += N[j] * velocity[0];
        //                             geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += N[j] * velocity[1];
        //                             geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += N[j] * velocity[2];
        //                             for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
        //                             {
        //                                 geom[j].FastGetSolutionStepValue(R_NODE_DISTANCE) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_DISTANCE);
        //                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
        //                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
        //                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z);
        //                             }
        //                         }
        //                         geom[j].UnSetLock();
        //                     }
        //                 }
        //             }
        //         }
        //     }
        //     //right hand side, velocity and distance, for the current iteration has been calculated
        //     double sumprojecteddeltadistancenorm = 0.0;
        //     double sumprojecteddeltavelocitynorm = 0.0;
        //     double sumdistancenorm = 0.0;
        //     double sumprojectedvelocitynorm = 0.0;
        //     #pragma omp parallel for
        //     for (int kkk = 0; kkk < number_of_threads; kkk++)
        //     {
        //         for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
        //         {
        //             ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
        //             inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) = inode->FastGetSolutionStepValue(R_NODE_DISTANCE) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
        //             inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_X) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
        //             inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Y) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
        //             inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Z) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);

        //             inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE);
        //             inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X);
        //             inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y);
        //             inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z);

        //             double projecteddeltadistancenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE));
        //             double projecteddeltavelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z));

        //             double distancenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) * inode->FastGetSolutionStepValue(PROJECTED_DISTANCE));
        //             double projectedvelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z));
        //             #pragma omp critical
        //             {
        //                 sumprojecteddeltadistancenorm += projecteddeltadistancenorm;
        //                 sumprojecteddeltavelocitynorm += projecteddeltavelocitynorm;
        //                 sumdistancenorm += distancenorm;
        //                 sumprojectedvelocitynorm += projectedvelocitynorm;
        //             }
        //         }
        //     }
        //     projecteddeltadistanceratio = sumprojecteddeltadistancenorm / sumdistancenorm;
        //     projecteddeltavelocityratio = sumprojecteddeltavelocitynorm / sumprojectedvelocitynorm;
        //     // if (projecteddeltadistanceratio<0.00001 and projecteddeltavelocityratio<0.00001)
        //     //    if (projecteddeltadistanceratio<0.0001 and projecteddeltavelocityratio<0.0001)
        //     if (projecteddeltadistanceratio < 0.000001 and projecteddeltavelocityratio < 0.000001)
        //     {
        //         continue_distance_iteration = false;
        //         continue_velocity_iteration = false;
        //     }

        //     KRATOS_WATCH(whileloopcounter);
        //     KRATOS_WATCH(projecteddeltadistanceratio);
        //     KRATOS_WATCH(projecteddeltavelocityratio);
        // }
        // KRATOS_WATCH(whileloopcounter);
        // KRATOS_WATCH(projecteddeltadistanceratio);
        // KRATOS_WATCH(projecteddeltavelocityratio);

        // #pragma omp parallel for
        // for (int kkk = 0; kkk < number_of_threads; kkk++)
        // {
        //     for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
        //     {
        //         ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
        //         if (inode->IsFixed(DISTANCE))
        //             inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = inode->GetSolutionStepValue(DISTANCE, 1);
        //     }
        // }
        // KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerianSecondOrderDistance_Nitsche() //explicit
    {
//         KRATOS_TRY
//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         const double time = CurrentProcessInfo[TIME];
//         const double delta_t = CurrentProcessInfo[DELTA_TIME];
//         const double threshold = 0.0 / (double(TDim) + 1.0);
//         std::cout << "Projecting info to mesh with Least Squares-Second Order" << std::endl;
//         vector<unsigned int> node_partition;
//         vector<unsigned int> particle_partition;
//         vector<unsigned int> element_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = 0.0;
//                 inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) = 0.0;
//                 inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE) = 0.0;
//             }
//         }
// //loop on particles and calculate the lumped mass matrix values for the nodes
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {
//                 const unsigned int max_results = 10000;
//                 ResultContainerType results(max_results);
//                 array_1d<double, 3> position;
//                 array_1d<double, TDim + 1> N;
//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 position = pparticle.Coordinates(); //initial coordinates
//                 const float particle_distance = pparticle.GetDistance();
//                 array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                 array_1d<double, 3> last_useful_vel;
//                 if ((results.size()) != max_results)
//                     results.resize(max_results);
//                 ResultIteratorType result_begin = results.begin();
//                 //we have to find on which element we are (using bins).
//                 bool &erase_flag = pparticle.GetEraseFlag();
//                 Element::Pointer pelement;
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                 Geometry<Node<3>> &geom = pelement->GetGeometry();
//                 if (isfound == true && pparticle.GetEraseFlag() == false)
//                 {
//                     array_1d<double, 3> &position = pparticle.Coordinates();
//                     const array_1d<float, 3> &velocity = pparticle.GetVelocity();
//                     const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
//                     for (int j = 0; j != (TDim + 1); j++)                     //going through the 3/4 nodes of the element
//                     {
//                         geom[j].SetLock();
//                         double weight = N(j);
//                         if (weight < threshold)
//                             weight = 1e-10;
//                         if (weight < 0.0)
//                         {
//                             KRATOS_WATCH(weight)
//                         }
//                         else
//                         {
//                             for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
//                             {
//                                 geom[j].FastGetSolutionStepValue(LUMPED_MASS_VALUE) += N[j] * N[k];
//                             }
//                         }
//                         geom[j].UnSetLock();
//                     }
//                 }
//             }
//         }
//         //iteration for the calculation of the projected velocity and distance is starting
//         bool continue_distance_iteration = true;
//         double eps_distance = 0.0;
//         int whileloopcounter = 0;
//         double projecteddeltadistanceratio = 0.0;
//         while (continue_distance_iteration == true)
//         {
//             whileloopcounter++;
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//                 {
//                     ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                     inode->FastGetSolutionStepValue(R_NODE_DISTANCE) = 0.0;
//                 }
//             }
// //Right hand side is changing in every iteration, therefore they have to be set to zero at every iteration
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//                 {
//                     const unsigned int max_results = 10000;
//                     ResultContainerType results(max_results);
//                     array_1d<double, 3> position;
//                     array_1d<double, TDim + 1> N;
//                     PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                     position = pparticle.Coordinates(); //initial coordinates
//                     const float particle_distance = pparticle.GetDistance();
//                     array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                     array_1d<double, 3> last_useful_vel;
//                     double sum_Ns_without_other_phase_nodes;
//                     ResultIteratorType result_begin = results.begin();
//                     if ((results.size()) != max_results)
//                         results.resize(max_results);
//                     bool &erase_flag = pparticle.GetEraseFlag();
//                     Element::Pointer pelement;
//                     bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                     if (isfound == true && pparticle.GetEraseFlag() == false)
//                     {
//                         Geometry<Node<3>> &geom = pelement->GetGeometry();
//                         array_1d<double, 3> &position = pparticle.Coordinates();
//                         const array_1d<float, 3> &velocity = pparticle.GetVelocity();
//                         const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
//                         for (int j = 0; j != (TDim + 1); j++)                     //going through the 3/4 nodes of the element
//                         {
//                             geom[j].SetLock();
//                             double weight = N(j);
//                             if (weight < threshold)
//                                 weight = 1e-10;
//                             if (weight < 0.0)
//                             {
//                                 KRATOS_WATCH(weight)
//                             }
//                             else
//                             {
//                                 geom[j].FastGetSolutionStepValue(R_NODE_DISTANCE) += N[j] * particle_distance;
//                                 for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
//                                 {
//                                     geom[j].FastGetSolutionStepValue(R_NODE_DISTANCE) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_DISTANCE);
//                                 }
//                             }
//                             geom[j].UnSetLock();
//                         }
//                     }
//                 }
//             }
//             //right hand side, velocity and distance, for the current iteration has been calculated
//             double sumprojecteddeltadistancenorm = 0.0;
//             double sumprojecteddeltavelocitynorm = 0.0;
//             double sumdistancenorm = 0.0;
//             double sumprojectedvelocitynorm = 0.0;
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//                 {
//                     ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                     inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) = inode->FastGetSolutionStepValue(R_NODE_DISTANCE) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                     inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE);

//                     double projecteddeltadistancenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_DISTANCE));
//                     double distancenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) * inode->FastGetSolutionStepValue(PROJECTED_DISTANCE));

// #pragma omp critical
//                     {
//                         sumprojecteddeltadistancenorm += projecteddeltadistancenorm;
//                         sumdistancenorm += distancenorm;
//                     }
//                 }
//             }
//             projecteddeltadistanceratio = sumprojecteddeltadistancenorm / sumdistancenorm;
//             if (projecteddeltadistanceratio < 0.000000001)
//                 continue_distance_iteration = false;
//         }
//         KRATOS_WATCH(whileloopcounter);
//         KRATOS_WATCH(projecteddeltadistanceratio);

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 if (inode->IsFixed(DISTANCE))
//                     inode->FastGetSolutionStepValue(PROJECTED_DISTANCE) = inode->GetSolutionStepValue(DISTANCE, 1);
//             }
//         }
//         KRATOS_CATCH("")
    }

    void TransferLagrangianToEulerianSecondOrderVelocity_Nitsche() //explicit
    {
//         KRATOS_TRY
//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         const double time = CurrentProcessInfo[TIME];
//         const double delta_t = CurrentProcessInfo[DELTA_TIME];
//         const double threshold = 0.0 / (double(TDim) + 1.0);
//         std::cout << "Projecting info to mesh with Least Squares-Second Order" << std::endl;
//         vector<unsigned int> node_partition;
//         vector<unsigned int> particle_partition;
//         vector<unsigned int> element_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
//         ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//             {
//                 ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                 inode->FastGetSolutionStepValue(PROJECTED_VELOCITY) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE) = ZeroVector(3);
//                 inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE) = 0.0;
//                 inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE) = 0.0;
//             }
//         }
// //loop on particles and calculate the lumped mass matrix values for the nodes
// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {
//                 const unsigned int max_results = 10000;
//                 ResultContainerType results(max_results);
//                 array_1d<double, 3> position;
//                 array_1d<double, TDim + 1> N;
//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 position = pparticle.Coordinates(); //initial coordinates
//                 const float particle_distance = pparticle.GetDistance();
//                 array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                 array_1d<double, 3> last_useful_vel;
//                 if ((results.size()) != max_results)
//                     results.resize(max_results);
//                 ResultIteratorType result_begin = results.begin();
//                 //we have to find on which element we are (using bins).
//                 bool &erase_flag = pparticle.GetEraseFlag();
//                 Element::Pointer pelement;
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                 Geometry<Node<3>> &geom = pelement->GetGeometry();
//                 if (isfound == true && pparticle.GetEraseFlag() == false)
//                 {
//                     bool interface_element = CheckIfInterfaceElement(pelement);
//                     array_1d<double, 3> &position = pparticle.Coordinates();
//                     const array_1d<float, 3> &velocity = pparticle.GetVelocity();
//                     const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
//                     for (int j = 0; j != (TDim + 1); j++)                     //going through the 3/4 nodes of the element
//                     {
//                         geom[j].SetLock();
//                         double weight = N(j);
//                         if (weight < threshold)
//                             weight = 1e-10;
//                         if (weight < 0.0)
//                         {
//                             KRATOS_WATCH(weight)
//                         }
//                         else
//                         {
//                             if (interface_element)
//                             {
//                                 if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                 {
//                                     for (int k = 0; k != (TDim + 1); k++)
//                                     {
//                                         geom[j].FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE) += N[j] * N[k];
//                                     }
//                                 }
//                                 else
//                                 {
//                                     for (int k = 0; k != (TDim + 1); k++)
//                                     {
//                                         geom[j].FastGetSolutionStepValue(LUMPED_MASS_VALUE) += N[j] * N[k];
//                                     }
//                                 }
//                             }
//                             else
//                             {
//                                 for (int k = 0; k != (TDim + 1); k++)
//                                 {
//                                     geom[j].FastGetSolutionStepValue(LUMPED_MASS_VALUE) += N[j] * N[k];
//                                 }
//                             }
//                         }
//                         geom[j].UnSetLock();
//                     }
//                 }
//             }
//         }
//         //iteration for the calculation of the projected velocity and distance is starting
//         bool continue_velocity_iteration = true;
//         bool continue_distance_iteration = true;
//         double eps_velocity = 0.0;
//         double eps_distance = 0.0;
//         int whileloopcounter = 0;
//         double projecteddeltadistanceratio = 0.0;
//         double projecteddeltavelocityratio = 0.0;
//         while (continue_velocity_iteration == true and continue_distance_iteration == true)
//         {
//             whileloopcounter++;
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//                 {
//                     ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                     inode->FastGetSolutionStepValue(R_NODE_VELOCITY) = ZeroVector(3);
//                     inode->FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE) = ZeroVector(3);
//                 }
//             }
// //Right hand side is changing in every iteration, therefore they have to be set to zero at every iteration
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//                 {
//                     const unsigned int max_results = 10000;
//                     ResultContainerType results(max_results);
//                     array_1d<double, 3> position;
//                     array_1d<double, TDim + 1> N;
//                     PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                     position = pparticle.Coordinates(); //initial coordinates
//                     const float particle_distance = pparticle.GetDistance();
//                     array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
//                     array_1d<double, 3> last_useful_vel;
//                     double sum_Ns_without_other_phase_nodes;
//                     ResultIteratorType result_begin = results.begin();
//                     if ((results.size()) != max_results)
//                         results.resize(max_results);
//                     bool &erase_flag = pparticle.GetEraseFlag();
//                     Element::Pointer pelement;
//                     bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();
//                     if (isfound == true && pparticle.GetEraseFlag() == false)
//                     {
//                         bool interface_element = CheckIfInterfaceElement(pelement);
//                         array_1d<double, 3> &position = pparticle.Coordinates();
//                         const array_1d<float, 3> &velocity = pparticle.GetVelocity();
//                         const float &particle_distance = pparticle.GetDistance(); // -1 if water, +1 if air
//                         for (int j = 0; j != (TDim + 1); j++)                     //going through the 3/4 nodes of the element
//                         {
//                             geom[j].SetLock();
//                             double weight = N(j);
//                             if (weight < threshold)
//                                 weight = 1e-10;
//                             if (weight < 0.0)
//                             {
//                                 KRATOS_WATCH(weight)
//                             }
//                             else
//                             {
//                                 if (interface_element)
//                                 {
//                                     if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                     {
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_X) += N[j] * velocity[0];
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Y) += N[j] * velocity[1];
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Z) += N[j] * velocity[2];
//                                         for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
//                                         {
//                                             if (particle_distance * geom[k].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                             {
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z);
//                                             }
//                                             else
//                                             {
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z);
//                                             }
//                                         }
//                                     }
//                                     else
//                                     {
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += N[j] * velocity[0];
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += N[j] * velocity[1];
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += N[j] * velocity[2];
//                                         for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
//                                         {
//                                             if (particle_distance * geom[k].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                             {
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z);
//                                             }
//                                             else
//                                             {
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
//                                                 geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z);
//                                             }
//                                         }
//                                     }
//                                 }
//                                 else
//                                 {
//                                     geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += N[j] * velocity[0];
//                                     geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += N[j] * velocity[1];
//                                     geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += N[j] * velocity[2];
//                                     for (int k = 0; k != (TDim + 1); k++) //x,y,(z)
//                                     {
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_X) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Y) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
//                                         geom[j].FastGetSolutionStepValue(R_NODE_VELOCITY_Z) += -N[j] * N[k] * geom[k].FastGetSolutionStepValue(PROJECTED_VELOCITY_Z);
//                                     }
//                                 }
//                             }
//                             geom[j].UnSetLock();
//                         }
//                     }
//                 }
//             }
//             //right hand side, velocity and distance, for the current iteration has been calculated
//             double sumprojecteddeltadistancenorm = 0.0;
//             double sumprojecteddeltavelocitynorm = 0.0;
//             double sumdistancenorm = 0.0;
//             double sumprojectedvelocitynorm = 0.0;
// #pragma omp parallel for
//             for (int kkk = 0; kkk < number_of_threads; kkk++)
//             {
//                 for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
//                 {
//                     ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
//                     if (inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE) == 0.0) //not a dublicated node
//                     {
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_X) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Y) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Z) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);

//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z);

//                         double projecteddeltavelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z));
//                         double projectedvelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z));

// #pragma omp critical
//                         {
//                             sumprojecteddeltavelocitynorm += projecteddeltavelocitynorm;
//                             sumprojectedvelocitynorm += projectedvelocitynorm;
//                         }
//                     }
//                     else //a dublicated node
//                     {
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_X) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Y) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_Z) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_X) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_X) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Y) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Y) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE);
//                         inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Z) = inode->FastGetSolutionStepValue(R_NODE_VELOCITY_NITSCHE_Z) / inode->FastGetSolutionStepValue(LUMPED_MASS_VALUE_NITSCHE);

//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_X);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Y);
//                         inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z) += inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Z);

//                         double projecteddeltavelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_Z));
//                         double projectedvelocitynorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_X) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Y) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_Z));
//                         double projecteddeltavelocitynitschenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_X) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_X) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Y) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Y) + inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Z) * inode->FastGetSolutionStepValue(PROJECTED_DELTA_VELOCITY_NITSCHE_Z));
//                         double projectedvelocitynitschenorm = sqrt(inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_X) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Y) + inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z) * inode->FastGetSolutionStepValue(PROJECTED_VELOCITY_NITSCHE_Z));
// #pragma omp critical
//                         {

//                             sumprojecteddeltavelocitynorm += projecteddeltavelocitynorm;
//                             sumprojectedvelocitynorm += projectedvelocitynorm;
//                             sumprojecteddeltavelocitynorm += projecteddeltavelocitynitschenorm;
//                             sumprojectedvelocitynorm += projectedvelocitynitschenorm;
//                         }
//                     }
//                 }
//             }
//             projecteddeltavelocityratio = sumprojecteddeltavelocitynorm / sumprojectedvelocitynorm;
//             if (projecteddeltavelocityratio < 0.000000001)
//                 continue_velocity_iteration = false;
//         }
//         KRATOS_WATCH(whileloopcounter);
//         KRATOS_WATCH(projecteddeltavelocityratio);

//         KRATOS_CATCH("")
    }

    void AccelerateParticlesWithoutMovingUsingDeltaVelocity()
    {
        KRATOS_TRY
        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
        Vector N(TDim + 1);
        const int max_results = 1000;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        KRATOS_INFO("MoveParticleUtilityFullBinsPfem2") << "Accelerating (Correcting the Velocity of) the Particles " << std::endl;

        vector<unsigned int> element_partition;
        vector<unsigned int> particle_partition;
        #ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
        #else
        int number_of_threads = 1;
        #endif
        number_of_threads = 1;
        omp_set_num_threads(number_of_threads);
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
        OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);

        const double &delta_t = CurrentProcessInfo[DELTA_TIME];
        const array_1d<double, 3> &gravity = CurrentProcessInfo[GRAVITY];

        #pragma omp parallel for firstprivate(results,N)
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
            {
                PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
                bool &erase_flag = pparticle.GetEraseFlag();
                if (erase_flag == false)
                {
                    array_1d<double, 3> &position = pparticle.Coordinates();
                    Element::Pointer pelement;
                    typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
                    bool is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
                    if (is_found == true)
                    {
                        Geometry<Node<3>> &geom = pelement->GetGeometry();
                        AccelerateParticleUsingDeltaVelocity(pparticle, pelement, geom, delta_t, gravity, N); //'lite' version, we pass by reference the geometry, so much cheaper
                    }
                    else
                        pparticle.GetEraseFlag() = true;
                }
            }
        }

        KRATOS_CATCH("")
    }

    void AccelerateParticlesWithoutMovingUsingMeshDeltaVelocity()
    {
//         KRATOS_TRY
//         //std::cout << "updating particles" << std::endl;
//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

//         //          const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
//         //(flag managed only by MoveParticles
//         //KRATOS_WATCH(offset)
//         //          ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
//         const unsigned int max_results = 1000;

//         vector<unsigned int> element_partition;
//         vector<unsigned int> particle_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);

//         const double &delta_t = CurrentProcessInfo[DELTA_TIME];
//         const array_1d<double, 3> &gravity = CurrentProcessInfo[GRAVITY];

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             //initial coordinates
//             array_1d<double, 3> position;
//             array_1d<float, 3> particle_velocity;
//             float particle_distance;
//             //               array_1d<double,3> last_useful_vel;
//             //               array_1d<double,3> vel;

//             ResultContainerType results(max_results);
//             if ((results.size()) != max_results)
//                 results.resize(max_results);
//             array_1d<double, TDim + 1> N;

//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {
//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 //                 Element::Pointer pelement(*ielembegin.base());
//                 Element::Pointer pelement;

//                 particle_distance = pparticle.GetDistance();
//                 particle_velocity = pparticle.GetVelocity();
//                 position = pparticle.Coordinates(); //initial coordinates

//                 //                 double sum_Ns_without_other_phase_nodes;

//                 //                 double only_integral  = 0.0 ;
//                 ResultIteratorType result_begin = results.begin();
//                 //we have to find on which particle we are (using bins).
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                 if (isfound == true)
//                 {
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();

//                     AccelerateParticleUsingMeshDeltaVelocity(pparticle, pelement, geom, delta_t, gravity, N); //'lite' version, we pass by reference the geometry, so much cheaper
//                 }
//                 else
//                     pparticle.GetEraseFlag() = true;
//             }
//         }
//         KRATOS_CATCH("")
    }

    void AccelerateParticlesWithoutMovingUsingDeltaVelocity_Nitsche()
    {
//         KRATOS_TRY
//         //std::cout << "updating particles" << std::endl;
//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

//         //          const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
//         //(flag managed only by MoveParticles
//         //KRATOS_WATCH(offset)
//         //          ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();
//         const unsigned int max_results = 1000;

//         vector<unsigned int> element_partition;
//         vector<unsigned int> particle_partition;
// #ifdef _OPENMP
//         int number_of_threads = omp_get_max_threads();
// #else
//         int number_of_threads = 1;
// #endif
//         OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);
//         OpenMPUtils::CreatePartition(number_of_threads, mparticles_vector.size(), particle_partition);

//         const double &delta_t = CurrentProcessInfo[DELTA_TIME];
//         const array_1d<double, 3> &gravity = CurrentProcessInfo[GRAVITY];

// #pragma omp parallel for
//         for (int kkk = 0; kkk < number_of_threads; kkk++)
//         {
//             //initial coordinates
//             array_1d<double, 3> position;
//             array_1d<float, 3> particle_velocity;
//             float particle_distance;
//             //               array_1d<double,3> last_useful_vel;
//             //               array_1d<double,3> vel;

//             ResultContainerType results(max_results);
//             if ((results.size()) != max_results)
//                 results.resize(max_results);
//             array_1d<double, TDim + 1> N = ZeroVector(TDim + 1);

//             for (int ii = particle_partition[kkk]; ii < particle_partition[kkk + 1]; ii++)
//             {
//                 PFEM_Particle_Fluid &pparticle = mparticles_vector[ii];
//                 //                 Element::Pointer pelement(*ielembegin.base());
//                 Element::Pointer pelement;

//                 particle_distance = pparticle.GetDistance();
//                 particle_velocity = pparticle.GetVelocity();
//                 position = pparticle.Coordinates(); //initial coordinates

//                 //                 double sum_Ns_without_other_phase_nodes;

//                 //                 double only_integral  = 0.0 ;
//                 ResultIteratorType result_begin = results.begin();
//                 //we have to find on which particle we are (using bins).
//                 bool isfound = FindPositionUsingBins(position, N, pelement, result_begin, max_results);
//                 if (isfound == true)
//                 {
//                     Geometry<Node<3>> &geom = pelement->GetGeometry();

//                     AccelerateParticleUsingDeltaVelocity_Nitsche(pparticle, pelement, geom, delta_t, gravity, N); //'lite' version, we pass by reference the geometry, so much cheaper
//                 }
//                 else
//                     pparticle.GetEraseFlag() = true;
//             }
//         }
//         KRATOS_CATCH("")
    }
    

    //**************************************************************************************************************
    //**************************************************************************************************************

    template <class TDataType>
    void AddUniqueWeakPointer(GlobalPointersVector<TDataType> &v, const typename TDataType::WeakPointer candidate)
    {
        typename GlobalPointersVector<TDataType>::iterator i = v.begin();
        typename GlobalPointersVector<TDataType>::iterator endit = v.end();
        while (i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if (i == endit)
        {
            v.push_back(candidate);
        }
    }

    //**************************************************************************************************************
    //**************************************************************************************************************

    void PreReseed(int minimum_number_of_particles)
    {
        KRATOS_TRY

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        const int max_results = 10000;

        //tools for the paralelization
        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        omp_set_num_threads(number_of_threads);
        vector<unsigned int> elem_partition;
        int number_of_rows = mr_model_part.Elements().size();
        elem_partition.resize(number_of_threads + 1);
        int elem_partition_size = number_of_rows / number_of_threads;
        elem_partition[0] = 0;
        elem_partition[number_of_threads] = number_of_rows;
        //KRATOS_WATCH(elem_partition_size);
        for (unsigned int i = 1; i < number_of_threads; i++)
            elem_partition[i] = elem_partition[i - 1] + elem_partition_size;

        const bool local_use_mesh_velocity_to_convect = muse_mesh_velocity_to_convect;
        #pragma omp parallel firstprivate(elem_partition)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() + elem_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k + 1];
            //ModelPart::NodesContainerType local_list=aux[k];
            //PointerVectorSet<PFEM_Particle_Fluid, IndexedObject> & list=aux[k];
            //KRATOS_WATCH(k);
            unsigned int freeparticle = 0; //we start with the first position in the particles array
            //int local_id=1;
            for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
            {
                ResultContainerType results(max_results);
                results.resize(max_results);
                BoundedMatrix<double, (TDim + 1), 3> pos;
                BoundedMatrix<double, (TDim + 1), (TDim + 1)> N;
                //const int & elem_id = ielem->Id();
                int &number_of_particles_in_elem = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);
                if (number_of_particles_in_elem < (minimum_number_of_particles)) // && (ielem->GetGeometry())[0].Y()<0.10 )
                {
                    //KRATOS_WATCH("elem with little particles")
                    Geometry<Node<3>> &geom = ielem->GetGeometry();
                    ComputeGaussPointPositionsForPreReseed(geom, pos, N);
                    //double conductivity = ielem->GetProperties()[CONDUCTIVITY];
                    //KRATOS_WATCH(conductivity);
                    for (unsigned int j = 0; j < (pos.size1()); j++) //i am dropping the last one, the one in the middle of the element
                    {
                        bool keep_looking = true;
                        while (keep_looking)
                        {
                            if (mparticles_vector[freeparticle].GetEraseFlag() == true)
                            {
                                #pragma omp critical
                                {
                                    if (mparticles_vector[freeparticle].GetEraseFlag() == true)
                                    {
                                        mparticles_vector[freeparticle].GetEraseFlag() = false;
                                        keep_looking = false;
                                    }
                                }
                                if (keep_looking == false)
                                    break;
                                else
                                    freeparticle++;
                            }
                            else
                            {
                                freeparticle++;
                            }
                        }

                        PFEM_Particle_Fluid pparticle(pos(j, 0), pos(j, 1), pos(j, 2));

                        Vector aux2_N(TDim + 1);
                        bool is_found = CalculatePosition(geom, pos(j, 0), pos(j, 1), pos(j, 2), aux2_N);
                        if (is_found == false)
                        {
                            KRATOS_WATCH(aux2_N);
                        }

                        pparticle.GetEraseFlag() = false;

                        ResultIteratorType result_begin = results.begin();
                        Element::Pointer pelement(*ielem.base());
                        MoveParticle_inverse_way(pparticle, pelement, result_begin, max_results, local_use_mesh_velocity_to_convect);

                        //and we copy it to the array:
                        mparticles_vector[freeparticle] = pparticle;

                        pparticle.GetEraseFlag() = false;

                        number_of_particles_in_elem++;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    //**************************************************************************************************************
    //**************************************************************************************************************

    void PostReseed(int minimum_number_of_particles, double mass_correction_factor) //pooyan's way
    {
        KRATOS_TRY

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        if (mass_correction_factor > 0.5)
            mass_correction_factor = 0.5;
        if (mass_correction_factor < -0.5)
            mass_correction_factor = -0.5;
        //mass_correction_factor=0.0;

        //ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
        //const double delta_t = CurrentProcessInfo[DELTA_TIME];
        //array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
        //const int max_results = 1000;

        const double threshold = mass_correction_factor * 0.5;

        //TOOLS FOR THE PARALELIZATION
        //int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        //KRATOS_WATCH(number_of_threads);
        vector<unsigned int> elem_partition;
        int number_of_rows = mr_model_part.Elements().size();
        //KRATOS_WATCH(number_of_threads);
        //KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
        elem_partition.resize(number_of_threads + 1);
        int elem_partition_size = number_of_rows / number_of_threads;
        elem_partition[0] = 0;
        elem_partition[number_of_threads] = number_of_rows;
        //KRATOS_WATCH(elem_partition_size);
        for (unsigned int i = 1; i < number_of_threads; i++)
            elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
            //typedef Node < 3 > PointType;
            //std::vector<ModelPart::NodesContainerType> aux;// aux;
            //aux.resize(number_of_threads);

            //ModelPart::NodesContainerType::iterator it_begin_particle_model_part = mr_linea_model_part.NodesBegin();
            //ModelPart::NodesContainerType::iterator it_end_particle_model_part = mr_linea_model_part.NodesEnd();

        #pragma omp parallel firstprivate(elem_partition) // firstprivate(results)//we will add the nodes in different parts of aux and later assemple everything toghether, remaming particles ids to get consecutive ids
        {
            unsigned int reused_particles = 0;

            unsigned int freeparticle = 0; //we start by the first position;

            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() + elem_partition[k];
            ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k + 1];

            BoundedMatrix<double, (3 + 2 * TDim), 3> pos; //7 particles (2D) or 9 particles (3D)
            BoundedMatrix<double, (3 + 2 * TDim), (TDim + 1)> N;

            array_1d<double, 3> vel_complete, vel_without_air_nodes;
            double sum_Ns_without_air_nodes;
            double mesh_distance;

            array_1d<double, (3 + 2 * TDim)> distances;
            array_1d<int, (3 + 2 * TDim)> positions;
            array_1d<bool, (3 + 2 * TDim)> is_water_particle; //for both

            unsigned int number_of_reseeded_particles;
            //unsigned int number_of_water_reseeded_particles;

            //array_1d<double, 3 > nodes_distances;

            //int local_id=1;
            for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
            {
                //results.resize(max_results);

                int &number_of_particles_in_elem = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);

                Geometry<Node<3>> &geom = ielem->GetGeometry();
                if ((number_of_particles_in_elem < (minimum_number_of_particles))) // && (geom[0].Y()<0.10) ) || (number_of_water_particles_in_elem>2 && number_of_particles_in_elem<(minimum_number_of_particles) ) )
                {

                    //bool reseed_more=false;
                    number_of_reseeded_particles = 0;

                    //reseed_more=true;
                    number_of_reseeded_particles = 3 + 2 * TDim;
                    ComputeGaussPointPositionsForPostReseed(geom, pos, N);

                    distances = ZeroVector(3 + 2 * TDim);

                    bool has_water_node = false;
                    bool has_air_node = false;
                    double mean_element_distance = 0.0;

                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        mean_element_distance += (1.0 / double(TDim + 1)) * (geom[j].FastGetSolutionStepValue(DISTANCE));
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0)
                            has_water_node = true;
                        else
                            has_air_node = true;
                    }

                    //first we check the particle distance according to the nodal values
                    for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
                    {
                        positions[j] = j + 1; //just creating a vector from 1 to 7 or whathever our lenght is (7 for 2d, 9 for 3d)
                        for (unsigned int l = 0; l < (TDim + 1); l++)
                        {
                            distances[j] += N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
                        }
                    }

                    if ((has_air_node && has_water_node)) //for slit elements we use the distance function
                    {

                        for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
                        {
                            if (distances[j] > threshold)
                                is_water_particle[j] = false;
                            else
                                is_water_particle[j] = true;
                        }
                    }
                    else if (has_air_node)
                    {

                        double water_fraction = 0.5 - 0.5 * (mean_element_distance);
                        if (water_fraction > 0.9 && mass_correction_factor < 0.0) //to avoid seeding air particles when we are in a pure water element
                            mass_correction_factor = 0.0;
                        unsigned int number_of_water_reseeded_particles = double(number_of_reseeded_particles) * (1.01 + mass_correction_factor * 1.0) * water_fraction;

                        BubbleSort(distances, positions, number_of_reseeded_particles); //ok. now we have the particles ordered from the "watermost" to "airmost". therefore we will fill the water particles and later the air ones using that order

                        for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
                        {
                            int array_position = positions[j] - 1;
                            if (array_position > 3 && number_of_reseeded_particles == 4)
                            {
                                KRATOS_WATCH("error in reseeding")
                            }

                            if ((j + 1) <= number_of_water_reseeded_particles) //means it is a water particle
                                is_water_particle[array_position] = true;
                            else
                                is_water_particle[array_position] = false;
                        }
                    }
                    else //only water particles
                    {

                        for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
                            is_water_particle[j] = true;
                    }

                    bool fix_distance = false;
                    unsigned int node_with_fixed_distance = 0;
                    for (unsigned int j = 0; j < (TDim + 1); j++) //we go over the 3/4 nodes:
                    {
                        if ((geom[j].IsFixed(DISTANCE)))
                        {
                            fix_distance = true;
                            node_with_fixed_distance = j;
                        }
                    }

                    // so now if the 3 were fixed, we assign the sign of the first node to all the particles:
                    if (fix_distance)
                    {
                        bool is_water_for_all_particles = true;
                        if ((geom[node_with_fixed_distance].FastGetSolutionStepValue(DISTANCE)) > 0.0)
                            is_water_for_all_particles = false;

                        for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
                            is_water_particle[j] = is_water_for_all_particles;
                    }

                    for (unsigned int j = 0; j < number_of_reseeded_particles; j++)
                    {
                        //now we have to find an empty space ( a particle that was about to be deleted) in the particles model part. once found. there will be our renewed particle:
                        bool keep_looking = true;
                        while (keep_looking)
                        {
                            if (mparticles_vector[freeparticle].GetEraseFlag() == true)
                            {

                                #pragma omp critical
                                {

                                    if (mparticles_vector[freeparticle].GetEraseFlag() == true)
                                    {

                                        mparticles_vector[freeparticle].GetEraseFlag() = false;

                                        keep_looking = false;
                                    }
                                }
                                if (keep_looking == false)
                                {

                                    break;
                                }
                                else
                                {

                                    freeparticle++;
                                }
                            }
                            else
                            {
                                freeparticle++;
                            }
                        }

                        PFEM_Particle_Fluid pparticle(pos(j, 0), pos(j, 1), pos(j, 2));

                        array_1d<float, 3> &vel = pparticle.GetVelocity();
                        float &distance = pparticle.GetDistance();

                        Vector aux_N(TDim + 1);
                        bool is_found = CalculatePosition(geom, pos(j, 0), pos(j, 1), pos(j, 2), aux_N);
                        if (is_found == false)
                        {
                            KRATOS_WATCH(aux_N);
                            KRATOS_WATCH(j)
                            KRATOS_WATCH(ielem->Id())
                        }

                        noalias(vel_complete) = ZeroVector(3);
                        noalias(vel_without_air_nodes) = ZeroVector(3);
                        sum_Ns_without_air_nodes = 0.0;

                        noalias(vel) = ZeroVector(3);
                        distance = 0.0;
                        mesh_distance = 0.0;
                        //oxygen = 0.0;

                        for (unsigned int l = 0; l < (TDim + 1); l++)
                        {
                            noalias(vel_complete) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
                            mesh_distance += N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
                            if ((geom[l].FastGetSolutionStepValue(DISTANCE)) < 0.0)
                            {
                                sum_Ns_without_air_nodes += N(j, l);
                                noalias(vel_without_air_nodes) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
                            }
                        }

                        ///COMMENT TO GET A CONTINOUS DISTANCE FUNCTION FIELD
                        if (is_water_particle[j])
                        {
                            distance = -1.0;
                        }
                        else
                        {
                            //if (mesh_distance<2.0)
                            distance = 1.0;
                            //else
                            //  distance=3.0;
                        }

                        if (distance < 0.0 && sum_Ns_without_air_nodes > 0.01)
                            vel = vel_without_air_nodes / sum_Ns_without_air_nodes;
                        else
                            vel = vel_complete;

                        pparticle.GetEraseFlag() = false;

                        mparticles_vector[freeparticle] = pparticle;
                        //                          element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
                        number_of_particles_in_elem++;

                        if (keep_looking)
                        {
                            KRATOS_THROW_ERROR(std::logic_error, "FINISHED THE LIST AND COULDNT FIND A FREE CELL FOR THE NEW PARTICLE!", "");
                        }
                        else
                        {
                            reused_particles++;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void PostReseed_Nitsche(int minimum_number_of_particles, double mass_correction_factor) //pooyan's way
    {
//         KRATOS_TRY

//         ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET];

//         if (mass_correction_factor > 0.5)
//             mass_correction_factor = 0.5;
//         if (mass_correction_factor < -0.5)
//             mass_correction_factor = -0.5;
//         //mass_correction_factor=0.0;

//         //ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
//         //const double delta_t = CurrentProcessInfo[DELTA_TIME];
//         //array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
//         //const int max_results = 1000;

//         const double threshold = mass_correction_factor * 0.5;

//         //TOOLS FOR THE PARALELIZATION
//         //int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
//         unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
//         //KRATOS_WATCH(number_of_threads);
//         vector<unsigned int> elem_partition;
//         int number_of_rows = mr_model_part.Elements().size();
//         //KRATOS_WATCH(number_of_threads);
//         //KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
//         elem_partition.resize(number_of_threads + 1);
//         int elem_partition_size = number_of_rows / number_of_threads;
//         elem_partition[0] = 0;
//         elem_partition[number_of_threads] = number_of_rows;
//         //KRATOS_WATCH(elem_partition_size);
//         for (unsigned int i = 1; i < number_of_threads; i++)
//             elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
//             //typedef Node < 3 > PointType;
//             //std::vector<ModelPart::NodesContainerType> aux;// aux;
//             //aux.resize(number_of_threads);

//             //ModelPart::NodesContainerType::iterator it_begin_particle_model_part = mr_linea_model_part.NodesBegin();
//             //ModelPart::NodesContainerType::iterator it_end_particle_model_part = mr_linea_model_part.NodesEnd();

// #pragma omp parallel firstprivate(elem_partition) // firstprivate(results)//we will add the nodes in different parts of aux and later assemple everything toghether, remaming particles ids to get consecutive ids
//         {
//             unsigned int reused_particles = 0;

//             unsigned int freeparticle = 0; //we start by the first position;

//             int k = OpenMPUtils::ThisThread();
//             ModelPart::ElementsContainerType::iterator it_begin = mr_model_part.ElementsBegin() + elem_partition[k];
//             ModelPart::ElementsContainerType::iterator it_end = mr_model_part.ElementsBegin() + elem_partition[k + 1];

//             BoundedMatrix<double, (3 + 2 * TDim), 3> pos; //7 particles (2D) or 9 particles (3D)
//             BoundedMatrix<double, (3 + 2 * TDim), (TDim + 1)> N;

//             array_1d<double, 3> vel_complete, vel_without_air_nodes;
//             double sum_Ns_without_air_nodes;
//             double mesh_distance;

//             array_1d<double, (3 + 2 * TDim)> distances;
//             array_1d<int, (3 + 2 * TDim)> positions;
//             array_1d<bool, (3 + 2 * TDim)> is_water_particle; //for both

//             unsigned int number_of_reseeded_particles;
//             //unsigned int number_of_water_reseeded_particles;

//             //array_1d<double, 3 > nodes_distances;

//             //int local_id=1;
//             for (ModelPart::ElementsContainerType::iterator ielem = it_begin; ielem != it_end; ielem++)
//             {
//                 //results.resize(max_results);

//                 int &number_of_particles_in_elem = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);
//                 ParticlePointerVector &element_particle_pointers = (ielem->GetValue(FLUID_PARTICLE_POINTERS));

//                 Geometry<Node<3>> &geom = ielem->GetGeometry();
//                 if ((number_of_particles_in_elem < (minimum_number_of_particles))) // && (geom[0].Y()<0.10) ) || (number_of_water_particles_in_elem>2 && number_of_particles_in_elem<(minimum_number_of_particles) ) )
//                 {

//                     //bool reseed_more=false;
//                     number_of_reseeded_particles = 0;

//                     //reseed_more=true;
//                     number_of_reseeded_particles = 3 + 2 * TDim;
//                     ComputeGaussPointPositionsForPostReseed(geom, pos, N);

//                     distances = ZeroVector(3 + 2 * TDim);

//                     bool has_water_node = false;
//                     bool has_air_node = false;
//                     double mean_element_distance = 0.0;

//                     for (unsigned int j = 0; j < (TDim + 1); j++)
//                     {
//                         mean_element_distance += (1.0 / double(TDim + 1)) * (geom[j].FastGetSolutionStepValue(DISTANCE));
//                         if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0)
//                             has_water_node = true;
//                         else
//                             has_air_node = true;
//                     }

//                     //first we check the particle distance according to the nodal values
//                     for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
//                     {
//                         positions[j] = j + 1; //just creating a vector from 1 to 7 or whathever our lenght is (7 for 2d, 9 for 3d)
//                         for (unsigned int l = 0; l < (TDim + 1); l++)
//                         {
//                             distances[j] += N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
//                         }
//                     }

//                     if ((has_air_node && has_water_node)) //for slit elements we use the distance function
//                     {
//                         for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
//                         {
//                             if (distances[j] > threshold)
//                                 is_water_particle[j] = false;
//                             else
//                                 is_water_particle[j] = true;
//                         }
//                     }
//                     else if (has_air_node)
//                     {
//                         double water_fraction = 0.5 - 0.5 * (mean_element_distance);
//                         if (water_fraction > 0.9 && mass_correction_factor < 0.0) //to avoid seeding air particles when we are in a pure water element
//                             mass_correction_factor = 0.0;
//                         unsigned int number_of_water_reseeded_particles = double(number_of_reseeded_particles) * (1.01 + mass_correction_factor * 1.0) * water_fraction;

//                         BubbleSort(distances, positions, number_of_reseeded_particles); //ok. now we have the particles ordered from the "watermost" to "airmost". therefore we will fill the water particles and later the air ones using that order

//                         for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
//                         {
//                             int array_position = positions[j] - 1;
//                             if (array_position > 3 && number_of_reseeded_particles == 4)
//                             {
//                                 KRATOS_WATCH("error in reseeding")
//                             }

//                             if ((j + 1) <= number_of_water_reseeded_particles) //means it is a water particle
//                                 is_water_particle[array_position] = true;
//                             else
//                                 is_water_particle[array_position] = false;
//                         }
//                     }
//                     else //only water particles
//                     {
//                         for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
//                             is_water_particle[j] = true;
//                     }

//                     bool fix_distance = false;
//                     unsigned int node_with_fixed_distance = 0;
//                     for (unsigned int j = 0; j < (TDim + 1); j++) //we go over the 3/4 nodes:
//                     {
//                         if ((geom[j].IsFixed(DISTANCE)))
//                         {
//                             fix_distance = true;
//                             node_with_fixed_distance = j;
//                         }
//                     }
//                     // so now if the 3 were fixed, we assign the sign of the first node to all the particles:
//                     if (fix_distance)
//                     {
//                         bool is_water_for_all_particles = true;
//                         if ((geom[node_with_fixed_distance].FastGetSolutionStepValue(DISTANCE)) > 0.0)
//                             is_water_for_all_particles = false;

//                         for (unsigned int j = 0; j < number_of_reseeded_particles; j++) //first we order particles
//                             is_water_particle[j] = is_water_for_all_particles;
//                     }

//                     for (unsigned int j = 0; j < number_of_reseeded_particles; j++)
//                     {
//                         //now we have to find an empty space ( a particle that was about to be deleted) in the particles model part. once found. there will be our renewed particle:
//                         bool keep_looking = true;
//                         while (keep_looking)
//                         {
//                             if (mparticles_vector[freeparticle].GetEraseFlag() == true)
//                             {
// #pragma omp critical
//                                 {
//                                     if (mparticles_vector[freeparticle].GetEraseFlag() == true)
//                                     {
//                                         mparticles_vector[freeparticle].GetEraseFlag() = false;
//                                         keep_looking = false;
//                                     }
//                                 }
//                                 if (keep_looking == false)
//                                     break;

//                                 else
//                                     freeparticle++;
//                             }
//                             else
//                             {
//                                 freeparticle++;
//                             }
//                         }

//                         PFEM_Particle_Fluid pparticle(pos(j, 0), pos(j, 1), pos(j, 2));

//                         array_1d<float, 3> &vel = pparticle.GetVelocity();
//                         float &distance = pparticle.GetDistance();

//                         array_1d<double, TDim + 1> aux_N;
//                         bool is_found = CalculatePosition(geom, pos(j, 0), pos(j, 1), pos(j, 2), aux_N);
//                         if (is_found == false)
//                         {
//                             KRATOS_WATCH(aux_N);
//                             KRATOS_WATCH(j)
//                             KRATOS_WATCH(ielem->Id())
//                         }

//                         noalias(vel_complete) = ZeroVector(3);
//                         noalias(vel_without_air_nodes) = ZeroVector(3);
//                         sum_Ns_without_air_nodes = 0.0;

//                         noalias(vel) = ZeroVector(3);
//                         distance = 0.0;
//                         mesh_distance = 0.0;
//                         //oxygen = 0.0;

//                         ///COMMENT TO GET A CONTINOUS DISTANCE FUNCTION FIELD
//                         if (is_water_particle[j])
//                         {
//                             distance = -1.0;
//                         }
//                         else
//                         {
//                             //if (mesh_distance<2.0)
//                             distance = 1.0;
//                             //else
//                             //  distance=3.0;
//                         }

//                         bool interface_element = CheckIfInterfaceElement(ielem);
//                         if (interface_element) //the particle is inside an interface element
//                         {
//                             for (unsigned int l = 0; l < (TDim + 1); l++)
//                             {
//                                 if (distance * geom[l].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
//                                     noalias(vel) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY_NITSCHE);
//                                 else
//                                     noalias(vel) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
//                             }
//                         }
//                         else
//                         {
//                             for (unsigned int l = 0; l < (TDim + 1); l++)
//                             {
//                                 noalias(vel_complete) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
//                                 mesh_distance += N(j, l) * geom[l].FastGetSolutionStepValue(DISTANCE);
//                                 if ((geom[l].FastGetSolutionStepValue(DISTANCE)) < 0.0)
//                                 {
//                                     sum_Ns_without_air_nodes += N(j, l);
//                                     noalias(vel_without_air_nodes) += N(j, l) * geom[l].FastGetSolutionStepValue(VELOCITY);
//                                 }
//                             }

//                             if (distance < 0.0 && sum_Ns_without_air_nodes > 0.01)
//                                 vel = vel_without_air_nodes / sum_Ns_without_air_nodes;
//                             else
//                                 vel = vel_complete;
//                         }

//                         pparticle.GetEraseFlag() = false;

//                         mparticles_vector[freeparticle] = pparticle;
//                         //                          element_particle_pointers(offset+number_of_particles_in_elem) = &mparticles_vector[freeparticle];
//                         number_of_particles_in_elem++;

//                         if (keep_looking)
//                         {
//                             KRATOS_THROW_ERROR(std::logic_error, "FINISHED THE LIST AND COULDNT FIND A FREE CELL FOR THE NEW PARTICLE!", "");
//                         }
//                         else
//                         {
//                             reused_particles++;
//                         }
//                     }
//                 }
//             }
//         }

//         KRATOS_CATCH("")
    }

    void ExecuteParticlesPritingTool(ModelPart &lagrangian_model_part, int input_filter_factor)
    {
        KRATOS_TRY
        //mfilter_factor; //we will only print one out of every "filter_factor" particles of the total particle list

        //          if(mparticle_printing_tool_initialized==false)
        //          {
        mfilter_factor = input_filter_factor;

        //              if(lagrangian_model_part.NodesBegin()-lagrangian_model_part.NodesEnd()>0)
        //                  KRATOS_THROW_ERROR(std::logic_error, "AN EMPTY MODEL PART IS REQUIRED FOR THE PRINTING OF PARTICLES", "");

        lagrangian_model_part.Nodes().clear();

        lagrangian_model_part.AddNodalSolutionStepVariable(VELOCITY);
        lagrangian_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
        lagrangian_model_part.AddNodalSolutionStepVariable(DISTANCE);
        lagrangian_model_part.AddNodalSolutionStepVariable(ELEMENT_ID);
        int counter = 0;
        for (unsigned int i = 0; i != (mmaximum_number_of_particles * mnelems); i++)
        {
            PFEM_Particle_Fluid &pparticle = mparticles_vector[i];
            if (pparticle.GetEraseFlag() == false)
            {
                Node<3>::Pointer pnode = lagrangian_model_part.CreateNewNode(counter + mlast_node_id + 1, 0.0, 0.0, 0.0); //recordar que es el nueevo model part!!
                //pnode->SetBufferSize(mr_model_part.NodesBegin()->GetBufferSize());
                pnode->SetBufferSize(1);
                counter++;
            }
        }
        mparticle_printing_tool_initialized = true;
        //          }

        //resetting data of the unused particles
        const double inactive_particle_position = -10.0;
        array_1d<double, 3> inactive_particle_position_vector;
        inactive_particle_position_vector(0) = inactive_particle_position;
        inactive_particle_position_vector(1) = inactive_particle_position;
        inactive_particle_position_vector(2) = inactive_particle_position;
        ModelPart::NodesContainerType::iterator inodebegin = lagrangian_model_part.NodesBegin();
        for (unsigned int ii = 0; ii < lagrangian_model_part.Nodes().size(); ii++)
        {
            ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
            inode->FastGetSolutionStepValue(DISTANCE) = 0.0;
            inode->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
            inode->FastGetSolutionStepValue(DISPLACEMENT) = inactive_particle_position_vector;
            inode->FastGetSolutionStepValue(ELEMENT_ID) = 0.0;
        }

        counter = 0;
        //ModelPart::NodesContainerType::iterator it_begin = lagrangian_model_part.NodesBegin();
        for (int i = 0; i != mmaximum_number_of_particles * mnelems; i++)
        {
            PFEM_Particle_Fluid &pparticle = mparticles_vector[i];
            if (pparticle.GetEraseFlag() == false)
            //                 if(pparticle.GetEraseFlag()==false && i%mfilter_factor==0 && pparticle.X()>-1.0 && pparticle.X()<1.0 && pparticle.Y()>-1.0 && pparticle.Y()>1.0)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + counter; //copying info from the particle to the (printing) node.
                inode->FastGetSolutionStepValue(DISTANCE) = pparticle.GetDistance();
                inode->FastGetSolutionStepValue(VELOCITY) = pparticle.GetVelocity();
                inode->FastGetSolutionStepValue(DISPLACEMENT) = pparticle.Coordinates();
                counter++;
            }
        }

        KRATOS_CATCH("")
    }

    void ExecuteParticlesPritingToolForDroppletsOnly(ModelPart &lagrangian_model_part, int input_filter_factor)
    {
        KRATOS_TRY
        //mfilter_factor; //we will only print one out of every "filter_factor" particles of the total particle list
        const int first_particle_id = 1000000;
        if (mparticle_printing_tool_initialized == false)
        {
            mfilter_factor = input_filter_factor;

            if (lagrangian_model_part.NodesBegin() - lagrangian_model_part.NodesEnd() > 0)
                KRATOS_THROW_ERROR(std::logic_error, "AN EMPTY MODEL PART IS REQUIRED FOR THE PRINTING OF PARTICLES", "");

            lagrangian_model_part.AddNodalSolutionStepVariable(VELOCITY);
            lagrangian_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
            lagrangian_model_part.AddNodalSolutionStepVariable(DISTANCE);

            for (unsigned int i = 0; i != ((mmaximum_number_of_particles * mnelems) / mfilter_factor) + mfilter_factor; i++)
            {
                Node<3>::Pointer pnode = lagrangian_model_part.CreateNewNode(i + first_particle_id + 1, 0.0, 0.0, 0.0); //recordar que es el nueevo model part!!
                //pnode->SetBufferSize(mr_model_part.NodesBegin()->GetBufferSize());
                pnode->SetBufferSize(1);
            }
            mparticle_printing_tool_initialized = true;
        }

        //resetting data of the unused particles
        const double inactive_particle_position = -10.0;
        array_1d<double, 3> inactive_particle_position_vector;
        inactive_particle_position_vector(0) = inactive_particle_position;
        inactive_particle_position_vector(1) = inactive_particle_position;
        inactive_particle_position_vector(2) = inactive_particle_position;
        ModelPart::NodesContainerType::iterator inodebegin = lagrangian_model_part.NodesBegin();
        for (unsigned int ii = 0; ii < lagrangian_model_part.Nodes().size(); ii++)
        {
            ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
            inode->FastGetSolutionStepValue(DISTANCE) = 0.0;
            inode->FastGetSolutionStepValue(VELOCITY) = ZeroVector(3);
            inode->FastGetSolutionStepValue(DISPLACEMENT) = inactive_particle_position_vector;
        }
        const int max_number_of_printed_particles = lagrangian_model_part.Nodes().size();

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
                                                                               //(flag managed only by MoveParticles
        //KRATOS_WATCH(offset)
        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();

        int counter = 0;
        for (unsigned int ii = 0; ii < mr_model_part.Elements().size(); ii++)
        {
            ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
            Element::Pointer pelement(*ielem.base());
            Geometry<Node<3>> &geom = ielem->GetGeometry();
            //double mean_elem_dist=0.0;
            bool pure_air_elem = true;
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                if (geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0)
                    pure_air_elem = false;
                //mean_elem_dist += geom[j].FastGetSolutionStepValue(DISTANCE);
            }
            //if (mean_elem_dist>0.0) //only air elements
            if (pure_air_elem == true)
            {
                ParticlePointerVector &element_particle_pointers = (ielem->GetValue(FLUID_PARTICLE_POINTERS));
                int &number_of_particles_in_elem = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);
                //std::cout << "elem " << ii << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;

                for (int iii = 0; iii < number_of_particles_in_elem; iii++)
                {
                    //KRATOS_WATCH(iii)
                    if (iii > mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
                        break;

                    PFEM_Particle_Fluid &pparticle = element_particle_pointers[offset + iii];

                    bool erase_flag = pparticle.GetEraseFlag();
                    if (erase_flag == false && pparticle.GetDistance() < 0.0)
                    {
                        ModelPart::NodesContainerType::iterator inode = inodebegin + counter; //copying info from the particle to the (printing) node.
                        inode->FastGetSolutionStepValue(DISTANCE) = pparticle.GetDistance();
                        inode->FastGetSolutionStepValue(VELOCITY) = pparticle.GetVelocity();
                        inode->FastGetSolutionStepValue(DISPLACEMENT) = pparticle.Coordinates();
                        counter++;
                    }
                }
            }

            if (counter > (max_number_of_printed_particles - 30)) //we are approaching the end of the model part. so we stop before it's too late
                break;
        }

        KRATOS_CATCH("")
    }

    void AssignNodalVelocityUsingInletConditions(const double inlet_vel)
    {
        KRATOS_TRY

        //first we are going to delete all the velocities!
        ModelPart::ConditionsContainerType::iterator iconditionbegin = mr_model_part.ConditionsBegin();
        vector<unsigned int> condition_partition;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif

        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Conditions().size(), condition_partition);

#pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = condition_partition[kkk]; ii < condition_partition[kkk + 1]; ii++)
            {
                ModelPart::ConditionsContainerType::iterator icondition = iconditionbegin + ii;
                if (icondition->GetValue(IS_INLET) > 0.5)
                {
                    Geometry<Node<3>> &geom = icondition->GetGeometry();
                    array_1d<double, 3> normal = ZeroVector(3);
                    this->CalculateNormal(geom, normal);
                    const double normal_lenght = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
                    const array_1d<double, 3> velocity = -inlet_vel / normal_lenght * normal;
                    for (unsigned int l = 0; l < (TDim); l++)
                    {
                        geom[l].SetLock();
                        geom[l].FastGetSolutionStepValue(VELOCITY) = velocity;
                        geom[l].UnSetLock();
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void RotateParticlesAndDomainVelocities(array_1d<double, 3> rotations)
    {
        KRATOS_TRY

        if (fabs(rotations[0]) > 0.000000001 || fabs(rotations[1]) > 0.000000001)
            KRATOS_THROW_ERROR(std::invalid_argument, "ROTATIONS ONLY IMPLEMENTED AROUND Z AXIS! (xy plane) ", "");

        const double cosinus_theta = cos(rotations[2]);
        const double sinus_theta = sin(rotations[2]);

        //std::cout << "updating particles" << std::endl;
        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        const int offset = CurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET]; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
                                                                               //(flag managed only by MoveParticles
        //KRATOS_WATCH(offset)
        ModelPart::ElementsContainerType::iterator ielembegin = mr_model_part.ElementsBegin();

        vector<unsigned int> element_partition;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Elements().size(), element_partition);

#pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = element_partition[kkk]; ii < element_partition[kkk + 1]; ii++)
            {
                //const int & elem_id = ielem->Id();
                ModelPart::ElementsContainerType::iterator ielem = ielembegin + ii;
                Element::Pointer pelement(*ielem.base());

                ParticlePointerVector &element_particle_pointers = (ielem->GetValue(FLUID_PARTICLE_POINTERS));
                int &number_of_particles_in_elem = ielem->GetValue(NUMBER_OF_FLUID_PARTICLES);
                //std::cout << "elem " << ii << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;

                for (int iii = 0; iii < number_of_particles_in_elem; iii++)
                {
                    //KRATOS_WATCH(iii)
                    if (iii > mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
                        break;

                    PFEM_Particle_Fluid &pparticle = element_particle_pointers[offset + iii];

                    bool erase_flag = pparticle.GetEraseFlag();
                    if (erase_flag == false)
                    {
                        array_1d<float, 3> &vel = pparticle.GetVelocity();
                        const float vel_x = vel[0];
                        const float vel_y = vel[1];
                        vel[0] = cosinus_theta * vel_x + sinus_theta * vel_y;
                        vel[1] = cosinus_theta * vel_y - sinus_theta * vel_x;
                    }
                }
            }
        }

        ModelPart::NodesContainerType::iterator inodebegin = mr_model_part.NodesBegin();
        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, mr_model_part.Nodes().size(), node_partition);

#pragma omp parallel for
        for (int kkk = 0; kkk < number_of_threads; kkk++)
        {
            for (unsigned int ii = node_partition[kkk]; ii < node_partition[kkk + 1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + ii;
                if (inode->IsFixed(VELOCITY_X) == false)
                {
                    array_1d<double, 3> &vel = inode->FastGetSolutionStepValue(VELOCITY);
                    const double vel_x = vel[0];
                    const double vel_y = vel[1];
                    vel[0] = cosinus_theta * vel_x + sinus_theta * vel_y;
                    vel[1] = cosinus_theta * vel_y - sinus_theta * vel_x;
                }
            }
        }
        KRATOS_CATCH("")
    }

protected:
private:
    void Check()
    {
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing DISTANCE variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing VELOCITY variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing PRESSURE variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(PROJECTED_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing PROJECTED_VELOCITY variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(DELTA_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing DELTA_VELOCITY variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing MESH_VELOCITY variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(YP) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing YP variable on solution step data", "");
        if (mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false)
            KRATOS_THROW_ERROR(std::invalid_argument, "missing NORMAL variable on solution step data", "");
    }

    ///this function moves a particle according to the "velocity" given
    ///by "rVariable". The movement is performed in nsubsteps, during a total time
    ///of Dt
    void MoveParticle(array_1d<double, 3>& position,
                      PFEM_Particle_Fluid& pparticle,
                      Element::Pointer& pelement,
                      typename BinBasedFastPointLocator<TDim>::ResultIteratorType& result_begin,
                      const unsigned int& max_results,
                      const array_1d<double, 3>& mesh_displacement,
                      const bool& discriminate_streamlines,
                      Vector& N,
                      const double& delta_t,
                      const array_1d<double, 3>& gravity)
    {

        unsigned int nsubsteps;
        double substep_dt;

        bool KEEP_INTEGRATING = true;
        bool is_found;
        //bool have_air_node;
        //bool have_water_node;

        array_1d<double, 3> vel = ZeroVector(3);
        array_1d<double, 3> vel_without_other_phase_nodes = ZeroVector(3);

        const float particle_distance = pparticle.GetDistance();
        array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        array_1d<double, 3> last_useful_vel;
        double sum_Ns_without_other_phase_nodes;
        //bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
        double only_integral = 0.0;

        const Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in

        vel_without_other_phase_nodes = ZeroVector(3);
        sum_Ns_without_other_phase_nodes = 0.0;
        //distance=0.0;

        if (particle_distance < 0.0 && discriminate_streamlines == true)
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
                {
                    sum_Ns_without_other_phase_nodes += N[j];
                    noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                }
            }

            if (sum_Ns_without_other_phase_nodes > 0.01)
            {
                vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                //flying_water_particle=false;
            }
            else
            {
                vel = particle_velocity;
            }
        } 
        else // air particle or we are not following streamlines
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
            }
            //flying_water_particle=false;
        }

        //calculating substep to get +- courant(substep) = 0.1
        nsubsteps = 10.0 * (delta_t * pelement->GetValue(VELOCITY_OVER_ELEM_SIZE));
        if (nsubsteps < 1)
            nsubsteps = 1;
        substep_dt = delta_t / double(nsubsteps);

        only_integral = 1.0; // weight;//*double(nsubsteps);

        position += vel * substep_dt; //weight;

        ///*****
        last_useful_vel = vel;
        ///*****
        //DONE THE FIRST LOCATION OF THE PARTICLE, NOW WE PROCEED TO STREAMLINE INTEGRATION USING THE MESH VELOCITY
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        unsigned int check_from_element_number = 0;

        for (unsigned int i = 0; i < (nsubsteps - 1); i++) // this is for the substeps n+1. in the first one we already knew the position of the particle.
        {
            // is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
            is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
            if (is_found == true)
            {
                const Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                sum_Ns_without_other_phase_nodes = 0.0;

                if (particle_distance < 0.0 && discriminate_streamlines == true)
                {
                    vel_without_other_phase_nodes = ZeroVector(3);

                    for (unsigned int j = 0; j < TDim + 1; j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                        }
                    }

                    if (sum_Ns_without_other_phase_nodes > 0.01)
                    {
                        vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        //flying_water_particle=false;
                    }
                    else
                    {
                        particle_velocity += substep_dt * gravity;
                        vel = particle_velocity;
                    }
                }
                else //air particle or we are not discriminating streamlines
                {
                    vel_without_other_phase_nodes = ZeroVector(3);
                    vel = ZeroVector(3);
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                    }

                    //flying_water_particle=false;
                }

                only_integral += 1.0; //values saved for the current time step

                position += vel * substep_dt; //weight;
            }
            else
            {
                KEEP_INTEGRATING = false;
                break;
            }
        }
        if (KEEP_INTEGRATING == false)
        {
            pparticle.GetEraseFlag() = true;
        }
        else
        {
            // is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results);
            is_found = mpSearchStructure->FindPointOnMesh(position, N, pelement, result_begin, max_results);
            if (is_found == false)
                pparticle.GetEraseFlag() = true;
        }

    }
/*
    void MoveParticle_Nitsche(PFEM_Particle_Fluid &pparticle,
                              Element::Pointer &pelement,
                              ResultIteratorType result_begin,
                              const unsigned int max_results,
                              const array_1d<double, 3> mesh_displacement,
                              const bool discriminate_streamlines,
                              array_1d<double, TDim + 1> N,
                              double &delta_t,
                              const array_1d<double, 3> &gravity)
    {

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        unsigned int nsubsteps;
        double substep_dt;

        bool KEEP_INTEGRATING = false;
        bool is_found, interface_element;
        //bool have_air_node;
        //bool have_water_node;

        array_1d<double, 3> vel = ZeroVector(3);
        array_1d<double, 3> vel_without_other_phase_nodes = ZeroVector(3);
        array_1d<double, 3> position;
        array_1d<double, 3> mid_position;

        //we start with the first position, then it will enter the loop.
        position = pparticle.Coordinates(); //initial coordinates

        const float particle_distance = pparticle.GetDistance();
        array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        //double distance=0.0;
        array_1d<double, 3> last_useful_vel;
        double sum_Ns_without_other_phase_nodes;
        //double pressure=0.0;
        ///*****
        //bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
        double only_integral = 0.0;

        is_found = FindPositionUsingBins(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
        if (is_found == true)
        {
            interface_element = CheckIfInterfaceElement(pelement);
            KEEP_INTEGRATING = true;
            Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in

            vel_without_other_phase_nodes = ZeroVector(3);
            sum_Ns_without_other_phase_nodes = 0.0;
            //distance=0.0;

            if (interface_element) //the particle is inside an interface element
            {
                for (unsigned int j = 0; j < (TDim + 1); j++)
                {
                    if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE) * N[j];
                    else
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                }
            }
            else //the particle is inside a non-interface element
            {
                if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                        }
                        if (sum_Ns_without_other_phase_nodes > 0.01)
                        {
                            vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        }
                        else
                        {
                            vel = particle_velocity;
                        }
                    }
                }
                else // air particle or we are not following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                    }
                }
            }

            //calculating substep to get +- courant(substep) = 0.1
            nsubsteps = 10.0 * (delta_t * pelement->GetValue(VELOCITY_OVER_ELEM_SIZE));
            if (nsubsteps < 1)
                nsubsteps = 1;
            substep_dt = delta_t / double(nsubsteps);

            only_integral = 1.0; // weight;//*double(nsubsteps);

            position += vel * substep_dt; //weight;

            ///*****
            last_useful_vel = vel;
            ///*****

            //DONE THE FIRST LOCATION OF THE PARTICLE, NOW WE PROCEED TO STREAMLINE INTEGRATION USING THE MESH VELOCITY
            //////////////////////////////////////////////////////////////////////////////////////////////////////
            unsigned int check_from_element_number = 0;

            for (unsigned int i = 0; i < (nsubsteps - 1); i++) // this is for the substeps n+1. in the first one we already knew the position of the particle.
            {
                if (KEEP_INTEGRATING == true)
                {
                    is_found = FindPositionUsingBins(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
                    if (is_found == true)
                    {
                        interface_element = CheckIfInterfaceElement(pelement);
                        Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                        sum_Ns_without_other_phase_nodes = 0.0;
                        vel_without_other_phase_nodes = ZeroVector(3);
                        vel = ZeroVector(3);

                        if (interface_element) //the particle is inside an interface element
                        {
                            for (unsigned int j = 0; j < (TDim + 1); j++)
                            {
                                if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                                    noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE) * N[j];
                                else
                                    noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                            }
                        }
                        else //the particle is inside a non-interface element
                        {
                            if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                            {
                                for (unsigned int j = 0; j < (TDim + 1); j++)
                                {
                                    if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                                    {
                                        sum_Ns_without_other_phase_nodes += N[j];
                                        noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                                    }
                                    if (sum_Ns_without_other_phase_nodes > 0.01)
                                    {
                                        vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                                    }
                                    else
                                    {
                                        particle_velocity += substep_dt * gravity;
                                        vel = particle_velocity;
                                    }
                                }
                            }
                            else // air particle or we are not following streamlines
                            {
                                for (unsigned int j = 0; j < (TDim + 1); j++)
                                {
                                    noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                                }
                            }
                        }

                        only_integral += 1.0; //values saved for the current time step

                        position += vel * substep_dt; //weight;
                    }
                    else
                    {
                        KEEP_INTEGRATING = false;
                        break;
                    }
                }
                else
                    break;
            }
        }

        //if there's a mesh velocity, we add it at the end in a single step:
        position -= mesh_displacement;

        if (KEEP_INTEGRATING == false)
        {
            pparticle.GetEraseFlag() = true;
        }
        else
        {
            is_found = FindPositionUsingBins(position, N, pelement, result_begin, max_results); //we must save the pointer of the last element that we're in (inside the pointervector pelement)
            if (is_found == false)
                pparticle.GetEraseFlag() = true;
        }
        pparticle.Coordinates() = position;
    }

    void MoveParticlePimpleIterativeRK4Fullstep(PFEM_Particle_Fluid &pparticle,
                                                Element::Pointer &pelement,
                                                ResultIteratorType result_begin,
                                                const unsigned int max_results,
                                                const array_1d<double, 3> mesh_displacement,
                                                const bool discriminate_streamlines,
                                                array_1d<double, TDim + 1> N,
                                                double &delta_t,
                                                const array_1d<double, 3> &gravity)
    {

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        unsigned int nsubsteps;
        double substep_dt;

        bool KEEP_INTEGRATING = false;
        bool is_found, is_found1, is_found2, is_found3, is_found4;
        //bool have_air_node;
        //bool have_water_node;

        array_1d<double, 3> vel = ZeroVector(3);
        array_1d<double, 3> vel_without_other_phase_nodes = ZeroVector(3);
        array_1d<double, 3> position;
        array_1d<double, 3> k1 = ZeroVector(3);
        array_1d<double, 3> k2 = ZeroVector(3);
        array_1d<double, 3> k2_aux = ZeroVector(3);
        array_1d<double, 3> k3 = ZeroVector(3);
        array_1d<double, 3> k3_aux = ZeroVector(3);
        array_1d<double, 3> k4 = ZeroVector(3);
        array_1d<double, 3> k4_aux = ZeroVector(3);
        double alphatau = 0.0;

        //we start with the first position, then it will enter the loop.
        position = pparticle.Coordinates(); //initial coordinates

        const float particle_distance = pparticle.GetDistance();
        array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        //double distance=0.0;
        double sum_Ns_without_other_phase_nodes;
        //double pressure=0.0;
        ///*****
        //bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
        //RK4 first step
        is_found1 = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
        if (is_found1 == true)
        {
            KEEP_INTEGRATING = true;
            Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
            vel_without_other_phase_nodes = ZeroVector(3);
            sum_Ns_without_other_phase_nodes = 0.0;
            alphatau = 0.0;
            //distance=0.0;
            if (particle_distance < 0.0 && discriminate_streamlines == true)
            {
                for (unsigned int j = 0; j < (TDim + 1); j++)
                {
                    if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //ok. useful info!
                    {
                        sum_Ns_without_other_phase_nodes += N[j];
                        noalias(vel_without_other_phase_nodes) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                      if (use_mesh_velocity_to_convect)
                        //                          noalias(vel_without_other_phase_nodes) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }

                    noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                    //                  if (use_mesh_velocity_to_convect)
                    //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                }

                if (sum_Ns_without_other_phase_nodes > 0.01)
                {
                    vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                    //flying_water_particle=false;
                }
                else
                {
                    vel = particle_velocity;
                    //                  if (use_mesh_velocity_to_convect)
                    //                  {
                    //                      for(unsigned int j=0; j<(TDim+1); j++)
                    //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    //                  }
                }
            }
            else // air particle or we are not following streamlines
            {
                for (unsigned int j = 0; j < (TDim + 1); j++)
                {
                    noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                    //                  if (use_mesh_velocity_to_convect)
                    //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                }
                //flying_water_particle=false;
            }

            k1 = vel * delta_t; //weight;
            k2_aux = position + k1 / 2.0;

            //RK4 second step
            is_found2 = FindNodeOnMesh(k2_aux, N, pelement, result_begin, max_results);
            if (is_found2 == true)
            {
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                alphatau = delta_t / 2.0;
                //distance=0.0;

                if (particle_distance < 0.0 && discriminate_streamlines == true)
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //ok. useful info!
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                            //                      if (use_mesh_velocity_to_convect)
                            //                          noalias(vel_without_other_phase_nodes) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        }

                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }

                    if (sum_Ns_without_other_phase_nodes > 0.01)
                    {
                        vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        //flying_water_particle=false;
                    }
                    else
                    {
                        vel = particle_velocity;
                        //                  if (use_mesh_velocity_to_convect)
                        //                  {
                        //                      for(unsigned int j=0; j<(TDim+1); j++)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        //                  }
                    }
                }
                else // air particle or we are not following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }
                    //flying_water_particle=false;
                }
            }

            k2 = vel * delta_t; //weight;
            k3_aux = position + k2 / 2.0;

            //RK4 third step
            is_found3 = FindNodeOnMesh(k3_aux, N, pelement, result_begin, max_results);
            if (is_found3 == true)
            {
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                //distance=0.0;

                if (particle_distance < 0.0 && discriminate_streamlines == true)
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //ok. useful info!
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                            //                      if (use_mesh_velocity_to_convect)
                            //                          noalias(vel_without_other_phase_nodes) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        }

                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }

                    if (sum_Ns_without_other_phase_nodes > 0.01)
                    {
                        vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        //flying_water_particle=false;
                    }
                    else
                    {
                        vel = particle_velocity;
                        //                  if (use_mesh_velocity_to_convect)
                        //                  {
                        //                      for(unsigned int j=0; j<(TDim+1); j++)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        //                  }
                    }
                }
                else // air particle or we are not following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }
                    //flying_water_particle=false;
                }
            }

            k3 = vel * delta_t; //weight;
            k4_aux = position + k3;

            //RK4 fourth step
            is_found4 = FindNodeOnMesh(k4_aux, N, pelement, result_begin, max_results);
            if (is_found4 == true)
            {
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                alphatau = delta_t;
                //distance=0.0;

                if (particle_distance < 0.0 && discriminate_streamlines == true)
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //ok. useful info!
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                            //                      if (use_mesh_velocity_to_convect)
                            //                          noalias(vel_without_other_phase_nodes) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        }

                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }

                    if (sum_Ns_without_other_phase_nodes > 0.01)
                    {
                        vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        //flying_water_particle=false;
                    }
                    else
                    {
                        vel = particle_velocity;
                        //                  if (use_mesh_velocity_to_convect)
                        //                  {
                        //                      for(unsigned int j=0; j<(TDim+1); j++)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                        //                  }
                    }
                }
                else // air particle or we are not following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += ((1.0 - alphatau) * geom[j].FastGetSolutionStepValue(VELOCITY) * N[j] + (alphatau)*geom[j].FastGetSolutionStepValue(VELOCITY_PIMPLE) * N[j]);
                        //                  if (use_mesh_velocity_to_convect)
                        //                          noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY)*N[j];
                    }
                    //flying_water_particle=false;
                }
            }

            k4 = vel * delta_t; //weight;
        }

        if (is_found1 == true && is_found2 == true && is_found3 == true && is_found4 == true)
        {
            position += (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        else
            KEEP_INTEGRATING = false;

        //if there's a mesh velocity, we add it at the end in a single step:
        position -= mesh_displacement;

        if (KEEP_INTEGRATING == false)
        {
            pparticle.GetEraseFlag() = true;
        }
        else
        {
            is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //we must save the pointer of the last element that we're in (inside the pointervector pelement)
            if (is_found == false)
                pparticle.GetEraseFlag() = true;
        }

        pparticle.Coordinates() = position;
    }
*/
    // void MoveParticleRK4(PFEM_Particle_Fluid &pparticle,
    //                      Element::Pointer &pelement,
    //                      ResultIteratorType result_begin,
    //                      const unsigned int max_results,
    //                      const array_1d<double, 3> mesh_displacement,
    //                      const bool discriminate_streamlines,
    //                      array_1d<double, TDim + 1> &N,
    //                      double &delta_t,
    //                      const array_1d<double, 3> &gravity)
    // {

        
    //     ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

    //     unsigned int nsubsteps;
    //     double substep_dt;

    //     bool KEEP_INTEGRATING = true;
    //     bool is_found, is_found1, is_found2, is_found3, is_found4;
    //     //bool have_air_node;
    //     //bool have_water_node;

    //     array_1d<double, 3> vel = ZeroVector(3);
    //     array_1d<double, 3> vel_without_other_phase_nodes = ZeroVector(3);
    //     array_1d<double, 3> position;
    //     array_1d<double, 3> k1 = ZeroVector(3);
    //     array_1d<double, 3> k2 = ZeroVector(3);
    //     array_1d<double, 3> k2_aux = ZeroVector(3);
    //     array_1d<double, 3> k3 = ZeroVector(3);
    //     array_1d<double, 3> k3_aux = ZeroVector(3);
    //     array_1d<double, 3> k4 = ZeroVector(3);
    //     array_1d<double, 3> k4_aux = ZeroVector(3);

    //     //we start with the first position, then it will enter the loop.
    //     position = pparticle.Coordinates(); //initial coordinates

    //     const float particle_distance = pparticle.GetDistance();
    //     array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
    //     //double distance=0.0;
    //     double sum_Ns_without_other_phase_nodes;
    //     //double pressure=0.0;
    //     ///*****
    //     //bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
    //     //RK4 first step
    //     // is_found1 = FindNodeOnMesh(position, N ,pelement,result_begin,max_results); //good, now we know where this point is:

    //     Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
    //     vel_without_other_phase_nodes = ZeroVector(3);
    //     sum_Ns_without_other_phase_nodes = 0.0;
    //     //distance=0.0;
    //     if (particle_distance < 0.0 && discriminate_streamlines == true)
    //     {
    //         for (unsigned int j = 0; j < (TDim + 1); j++)
    //         {
    //             if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
    //             {
    //                 sum_Ns_without_other_phase_nodes += N[j];
    //                 noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //             }
    //         }

    //         if (sum_Ns_without_other_phase_nodes > 0.01)
    //         {
    //             vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
    //             //flying_water_particle=false;
    //         }
    //         else
    //         {
    //             vel = particle_velocity;
    //         }
    //     }
    //     else // air particle or we are not following streamlines
    //     {
    //         for (unsigned int j = 0; j < (TDim + 1); j++)
    //         {
    //             noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //         }
    //         //flying_water_particle=false;
    //     }

    //     k1 = vel * delta_t; //weight;
    //     k2_aux = position + k1 / 2.0;

    //     //RK4 second step
    //     is_found2 = FindNodeOnMesh(k2_aux, N, pelement, result_begin, max_results);
    //     if (is_found2 == true)
    //     {
    //         Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
    //         vel = ZeroVector(3);
    //         vel_without_other_phase_nodes = ZeroVector(3);
    //         sum_Ns_without_other_phase_nodes = 0.0;
    //         //distance=0.0;

    //         if (particle_distance < 0.0 && discriminate_streamlines == true)
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
    //                 {
    //                     sum_Ns_without_other_phase_nodes += N[j];
    //                     noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //                 }

    //             }

    //             if (sum_Ns_without_other_phase_nodes > 0.01)
    //             {
    //                 vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
    //                 //flying_water_particle=false;
    //             }
    //             else
    //             {
    //                 vel = particle_velocity;
    //             }
    //         }
    //         else // air particle or we are not following streamlines
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //             }
    //             //flying_water_particle=false;
    //         }
        

    //         k2 = vel * delta_t; //weight;
    //         k3_aux = position + k2 / 2.0;

        

    //      //RK4 third step
    //      is_found3 = FindNodeOnMesh(k3_aux, N, pelement, result_begin, max_results);
    //      if (is_found3 == true)
    //      {
    //         Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
    //         vel = ZeroVector(3);
    //         vel_without_other_phase_nodes = ZeroVector(3);
    //         sum_Ns_without_other_phase_nodes = 0.0;
    //         //distance=0.0;

    //         if (particle_distance < 0.0 && discriminate_streamlines == true)
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
    //                 {
    //                     sum_Ns_without_other_phase_nodes += N[j];
    //                     noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //                 }
    //             }

    //             if (sum_Ns_without_other_phase_nodes > 0.01)
    //             {
    //                 vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
    //                 //flying_water_particle=false;
    //             }
    //             else
    //             {
    //                 vel = particle_velocity;
    //             }
    //         }
    //         else // air particle or we are not following streamlines
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //             }
    //             //flying_water_particle=false;
    //         }
        

    //       k3 = vel * delta_t; //weight;
    //       k4_aux = position + k3;

    //        //RK4 fourth step
    //        is_found4 = FindNodeOnMesh(k4_aux, N, pelement, result_begin, max_results);
    //        if (is_found4 == true)
    //        {
    //           Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
    //           vel = ZeroVector(3);
    //           vel_without_other_phase_nodes = ZeroVector(3);
    //           sum_Ns_without_other_phase_nodes = 0.0;
    //           //distance=0.0;

    //           if (particle_distance < 0.0 && discriminate_streamlines == true)
    //           {
    //               for (unsigned int j = 0; j < (TDim + 1); j++)
    //               {
    //                   if ((geom[j].FastGetSolutionStepValue(DISTANCE,1)) < 0.0) //ok. useful info!
    //                   {
    //                       sum_Ns_without_other_phase_nodes += N[j];
    //                       noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //                   }
    //               }

    //               if (sum_Ns_without_other_phase_nodes > 0.01)
    //               {
    //                   vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
    //                   //flying_water_particle=false;
    //               }
    //               else
    //               {
    //                   vel = particle_velocity;
    //               }
    //           }
    //           else // air particle or we are not following streamlines
    //           {
    //               for (unsigned int j = 0; j < (TDim + 1); j++)
    //               {
    //                   noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
    //               }
    //               //flying_water_particle=false;
    //           }
    //           k4 = vel * delta_t; //weight;
    //        }
    //      }
    //     }
    //     if (is_found2 == true && is_found3 == true && is_found4 == true)
    //     {
    //         position += (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    //     }
    //     else
    //         KEEP_INTEGRATING = false;


    //     if (KEEP_INTEGRATING == false)
    //     {
    //         pparticle.GetEraseFlag() = true;
    //     }
    //     else
    //     {
    //         is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //we must save the pointer of the last element that we're in (inside the pointervector pelement)
    //         if (is_found == false)
    //             pparticle.GetEraseFlag() = true;
    //     }

    //     pparticle.Coordinates() = position;
    // }
/*
    void MoveParticleRK4_Nitsche(PFEM_Particle_Fluid &pparticle,
                                 Element::Pointer &pelement,
                                 ResultIteratorType result_begin,
                                 const unsigned int max_results,
                                 const array_1d<double, 3> mesh_displacement,
                                 const bool discriminate_streamlines,
                                 array_1d<double, TDim + 1> N,
                                 double &delta_t,
                                 const array_1d<double, 3> &gravity)
    {

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();

        unsigned int nsubsteps;
        double substep_dt;

        bool KEEP_INTEGRATING = false;
        bool is_found, is_found1, is_found2, is_found3, is_found4, interface_element;
        //bool have_air_node;
        //bool have_water_node;

        array_1d<double, 3> vel = ZeroVector(3);
        array_1d<double, 3> vel_without_other_phase_nodes = ZeroVector(3);
        array_1d<double, 3> position;
        array_1d<double, 3> k1 = ZeroVector(3);
        array_1d<double, 3> k2 = ZeroVector(3);
        array_1d<double, 3> k2_aux = ZeroVector(3);
        array_1d<double, 3> k3 = ZeroVector(3);
        array_1d<double, 3> k3_aux = ZeroVector(3);
        array_1d<double, 3> k4 = ZeroVector(3);
        array_1d<double, 3> k4_aux = ZeroVector(3);
        double alphatau = 0.0;

        //we start with the first position, then it will enter the loop.
        position = pparticle.Coordinates(); //initial coordinates

        const float particle_distance = pparticle.GetDistance();
        array_1d<float, 3> particle_velocity = pparticle.GetVelocity();
        //double distance=0.0;
        double sum_Ns_without_other_phase_nodes;
        //double pressure=0.0;
        ///*****
        //bool flying_water_particle=true; //if a water particle does not find a water element in its whole path, then we add the gravity*dt
        //RK4 first step
        //check if the particle is inside an interface element
        interface_element = CheckIfInterfaceElement(pelement);
        is_found1 = FindPositionUsingBins(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
        if (is_found1 == true)
        {
            interface_element = CheckIfInterfaceElement(pelement);
            KEEP_INTEGRATING = true;
            Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
            vel_without_other_phase_nodes = ZeroVector(3);
            sum_Ns_without_other_phase_nodes = 0.0;
            alphatau = 0.0;
            //distance=0.0;
            if (interface_element) //the particle is inside an interface element
            {
                for (unsigned int j = 0; j < (TDim + 1); j++)
                {
                    if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE, 1) * N[j];
                    else
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                }
            }
            else //the particle is inside a non-interface element
            {
                if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                        {
                            sum_Ns_without_other_phase_nodes += N[j];
                            noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                        }
                        if (sum_Ns_without_other_phase_nodes > 0.01)
                        {
                            vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                        }
                        else
                        {
                            vel = particle_velocity;
                        }
                    }
                }
                else // air particle or we are not following streamlines
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                    }
                }
            }

            k1 = vel * delta_t; //weight;
            k2_aux = position + k1 / 2.0;

            //RK4 second step
            is_found2 = FindPositionUsingBins(k2_aux, N, pelement, result_begin, max_results);
            if (is_found2 == true)
            {
                interface_element = CheckIfInterfaceElement(pelement);
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                alphatau = delta_t / 2.0;
                //distance=0.0;
                if (interface_element) //the particle is inside an interface element
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE, 1) * N[j];
                        else
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                    }
                }
                else //the particle is inside a non-interface element
                {
                    if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                            {
                                sum_Ns_without_other_phase_nodes += N[j];
                                noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                            }
                            if (sum_Ns_without_other_phase_nodes > 0.01)
                            {
                                vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                            }
                            else
                            {
                                particle_velocity += (delta_t / 2.0) * gravity;
                                vel = particle_velocity;
                            }
                        }
                    }
                    else // air particle or we are not following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                        }
                    }
                }
            }

            k2 = vel * delta_t; //weight;
            k3_aux = position + k2 / 2.0;

            //RK4 third step
            is_found3 = FindPositionUsingBins(k3_aux, N, pelement, result_begin, max_results);
            if (is_found3 == true)
            {
                interface_element = CheckIfInterfaceElement(pelement);
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                //distance=0.0;
                if (interface_element)
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE, 1) * N[j];
                        else
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                    }
                }
                else //the particle is inside a non-interface element
                {
                    if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                            {
                                sum_Ns_without_other_phase_nodes += N[j];
                                noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                            }
                            if (sum_Ns_without_other_phase_nodes > 0.01)
                            {
                                vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                            }
                            else
                            {
                                vel = particle_velocity;
                            }
                        }
                    }
                    else // air particle or we are not following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                        }
                    }
                }
            }

            k3 = vel * delta_t; //weight;
            k4_aux = position + k3;

            //RK4 fourth step
            is_found4 = FindPositionUsingBins(k4_aux, N, pelement, result_begin, max_results);
            if (is_found4 == true)
            {
                interface_element = CheckIfInterfaceElement(pelement);
                Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
                vel = ZeroVector(3);
                vel_without_other_phase_nodes = ZeroVector(3);
                sum_Ns_without_other_phase_nodes = 0.0;
                alphatau = delta_t;
                //distance=0.0;

                if (interface_element)
                {
                    for (unsigned int j = 0; j < (TDim + 1); j++)
                    {
                        if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY_NITSCHE, 1) * N[j];
                        else
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                    }
                }
                else //the particle is inside a non-interface element
                {
                    if (particle_distance < 0.0 && discriminate_streamlines == true) // water particle and following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0) //if the node is negative distance
                            {
                                sum_Ns_without_other_phase_nodes += N[j];
                                noalias(vel_without_other_phase_nodes) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                            }
                            if (sum_Ns_without_other_phase_nodes > 0.01)
                            {
                                vel = vel_without_other_phase_nodes / sum_Ns_without_other_phase_nodes;
                            }
                            else
                            {
                                particle_velocity += (delta_t / 2.0) * gravity;
                                vel = particle_velocity;
                            }
                        }
                    }
                    else // air particle or we are not following streamlines
                    {
                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY, 1) * N[j];
                        }
                    }
                }
            }

            k4 = vel * delta_t; //weight;
        }

        if (is_found1 == true && is_found2 == true && is_found3 == true && is_found4 == true)
        {
            position += (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        else
            KEEP_INTEGRATING = false;

        //if there's a mesh velocity, we add it at the end in a single step:
        position -= mesh_displacement;

        if (KEEP_INTEGRATING == false)
        {
            pparticle.GetEraseFlag() = true;
        }
        else
        {
            is_found = FindPositionUsingBins(position, N, pelement, result_begin, max_results); //we must save the pointer of the last element that we're in (inside the pointervector pelement)
            if (is_found == false)
                pparticle.GetEraseFlag() = true;
        }

        pparticle.Coordinates() = position;
    }
*/
    void AccelerateParticleUsingDeltaVelocity(
        PFEM_Particle_Fluid& pparticle,
        Element::Pointer& pelement,
        Geometry<Node<3>>& geom,
        const double& delta_t,
        const array_1d<double, 3>& gravity,
        Vector& N)
    {
        //we start with the first position, then it will enter the loop.
        array_1d<double, 3> coords = pparticle.Coordinates();
        float &particle_distance = pparticle.GetDistance();
        //double distance=0.0;
        array_1d<double, 3> delta_velocity = ZeroVector(3);

        array_1d<double, 3> delta_velocity_without_air = ZeroVector(3);
        array_1d<double, 3> delta_velocity_without_water = ZeroVector(3);
        double sum_Ns_without_water_nodes = 0.0;
        double sum_Ns_without_air_nodes = 0.0;
        if (particle_distance > 0.0) //no problem. air
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                //just for air
                if ((geom[j].FastGetSolutionStepValue(DISTANCE)) > 0.0)
                {
                    noalias(delta_velocity_without_water) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
                    sum_Ns_without_air_nodes += N[j];
                }
                //both air and water
                noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
            }

            if (sum_Ns_without_water_nodes > 0.01)
            {
                //delta_velocity = delta_velocity_without_water/sum_Ns_without_water_nodes ; //commented = using all the velocities always!
            }
            //else we use the complete field
        }
        else //water particle
        {

            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0)
                {
                    noalias(delta_velocity_without_air) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
                    sum_Ns_without_air_nodes += N[j];
                }
                noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
            }

            if (sum_Ns_without_air_nodes > 0.01)
            {
                delta_velocity = delta_velocity_without_air / sum_Ns_without_air_nodes;
            }
            else
            {
                if (mDENSITY_WATER > (10.0 * mDENSITY_AIR))
                {
                    delta_velocity = gravity * (1.0 - mDENSITY_AIR / mDENSITY_WATER) * delta_t;
                }
            }
        }
        pparticle.GetVelocity() = pparticle.GetVelocity() + delta_velocity;
    }

    // void AccelerateParticleUsingDeltaVelocity_Nitsche(
    //     PFEM_Particle_Fluid &pparticle,
    //     Element::Pointer &pelement,
    //     Geometry<Node<3>> &geom,
    //     const double &delta_t,
    //     const array_1d<double, 3> &gravity,
    //     array_1d<double, TDim + 1> &N)
    // {
    //     //      array_1d<double,TDim+1> N;

    //     ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
    //     //      const double delta_t = CurrentProcessInfo[DELTA_TIME];
    //     //      array_1d<double,3> gravity = CurrentProcessInfo[GRAVITY];

    //     //we start with the first position, then it will enter the loop.
    //     array_1d<double, 3> coords = pparticle.Coordinates();
    //     float &particle_distance = pparticle.GetDistance();
    //     //double distance=0.0;
    //     array_1d<double, 3> delta_velocity = ZeroVector(3);

    //     array_1d<double, 3> delta_velocity_without_air = ZeroVector(3);
    //     array_1d<double, 3> delta_velocity_without_water = ZeroVector(3);
    //     double sum_Ns_without_water_nodes = 0.0;
    //     double sum_Ns_without_air_nodes = 0.0;
    //     bool interface_element = CheckIfInterfaceElement(pelement);
    //     //      bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
    //     //      if(is_found == false)
    //     //      {
    //     //          KRATOS_WATCH(N)
    //     //          for (int j=0 ; j!=(TDim+1); j++)
    //     //                              if (N[j]<0.0 )
    //     //                                  N[j]=1e-10;
    //     //      }

    //     if (particle_distance > 0.0) //no problem. air
    //     {
    //         if (interface_element) //the particle is inside an interface element
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
    //                     noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY_NITSCHE) * N[j];
    //                 else
    //                     noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //             }
    //         }
    //         else
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 //just for air
    //                 if ((geom[j].FastGetSolutionStepValue(DISTANCE)) > 0.0)
    //                 {
    //                     noalias(delta_velocity_without_water) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //                     sum_Ns_without_air_nodes += N[j];
    //                 }
    //                 //both air and water
    //                 noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //             }

    //             if (sum_Ns_without_water_nodes > 0.01)
    //             {
    //                 //delta_velocity = delta_velocity_without_water/sum_Ns_without_water_nodes ; //commented = using all the velocities always!
    //             }
    //             //else we use the complete field
    //         }
    //     }
    //     else //water particle
    //     {
    //         if (interface_element) //the particle is inside an interface element
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 if (particle_distance * geom[j].FastGetSolutionStepValue(DISTANCE) < 0.0) //opposite signs->nitsche_dofs
    //                     noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY_NITSCHE) * N[j];
    //                 else
    //                     noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //             }
    //         }
    //         else
    //         {
    //             for (unsigned int j = 0; j < (TDim + 1); j++)
    //             {
    //                 if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0)
    //                 {
    //                     noalias(delta_velocity_without_air) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //                     sum_Ns_without_air_nodes += N[j];
    //                 }
    //                 noalias(delta_velocity) += geom[j].FastGetSolutionStepValue(DELTA_VELOCITY) * N[j];
    //             }

    //             if (sum_Ns_without_air_nodes > 0.01)
    //             {
    //                 delta_velocity = delta_velocity_without_air / sum_Ns_without_air_nodes;
    //             }
    //             else
    //             {
    //                 if (mDENSITY_WATER > (10.0 * mDENSITY_AIR))
    //                 {
    //                     delta_velocity = gravity * (1.0 - mDENSITY_AIR / mDENSITY_WATER) * delta_t;
    //                 }
    //             }
    //         }
    //     }
    //     pparticle.GetVelocity() = pparticle.GetVelocity() + delta_velocity;
    // }

    void AccelerateParticleUsingMeshDeltaVelocity(
        PFEM_Particle_Fluid &pparticle,
        Element::Pointer &pelement,
        Geometry<Node<3>> &geom,
        const double &delta_t,
        const array_1d<double, 3> &gravity,
        array_1d<double, TDim + 1> &N)
    {
        //         array_1d<double,TDim+1> N;

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        //         const double delta_t = CurrentProcessInfo[DELTA_TIME];
        //         array_1d<double,3> gravity = CurrentProcessInfo[GRAVITY];

        //we start with the first position, then it will enter the loop.
        array_1d<double, 3> coords = pparticle.Coordinates();
        float &particle_distance = pparticle.GetDistance();
        //double distance=0.0;
        array_1d<double, 3> mesh_delta_velocity = ZeroVector(3);

        array_1d<double, 3> mesh_delta_velocity_without_air = ZeroVector(3);
        array_1d<double, 3> mesh_delta_velocity_without_water = ZeroVector(3);
        double sum_Ns_without_water_nodes = 0.0;
        double sum_Ns_without_air_nodes = 0.0;

        //         bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
        //         if(is_found == false)
        //         {
        //             KRATOS_WATCH(N)
        //             for (int j=0 ; j!=(TDim+1); j++)
        //                                 if (N[j]<0.0 )
        //                                     N[j]=1e-10;
        //         }

        if (particle_distance > 0.0) //no problem. air
        {
            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                //just for air
                if ((geom[j].FastGetSolutionStepValue(DISTANCE)) > 0.0)
                {
                    noalias(mesh_delta_velocity_without_water) += geom[j].FastGetSolutionStepValue(MESH_DELTA_VELOCITY) * N[j];
                    sum_Ns_without_air_nodes += N[j];
                }
                //both air and water
                noalias(mesh_delta_velocity) += geom[j].FastGetSolutionStepValue(MESH_DELTA_VELOCITY) * N[j];
            }

            if (sum_Ns_without_water_nodes > 0.01)
            {
                //mesh_delta_velocity = mesh_delta_velocity_without_water/sum_Ns_without_water_nodes ; //commented = using all the velocities always!
            }
            //else we use the complete field
        }
        else //water particle
        {

            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                if ((geom[j].FastGetSolutionStepValue(DISTANCE)) < 0.0)
                {
                    noalias(mesh_delta_velocity_without_air) += geom[j].FastGetSolutionStepValue(MESH_DELTA_VELOCITY) * N[j];
                    sum_Ns_without_air_nodes += N[j];
                }
                noalias(mesh_delta_velocity) += geom[j].FastGetSolutionStepValue(MESH_DELTA_VELOCITY) * N[j];
            }

            if (sum_Ns_without_air_nodes > 0.01)
            {
                mesh_delta_velocity = mesh_delta_velocity_without_air / sum_Ns_without_air_nodes;
            }
            else
            {
                if (mDENSITY_WATER > (10.0 * mDENSITY_AIR))
                {
                    mesh_delta_velocity = gravity * (1.0 - mDENSITY_AIR / mDENSITY_WATER) * delta_t;
                }
            }
        }
        pparticle.GetVelocity() = pparticle.GetVelocity() + mesh_delta_velocity;
    }

    void MoveParticle_inverse_way(
        PFEM_Particle_Fluid &pparticle,
        Element::Pointer &pelement, //NOT A REFERENCE!! WE SHALL NOT OVERWRITE THE ELEMENT IT BELONGS TO!
        ResultIteratorType result_begin,
        const unsigned int max_results,
        const bool use_mesh_velocity_to_convect)
    {

        ProcessInfo &CurrentProcessInfo = mr_model_part.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        unsigned int nsubsteps;
        double substep_dt;

        bool KEEP_INTEGRATING = false;
        bool is_found;

        array_1d<double, 3> vel;
        array_1d<double, 3> particle_vel;
        array_1d<double, 3> position;
        array_1d<double, 3> mid_position;
        Vector N(TDim + 1);

        //we start with the first position, then it will enter the loop.
        position = pparticle.Coordinates(); // + (pparticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates

        float &distance = pparticle.GetDistance();
        double only_integral = 0.0;

        is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
        if (is_found == true)
        {
            KEEP_INTEGRATING = true;
            Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in
            vel = ZeroVector(3);
            particle_vel = ZeroVector(3);
            distance = 0.0;

            for (unsigned int j = 0; j < (TDim + 1); j++)
            {
                distance += geom[j].FastGetSolutionStepValue(DISTANCE) * N(j);
                noalias(particle_vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                if (use_mesh_velocity_to_convect)
                    noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY) * N[j];
            }
            //calculating substep to get +- courant(substep) = 1/4
            nsubsteps = 10.0 * (delta_t * pelement->GetValue(VELOCITY_OVER_ELEM_SIZE));
            if (nsubsteps < 1)
                nsubsteps = 1;
            substep_dt = delta_t / double(nsubsteps);

            only_integral = 1.0;          // weight;//*double(nsubsteps);
            position -= vel * substep_dt; //weight;

            for (unsigned int i = 0; i < (nsubsteps - 1); i++) // this is for the substeps n+1. in the first one we already knew the position of the particle.
            {
                if (KEEP_INTEGRATING == true)
                {
                    is_found = FindNodeOnMesh(position, N, pelement, result_begin, max_results); //good, now we know where this point is:
                    if (is_found == true)
                    {
                        Geometry<Node<3>> &geom = pelement->GetGeometry(); //the element we're in

                        vel = ZeroVector(3);
                        particle_vel = ZeroVector(3);
                        distance = 0.0;

                        for (unsigned int j = 0; j < (TDim + 1); j++)
                        {
                            noalias(particle_vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY) * N[j];
                            distance += geom[j].FastGetSolutionStepValue(DISTANCE) * N(j);
                            if (use_mesh_velocity_to_convect)
                                noalias(vel) -= geom[j].FastGetSolutionStepValue(MESH_VELOCITY) * N[j];
                        }

                        only_integral += 1.0;         //weight ; //values saved for the current time step
                        position -= vel * substep_dt; //weight;
                    }
                    else
                        KEEP_INTEGRATING = false;
                }
            }

            ///COMMENT TO GET A A CONTINOUS DISTANCE FUNCTION FIELD!!!!!
            if (distance > 0.0)
            {
                //if(distance<2.0)
                distance = 1.0;
                //else
                //  distance=3.0;
            }
            else
                distance = -1.0;

            pparticle.GetVelocity() = particle_vel;
        }
        //else {KRATOS_WATCH(position); }
    }

    void OverwriteParticleDataUsingTopographicDomain(
        PFEM_Particle_Fluid &pparticle,
        Element::Pointer &pelement,
        array_1d<double, 3> domains_offset,
        ResultIteratorType result_begin,
        const unsigned int max_results)
    {
        array_1d<double, TDim + 1> N;

        //we start with the first position, then it will enter the loop.
        array_1d<double, 3> coords = pparticle.Coordinates() + domains_offset;
        float &particle_distance = pparticle.GetDistance();
        bool is_found = FindNodeOnTopographicMesh(coords, N, pelement, result_begin, max_results); //good, now we know where this point is:

        if (is_found) //it is part of the solid topographic domain
        {
            particle_distance = -1.0;
        }
        else //it is outside the topographic domain, therefore it is air or whatever it means
        {
            particle_distance = 1.0;
        }

        pparticle.GetVelocity() = ZeroVector(3);
    }

    ///this function should find the element into which a given node is located
    ///and return a pointer to the element and the vector containing the
    ///shape functions that define the postion within the element
    ///if "false" is devolved the element is not found
    // bool FindNodeOnMeshFullBins(
    //     const array_1d<double, 3 >& rCoordinates,
    //     Vector& rNShapeFunction,
    //     typename EntityType::Pointer& pEntity,
    //     ResultIteratorType ItResultBegin,
    //     const SizeType max_results = 1000,
    //     const double Tolerance = 1.0e-5
    //     )
    // {
    //     // Ask to the container for the list of candidate entities
    //     const SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(typename BinsObjectDynamic<ConfigureType>::PointType{rCoordinates}, ItResultBegin, max_results);

    //     if (results_found > 0) {
    //         // Loop over the candidate entities and check if the particle falls within
    //         for (IndexType i = 0; i < static_cast<IndexType>(results_found); i++) {

    //             GeometryType& geom = (*(ItResultBegin + i))->GetGeometry();

    //             // Find local position
    //             array_1d<double, 3> point_local_coordinates;
    //             const bool is_found = geom.IsInside(rCoordinates, point_local_coordinates, Tolerance);
    //             geom.ShapeFunctionsValues(rNShapeFunction, point_local_coordinates);

    //             if (is_found) {
    //                 pEntity = (*(ItResultBegin + i));
    //                 return true;
    //             }
    //         }
    //     }

    //     // Not found case
    //     pEntity = nullptr;
    //     return false;
    // }

    bool FindNodeOnMesh(const array_1d<double, 3 >& rCoordinates,
                        Vector& rNShapeFunction,
                        typename EntityType::Pointer& pEntity,
                        ResultIteratorType result_begin,
                        const unsigned int max_results)
    {
        typedef std::size_t SizeType;

        //before using the bin to search for possible elements we check first the last element in which the particle was.
        Geometry<Node<3>> &geom_default = pEntity->GetGeometry(); //(*(i))->GetGeometry();
        bool is_found_1 = CalculatePosition(geom_default, rCoordinates[0], rCoordinates[1], rCoordinates[2], rNShapeFunction);
        if (is_found_1 == true) //that was easy!
        {
            return true;
        }

        //to begin with we check the neighbour elements; it is a bit more expensive
        GlobalPointersVector<Element> &neighb_elems = pEntity->GetValue(NEIGHBOUR_ELEMENTS);
        //the first we check is the one that has negative shape function, because it means it went outside in this direction:
        //commented, it is not faster than simply checking all the neighbours (branching) 
        //we check all the neighbour elements
        for (unsigned int i = 0; i != (neighb_elems.size()); i++)
        {

            Geometry<Node<3>> &geom = neighb_elems[i].GetGeometry();
            bool is_found_2 = CalculatePosition(geom, rCoordinates[0], rCoordinates[1], rCoordinates[2], rNShapeFunction);
            if (is_found_2)
            {
                pEntity = neighb_elems[i].shared_from_this();
                return true;
            }
        }

        //if checking all the neighbour elements did not work, we have to use the bins
        //ask to the container for the list of candidate elements
        bool is_found_3 = mpSearchStructure->FindPointOnMesh(rCoordinates, rNShapeFunction, pEntity, result_begin, max_results);
        if (is_found_3)
        { 
            return true;
        }

        return false;
    }

    bool CheckIfInterfaceElement(Element::Pointer &pelement)
    {

        array_1d<double, 3> distances;
        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < 3; i++)
        {
            distances[i] = pelement->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            if (distances[i] > 0.0)
                npos++;
            if (distances[i] < 0.0)
                nneg++;
        }
        if (npos == 3.0 or nneg == 3.0)
            return false;
        else
            return true;
    }

    bool CheckIfInterfaceElement(ModelPart::ElementsContainerType::iterator &pelement)
    {

        array_1d<double, 3> distances;
        unsigned int npos = 0, nneg = 0;
        for (unsigned int i = 0; i < 3; i++)
        {
            distances[i] = pelement->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            if (distances[i] > 0.0)
                npos++;
            if (distances[i] < 0.0)
                nneg++;
        }
        if (npos == 3.0 or nneg == 3.0)
            return false;
        else
            return true;
    }

    // // VERSION INCLUDING PREDEFINED ELEMENTS FOLLOWING A TRAJECTORY
    // bool FindNodeOnMesh(array_1d<double, 3> &position,
    //                     array_1d<double, TDim + 1> &N,
    //                     Element::Pointer &pelement,
    //                     GlobalPointersVector<Element> &elements_in_trajectory,
    //                     unsigned int &number_of_elements_in_trajectory,
    //                     unsigned int &check_from_element_number,
    //                     ResultIteratorType result_begin,
    //                     const unsigned int max_results)
    // {
    //     typedef std::size_t SizeType;

    //     const array_1d<double, 3> &coords = position;
    //     array_1d<double, TDim + 1> aux_N;
    //     //before using the bin to search for possible elements we check first the last element in which the particle was.
    //     Geometry<Node<3>> &geom_default = pelement->GetGeometry(); //(*(i))->GetGeometry();
    //     bool is_found_1 = CalculatePosition(geom_default, coords[0], coords[1], coords[2], N);
    //     if (is_found_1 == true)
    //     {
    //         return true; //that was easy!
    //     }

    //     //if it was not found in the first element, we can proceed to check in the following elements (in the trajectory defined by previous particles that started from the same element.
    //     for (unsigned int i = (check_from_element_number); i != number_of_elements_in_trajectory; i++)
    //     {
    //         Geometry<Node<3>> &geom = elements_in_trajectory[i].GetGeometry();
    //         bool is_found_2 = CalculatePosition(geom, coords[0], coords[1], coords[2], aux_N);
    //         if (is_found_2)
    //         {
    //             pelement = elements_in_trajectory[i].shared_from_this();
    //             N = aux_N;
    //             check_from_element_number = i + 1; //now i element matches pelement, so to avoid cheching twice the same element we send the counter to the following element.
    //             return true;
    //         }
    //     }

    //     //now we check the neighbour elements:
    //     GlobalPointersVector<Element> &neighb_elems = pelement->GetValue(NEIGHBOUR_ELEMENTS);
    //     //the first we check is the one that has negative shape function, because it means it went outside in this direction:
    //     //commented, it is not faster than simply checking all the neighbours (branching)
    //     /*
    //     unsigned int checked_element=0;
    //     for (unsigned int i=0;i!=(TDim+1);i++)
    //     {
    //         if (N[i]<0.0)
    //         {
    //             checked_element=i;
    //             Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
    //             bool is_found_2 = CalculatePosition(geom,coords[0],coords[1],coords[2],aux_N);
    //             if (is_found_2)
    //             {
    //                 pelement=Element::Pointer(((neighb_elems(i))));
    //                 N=aux_N;
    //                 return true;
    //             }
    //             break;
    //         }
    //     }
    //     */
    //     //we check all the neighbour elements
    //     for (unsigned int i = 0; i != (neighb_elems.size()); i++)
    //     {

    //         Geometry<Node<3>> &geom = neighb_elems[i].GetGeometry();
    //         bool is_found_2 = CalculatePosition(geom, coords[0], coords[1], coords[2], N);
    //         if (is_found_2)
    //         {
    //             pelement = neighb_elems[i].shared_from_this();
    //             if (number_of_elements_in_trajectory < 20)
    //             {
    //                 elements_in_trajectory(number_of_elements_in_trajectory) = pelement;
    //                 number_of_elements_in_trajectory++;
    //                 check_from_element_number = number_of_elements_in_trajectory; //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the elements_in_trajectory list. we are the particle that is adding elements to the list
    //             }
    //             return true;
    //         }
    //     }

    //     //if checking all the neighbour elements did not work, we have to use the bins
    //     //ask to the container for the list of candidate elements
    //     SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(Point{coords}, result_begin, max_results);

    //     if (results_found > 0)
    //     {
    //         //loop over the candidate elements and check if the particle falls within
    //         for (SizeType i = 0; i < results_found; i++)
    //         {
    //             Geometry<Node<3>> &geom = (*(result_begin + i))->GetGeometry();

    //             //find local position
    //             bool is_found = CalculatePosition(geom, coords[0], coords[1], coords[2], N);

    //             if (is_found == true)
    //             {
    //                 pelement = Element::Pointer((*(result_begin + i)));
    //                 if (number_of_elements_in_trajectory < 20)
    //                 {
    //                     elements_in_trajectory(number_of_elements_in_trajectory) = pelement;
    //                     number_of_elements_in_trajectory++;
    //                     check_from_element_number = number_of_elements_in_trajectory; //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the elements_in_trajectory list. we are the particle that is adding elements to the list
    //                 }
    //                 return true;
    //             }
    //         }
    //     }

    //     //not found case
    //     return false;
    // }

    ///this function should find the element into which a given node is located
    ///and return a pointer to the element and the vector containing the
    ///shape functions that define the postion within the element
    ///if "false" is devolved the element is not found
    bool FindNodeOnTopographicMesh(array_1d<double, 3> &position,
                                   array_1d<double, TDim + 1> &N,
                                   Element::Pointer &pelement,
                                   ResultIteratorType result_begin,
                                   const unsigned int max_results)
    {
        typedef std::size_t SizeType;

        const array_1d<double, 3> &coords = position;
        array_1d<double, TDim + 1> aux_N;
        //before using the bin to search for possible elements we check first the last element in which the particle was.

        //ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin()+last_element;
        Geometry<Node<3>> &geom_default = pelement->GetGeometry(); //(*(i))->GetGeometry();
        bool is_found_1 = CalculatePosition(geom_default, coords[0], coords[1], coords[2], N);
        if (is_found_1 == true)
        {
            //pelement = (*(i));
            return true;
        }

        //to begin with we check the neighbour elements:
        GlobalPointersVector<Element> &neighb_elems = pelement->GetValue(NEIGHBOUR_ELEMENTS);
        for (unsigned int i = 0; i != (neighb_elems.size()); i++)
        {

            Geometry<Node<3>> &geom = neighb_elems[i].GetGeometry();
            bool is_found_2 = CalculatePosition(geom, coords[0], coords[1], coords[2], N);
            if (is_found_2)
            {
                pelement = neighb_elems[i].shared_from_this();
                return true;
            }
        }

        //ask to the container for the list of candidate elements
        SizeType results_found = mpTopographicBinsObjectDynamic->SearchObjectsInCell(Point{coords}, result_begin, max_results);
        //KRATOS_WATCH(results_found)

        if (results_found > 0)
        {
            //loop over the candidate elements and check if the particle falls within
            for (SizeType i = 0; i < results_found; i++)
            {
                Geometry<Node<3>> &geom = (*(result_begin + i))->GetGeometry();

                //find local position
                bool is_found = CalculatePosition(geom, coords[0], coords[1], coords[2], N);

                if (is_found == true)
                {
                    pelement = Element::Pointer((*(result_begin + i)));
                    return true;
                }
            }
        }

        //not found case
        return false;
    }

    //***************************************
    //***************************************

    inline bool CalculatePosition(Geometry<Node<3>> &geom,
                                  const double xc, const double yc, const double zc,
                                  Vector& N)
    {

        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;
        if (area == 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
        //KRATOS_WATCH(N);

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }
    ////////////
    //using the pre loaded nodal coordinates
    inline bool CalculatePosition(const array_1d<double, 3 * (TDim + 1)> &nodes_positions,
                                  const double xc, const double yc, const double zc,
                                  Vector& N)
    {
        const double &x0 = nodes_positions[0];
        const double &y0 = nodes_positions[1];
        const double &x1 = nodes_positions[3];
        const double &y1 = nodes_positions[4];
        const double &x2 = nodes_positions[6];
        const double &y2 = nodes_positions[7];

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;
        if (area == 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
        //KRATOS_WATCH(N);

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }

    //***************************************
    //***************************************

    // inline bool CalculatePosition(Geometry<Node<3>> &geom,
    //                               const double xc, const double yc, const double zc,
    //                               array_1d<double, 4> &N)
    // {

    //     double x0 = geom[0].X();
    //     double y0 = geom[0].Y();
    //     double z0 = geom[0].Z();
    //     double x1 = geom[1].X();
    //     double y1 = geom[1].Y();
    //     double z1 = geom[1].Z();
    //     double x2 = geom[2].X();
    //     double y2 = geom[2].Y();
    //     double z2 = geom[2].Z();
    //     double x3 = geom[3].X();
    //     double y3 = geom[3].Y();
    //     double z3 = geom[3].Z();

    //     double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    //     double inv_vol = 0.0;
    //     if (vol < 0.000000000000000000000000000001)
    //     {
    //         KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
    //     }
    //     else
    //     {
    //         inv_vol = 1.0 / vol;
    //     }

    //     N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
    //     N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
    //     N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
    //     N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;

    //     if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
    //         N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
    //         //if the xc yc zc is inside the tetrahedron return true
    //         return true;

    //     return false;
    // }
    // ///////////////////
    // //using the pre loaded nodal coordinates
    // inline bool CalculatePosition(const array_1d<double, 3 * (TDim + 1)> &nodes_positions,
    //                               const double xc, const double yc, const double zc,
    //                               array_1d<double, 4> &N)
    // {

    //     const double &x0 = nodes_positions[0];
    //     const double &y0 = nodes_positions[1];
    //     const double &z0 = nodes_positions[2];
    //     const double &x1 = nodes_positions[3];
    //     const double &y1 = nodes_positions[4];
    //     const double &z1 = nodes_positions[5];
    //     const double &x2 = nodes_positions[6];
    //     const double &y2 = nodes_positions[7];
    //     const double &z2 = nodes_positions[8];
    //     const double &x3 = nodes_positions[9];
    //     const double &y3 = nodes_positions[10];
    //     const double &z3 = nodes_positions[11];

    //     double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

    //     double inv_vol = 0.0;
    //     if (vol < 0.000000000000000000000000000001)
    //     {
    //         KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
    //     }
    //     else
    //     {
    //         inv_vol = 1.0 / vol;
    //     }

    //     N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
    //     N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
    //     N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
    //     N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;

    //     if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
    //         N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
    //         //if the xc yc zc is inside the tetrahedron return true
    //         return true;

    //     return false;
    // }

    inline double CalculateVol(const double x0, const double y0,
                               const double x1, const double y1,
                               const double x2, const double y2)
    {
        return 0.5 * ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));
    }
    //***************************************
    //***************************************

    inline double CalculateVol(const double x0, const double y0, const double z0,
                               const double x1, const double y1, const double z1,
                               const double x2, const double y2, const double z2,
                               const double x3, const double y3, const double z3)
    {
        double x10 = x1 - x0;
        double y10 = y1 - y0;
        double z10 = z1 - z0;

        double x20 = x2 - x0;
        double y20 = y2 - y0;
        double z20 = z2 - z0;

        double x30 = x3 - x0;
        double y30 = y3 - y0;
        double z30 = z3 - z0;

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return detJ * 0.1666666666666666666667;
    }

    void ComputeGaussPointPositions_4(Geometry<Node<3>> &geom, BoundedMatrix<double, 7, 3> &pos, BoundedMatrix<double, 7, 3> &N)
    {
        double one_third = 1.0 / 3.0;
        double one_sixt = 0.15; //1.0 / 6.0;
        double two_third = 0.7; //2.0 * one_third;

        N(0, 0) = one_sixt;
        N(0, 1) = one_sixt;
        N(0, 2) = two_third;
        N(1, 0) = two_third;
        N(1, 1) = one_sixt;
        N(1, 2) = one_sixt;
        N(2, 0) = one_sixt;
        N(2, 1) = two_third;
        N(2, 2) = one_sixt;
        N(3, 0) = one_third;
        N(3, 1) = one_third;
        N(3, 2) = one_third;

        //first
        pos(0, 0) = one_sixt * geom[0].X() + one_sixt * geom[1].X() + two_third * geom[2].X();
        pos(0, 1) = one_sixt * geom[0].Y() + one_sixt * geom[1].Y() + two_third * geom[2].Y();
        pos(0, 2) = one_sixt * geom[0].Z() + one_sixt * geom[1].Z() + two_third * geom[2].Z();

        //second
        pos(1, 0) = two_third * geom[0].X() + one_sixt * geom[1].X() + one_sixt * geom[2].X();
        pos(1, 1) = two_third * geom[0].Y() + one_sixt * geom[1].Y() + one_sixt * geom[2].Y();
        pos(1, 2) = two_third * geom[0].Z() + one_sixt * geom[1].Z() + one_sixt * geom[2].Z();

        //third
        pos(2, 0) = one_sixt * geom[0].X() + two_third * geom[1].X() + one_sixt * geom[2].X();
        pos(2, 1) = one_sixt * geom[0].Y() + two_third * geom[1].Y() + one_sixt * geom[2].Y();
        pos(2, 2) = one_sixt * geom[0].Z() + two_third * geom[1].Z() + one_sixt * geom[2].Z();

        //fourth
        pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
        pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
        pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();
    }

    void ComputeGaussPointPositionsForPostReseed(Geometry<Node<3>> &geom, BoundedMatrix<double, 7, 3> &pos, BoundedMatrix<double, 7, 3> &N) //2d
    {
        double one_third = 1.0 / 3.0;
        double one_eight = 0.12;      //1.0 / 6.0;
        double three_quarters = 0.76; //2.0 * one_third;

        N(0, 0) = one_eight;
        N(0, 1) = one_eight;
        N(0, 2) = three_quarters;

        N(1, 0) = three_quarters;
        N(1, 1) = one_eight;
        N(1, 2) = one_eight;

        N(2, 0) = one_eight;
        N(2, 1) = three_quarters;
        N(2, 2) = one_eight;

        N(3, 0) = one_third;
        N(3, 1) = one_third;
        N(3, 2) = one_third;

        N(4, 0) = one_eight;
        N(4, 1) = 0.44;
        N(4, 2) = 0.44;

        N(5, 0) = 0.44;
        N(5, 1) = one_eight;
        N(5, 2) = 0.44;

        N(6, 0) = 0.44;
        N(6, 1) = 0.44;
        N(6, 2) = one_eight;

        //first
        pos(0, 0) = one_eight * geom[0].X() + one_eight * geom[1].X() + three_quarters * geom[2].X();
        pos(0, 1) = one_eight * geom[0].Y() + one_eight * geom[1].Y() + three_quarters * geom[2].Y();
        pos(0, 2) = one_eight * geom[0].Z() + one_eight * geom[1].Z() + three_quarters * geom[2].Z();

        //second
        pos(1, 0) = three_quarters * geom[0].X() + one_eight * geom[1].X() + one_eight * geom[2].X();
        pos(1, 1) = three_quarters * geom[0].Y() + one_eight * geom[1].Y() + one_eight * geom[2].Y();
        pos(1, 2) = three_quarters * geom[0].Z() + one_eight * geom[1].Z() + one_eight * geom[2].Z();

        //third
        pos(2, 0) = one_eight * geom[0].X() + three_quarters * geom[1].X() + one_eight * geom[2].X();
        pos(2, 1) = one_eight * geom[0].Y() + three_quarters * geom[1].Y() + one_eight * geom[2].Y();
        pos(2, 2) = one_eight * geom[0].Z() + three_quarters * geom[1].Z() + one_eight * geom[2].Z();

        //fourth
        pos(3, 0) = one_third * geom[0].X() + one_third * geom[1].X() + one_third * geom[2].X();
        pos(3, 1) = one_third * geom[0].Y() + one_third * geom[1].Y() + one_third * geom[2].Y();
        pos(3, 2) = one_third * geom[0].Z() + one_third * geom[1].Z() + one_third * geom[2].Z();

        //fifth
        pos(4, 0) = one_eight * geom[0].X() + 0.44 * geom[1].X() + 0.44 * geom[2].X();
        pos(4, 1) = one_eight * geom[0].Y() + 0.44 * geom[1].Y() + 0.44 * geom[2].Y();
        pos(4, 2) = one_eight * geom[0].Z() + 0.44 * geom[1].Z() + 0.44 * geom[2].Z();

        //sixth
        pos(5, 0) = 0.44 * geom[0].X() + one_eight * geom[1].X() + 0.44 * geom[2].X();
        pos(5, 1) = 0.44 * geom[0].Y() + one_eight * geom[1].Y() + 0.44 * geom[2].Y();
        pos(5, 2) = 0.44 * geom[0].Z() + one_eight * geom[1].Z() + 0.44 * geom[2].Z();

        //seventh
        pos(6, 0) = 0.44 * geom[0].X() + 0.44 * geom[1].X() + one_eight * geom[2].X();
        pos(6, 1) = 0.44 * geom[0].Y() + 0.44 * geom[1].Y() + one_eight * geom[2].Y();
        pos(6, 2) = 0.44 * geom[0].Z() + 0.44 * geom[1].Z() + one_eight * geom[2].Z();
    }

    void ComputeGaussPointPositionsForPostReseed(Geometry<Node<3>> &geom, BoundedMatrix<double, 9, 3> &pos, BoundedMatrix<double, 9, 4> &N) //3D
    {
        double one_quarter = 0.25;
        double small_fraction = 0.1; //1.0 / 6.0;
        double big_fraction = 0.7;   //2.0 * one_third;
        double mid_fraction = 0.3;   //2.0 * one_third;

        N(0, 0) = big_fraction;
        N(0, 1) = small_fraction;
        N(0, 2) = small_fraction;
        N(0, 3) = small_fraction;

        N(1, 0) = small_fraction;
        N(1, 1) = big_fraction;
        N(1, 2) = small_fraction;
        N(1, 3) = small_fraction;

        N(2, 0) = small_fraction;
        N(2, 1) = small_fraction;
        N(2, 2) = big_fraction;
        N(2, 3) = small_fraction;

        N(3, 0) = small_fraction;
        N(3, 1) = small_fraction;
        N(3, 2) = small_fraction;
        N(3, 3) = big_fraction;

        N(4, 0) = one_quarter;
        N(4, 1) = one_quarter;
        N(4, 2) = one_quarter;
        N(4, 3) = one_quarter;

        N(5, 0) = small_fraction;
        N(5, 1) = mid_fraction;
        N(5, 2) = mid_fraction;
        N(5, 3) = mid_fraction;

        N(6, 0) = mid_fraction;
        N(6, 1) = small_fraction;
        N(6, 2) = mid_fraction;
        N(6, 3) = mid_fraction;

        N(7, 0) = mid_fraction;
        N(7, 1) = mid_fraction;
        N(7, 2) = small_fraction;
        N(7, 3) = mid_fraction;

        N(8, 0) = mid_fraction;
        N(8, 1) = mid_fraction;
        N(8, 2) = mid_fraction;
        N(8, 3) = small_fraction;

        pos = ZeroMatrix(9, 3);
        for (unsigned int i = 0; i != 4; i++) //going through the 4 nodes
        {
            array_1d<double, 3> &coordinates = geom[i].Coordinates();
            for (unsigned int j = 0; j != 9; j++) //going through the 9 particles
            {
                for (unsigned int k = 0; k != 3; k++) //x,y,z
                    pos(j, k) += N(j, i) * coordinates[k];
            }
        }
    }

    void ComputeGaussPointPositionsForPreReseed(Geometry<Node<3>> &geom, BoundedMatrix<double, 3, 3> &pos, BoundedMatrix<double, 3, 3> &N) //2D
    {

        N(0, 0) = 0.5;
        N(0, 1) = 0.25;
        N(0, 2) = 0.25;

        N(1, 0) = 0.25;
        N(1, 1) = 0.5;
        N(1, 2) = 0.25;

        N(2, 0) = 0.25;
        N(2, 1) = 0.25;
        N(2, 2) = 0.5;

        //first
        pos(0, 0) = 0.5 * geom[0].X() + 0.25 * geom[1].X() + 0.25 * geom[2].X();
        pos(0, 1) = 0.5 * geom[0].Y() + 0.25 * geom[1].Y() + 0.25 * geom[2].Y();
        pos(0, 2) = 0.5 * geom[0].Z() + 0.25 * geom[1].Z() + 0.25 * geom[2].Z();

        //second
        pos(1, 0) = 0.25 * geom[0].X() + 0.5 * geom[1].X() + 0.25 * geom[2].X();
        pos(1, 1) = 0.25 * geom[0].Y() + 0.5 * geom[1].Y() + 0.25 * geom[2].Y();
        pos(1, 2) = 0.25 * geom[0].Z() + 0.5 * geom[1].Z() + 0.25 * geom[2].Z();

        //third
        pos(2, 0) = 0.25 * geom[0].X() + 0.25 * geom[1].X() + 0.5 * geom[2].X();
        pos(2, 1) = 0.25 * geom[0].Y() + 0.25 * geom[1].Y() + 0.5 * geom[2].Y();
        pos(2, 2) = 0.25 * geom[0].Z() + 0.25 * geom[1].Z() + 0.5 * geom[2].Z();
    }

    void ComputeGaussPointPositionsForPreReseed(Geometry<Node<3>> &geom, BoundedMatrix<double, 4, 3> &pos, BoundedMatrix<double, 4, 4> &N) //3D
    {
        //creating 4 particles, each will be closer to a node and equidistant to the other nodes

        N(0, 0) = 0.4;
        N(0, 1) = 0.2;
        N(0, 2) = 0.2;
        N(0, 3) = 0.2;

        N(1, 0) = 0.2;
        N(1, 1) = 0.4;
        N(1, 2) = 0.2;
        N(1, 3) = 0.2;

        N(2, 0) = 0.2;
        N(2, 1) = 0.2;
        N(2, 2) = 0.4;
        N(2, 3) = 0.2;

        N(3, 0) = 0.2;
        N(3, 1) = 0.2;
        N(3, 2) = 0.2;
        N(3, 3) = 0.4;

        pos = ZeroMatrix(4, 3);
        for (unsigned int i = 0; i != 4; i++) //going through the 4 nodes
        {
            array_1d<double, 3> &coordinates = geom[i].Coordinates();
            for (unsigned int j = 0; j != 4; j++) //going through the 4 particles
            {
                for (unsigned int k = 0; k != 3; k++) //x,y,z
                    pos(j, k) += N(j, i) * coordinates[k];
            }
        }
    }

    void ComputeGaussPointPositions_45(Geometry<Node<3>> &geom, BoundedMatrix<double, 45, 3> &pos, BoundedMatrix<double, 45, 3> &N) //2D
    {
        //std::cout << "NEW ELEMENT" << std::endl;
        unsigned int counter = 0;
        for (unsigned int i = 0; i != 9; i++)
        {
            for (unsigned int j = 0; j != (9 - i); j++)
            {
                N(counter, 0) = 0.05 + double(i) * 0.1;
                N(counter, 1) = 0.05 + double(j) * 0.1;
                N(counter, 2) = 1.0 - (N(counter, 1) + N(counter, 0));
                pos(counter, 0) = N(counter, 0) * geom[0].X() + N(counter, 1) * geom[1].X() + N(counter, 2) * geom[2].X();
                pos(counter, 1) = N(counter, 0) * geom[0].Y() + N(counter, 1) * geom[1].Y() + N(counter, 2) * geom[2].Y();
                pos(counter, 2) = N(counter, 0) * geom[0].Z() + N(counter, 1) * geom[1].Z() + N(counter, 2) * geom[2].Z();
                //std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
                counter++;
            }
        }
    }

    void ComputeGaussPointPositions_45(Geometry<Node<3>> &geom, BoundedMatrix<double, 60, 3> &pos, BoundedMatrix<double, 60, 4> &N) //3D
    {
        //std::cout << "NEW ELEMENT" << std::endl;
        KRATOS_THROW_ERROR(std::logic_error, "3d_45nodes is not complete", "");
        //             unsigned int counter=0;
        //             for (unsigned int i=0; i!=9;i++)
        //             {
        //                 for (unsigned int j=0; j!=(9-i);j++)
        //                 {
        //                     N(counter,0)=0.05+double(i)*0.1;
        //                     N(counter,1)=0.05+double(j)*0.1;
        //                     N(counter,2)=1.0 - ( N(counter,1)+ N(counter,0) ) ;
        //                     pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X();
        //                     pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y();
        //                     pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z();
        //                     //std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
        //                     counter++;
        //
        //                 }
        //             }
    }

    void ComputeGaussPointPositions_initial(Geometry<Node<3>> &geom, BoundedMatrix<double, 15, 3> &pos, BoundedMatrix<double, 15, 3> &N) //2D
    {
        //std::cout << "NEW ELEMENT" << std::endl;
        unsigned int counter = 0;
        for (unsigned int i = 0; i != 5; i++)
        {
            for (unsigned int j = 0; j != (5 - i); j++)
            {
                N(counter, 0) = 0.05 + double(i) * 0.2;
                N(counter, 1) = 0.05 + double(j) * 0.2;
                N(counter, 2) = 1.0 - (N(counter, 1) + N(counter, 0));
                pos(counter, 0) = N(counter, 0) * geom[0].X() + N(counter, 1) * geom[1].X() + N(counter, 2) * geom[2].X();
                pos(counter, 1) = N(counter, 0) * geom[0].Y() + N(counter, 1) * geom[1].Y() + N(counter, 2) * geom[2].Y();
                pos(counter, 2) = N(counter, 0) * geom[0].Z() + N(counter, 1) * geom[1].Z() + N(counter, 2) * geom[2].Z();
                //std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
                counter++;
            }
        }
    }

    void ComputeGaussPointPositions_initial(Geometry<Node<3>> &geom, BoundedMatrix<double, 20, 3> &pos, BoundedMatrix<double, 20, 4> &N) //3D
    {
        //std::cout << "NEW ELEMENT" << std::endl;
        //double total;
        double fraction_increment;
        unsigned int counter = 0;
        for (unsigned int i = 0; i != 4; i++) //going to build a particle "pyramid"(tetrahedra) by layers. the first layer will be made by a triangle of 4 base X 4 height. since it is a triangle, it means it will have 10 particles
        {
            //std::cout << "inside i" <<  i << std::endl;
            for (unsigned int j = 0; j != (4 - i); j++)
            {
                //std::cout << "inside j" << j << std::endl;
                for (unsigned int k = 0; k != (4 - i - j); k++)
                {
                    //std::cout << "inside k" << k << std::endl;
                    N(counter, 0) = 0.27 * (0.175 + double(i)); //this is our "surface" in which we will build each layer, so we must construct a triangle using what's left of the shape functions total (a total of 1)

                    //total = 1.0 - N(counter,0);
                    fraction_increment = 0.27; //

                    N(counter, 1) = fraction_increment * (0.175 + double(j));
                    N(counter, 2) = fraction_increment * (0.175 + double(k));
                    N(counter, 3) = 1.0 - (N(counter, 0) + N(counter, 1) + N(counter, 2));
                    pos(counter, 0) = N(counter, 0) * geom[0].X() + N(counter, 1) * geom[1].X() + N(counter, 2) * geom[2].X() + N(counter, 3) * geom[3].X();
                    pos(counter, 1) = N(counter, 0) * geom[0].Y() + N(counter, 1) * geom[1].Y() + N(counter, 2) * geom[2].Y() + N(counter, 3) * geom[3].Y();
                    pos(counter, 2) = N(counter, 0) * geom[0].Z() + N(counter, 1) * geom[1].Z() + N(counter, 2) * geom[2].Z() + N(counter, 3) * geom[3].Z();
                    //std::cout << N(counter,0) << " " << N(counter,1) << " " << N(counter,2) << " " << std::endl;
                    counter++;
                }
            }
        }
    }

    bool FindPositionUsingBins(array_1d<double, 3> &position,
                               array_1d<double, TDim + 1> &N,
                               Element::Pointer &pelement,
                               ResultIteratorType result_begin,
                               const unsigned int max_results)
    {
        typedef std::size_t SizeType;
        const array_1d<double, 3> &coords = position;
        ResultContainerType results_a(10000); //DWARNING-previously it was 10, but that was causing problems
        ResultIteratorType result_begin_a = results_a.begin();

        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(Point{coords}, result_begin_a, max_results);
        if (results_found > 0)
        {
            for (SizeType i = 0; i < results_found; i++)
            {
                Geometry<Node<3>> &geom = (*(result_begin_a + i))->GetGeometry();
                bool is_found = CalculatePosition(geom, coords[0], coords[1], coords[2], N);
                if (is_found == true)
                {
                    pelement = Element::Pointer((*(result_begin_a + i)));
                    return true;
                }
            }
        }

        return false;
    }

    // Bubble Sort Function for Descending Order
    void BubbleSort(array_1d<double, 7> &distances, array_1d<int, 7> &positions, unsigned int &arrange_number)
    {
        int i, j;
        bool flag = true; // set flag to 1 to start first pass
        double temp;      // holding variable
        int temp_position;
        int numLength = arrange_number;
        for (i = 1; (i <= numLength) && flag; i++)
        {
            flag = false;
            for (j = 0; j < (numLength - 1); j++)
            {
                if (distances[j + 1] < distances[j]) // descending order simply changes to >
                {
                    temp = distances[j]; // swap elements
                    distances[j] = distances[j + 1];
                    distances[j + 1] = temp;

                    temp_position = positions[j]; //swap positions
                    positions[j] = positions[j + 1];
                    positions[j + 1] = temp_position;

                    flag = true; // indicates that a swap occurred.
                }
            }
        }
        return; //arrays are passed to functions by address; nothing is returned
    }

    void BubbleSort(array_1d<double, 9> &distances, array_1d<int, 9> &positions, unsigned int &arrange_number)
    {
        int i, j;
        bool flag = true; // set flag to 1 to start first pass
        double temp;      // holding variable
        int temp_position;
        int numLength = arrange_number;
        for (i = 1; (i <= numLength) && flag; i++)
        {
            flag = false;
            for (j = 0; j < (numLength - 1); j++)
            {
                if (distances[j + 1] < distances[j]) // descending order simply changes to >
                {
                    temp = distances[j]; // swap elements
                    distances[j] = distances[j + 1];
                    distances[j + 1] = temp;

                    temp_position = positions[j]; //swap positions
                    positions[j] = positions[j + 1];
                    positions[j + 1] = temp_position;

                    flag = true; // indicates that a swap occurred.
                }
            }
        }
        return; //arrays are passed to functions by address; nothing is returned
    }

    template <class T>
    bool InvertMatrix(const T &input, T &inverse)
    {
        typedef permutation_matrix<std::size_t> pmatrix;

        // create a working copy of the input
        T A(input);

        // create a permutation matrix for the LU-factorization
        pmatrix pm(A.size1());

        // perform LU-factorization
        int res = lu_factorize(A, pm);
        if (res != 0)
            return false;

        // create identity matrix of "inverse"
        inverse.assign(identity_matrix<double>(A.size1()));

        // backsubstitute to get the inverse
        lu_substitute(A, pm, inverse);

        return true;
    }

    bool InvertMatrix3x3(const BoundedMatrix<double, TDim + 1, TDim + 1> &A, BoundedMatrix<double, TDim + 1, TDim + 1> &result)
    {
        double determinant = +A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) - A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) + A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
        double invdet = 1 / determinant;
        result(0, 0) = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) * invdet;
        result(1, 0) = -(A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) * invdet;
        result(2, 0) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * invdet;
        result(0, 1) = -(A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) * invdet;
        result(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * invdet;
        result(2, 1) = -(A(0, 0) * A(1, 2) - A(1, 0) * A(0, 2)) * invdet;
        result(0, 2) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)) * invdet;
        result(1, 2) = -(A(0, 0) * A(2, 1) - A(2, 0) * A(0, 1)) * invdet;
        result(2, 2) = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * invdet;

        return true;
    }

    ModelPart &mr_model_part;
    ModelPart *mtopographic_model_part_pointer;
    array_1d<double, 3> mcalculation_domain_complete_displacement;
    array_1d<double, 3> mcalculation_domain_added_displacement;
    bool mintialized_transfer_tool;
    bool muse_mesh_velocity_to_convect;
    int m_nparticles;
    int mnelems;
    double mDENSITY_WATER;
    double mDENSITY_AIR;

    //vector<double> mareas_vector; UNUSED SO COMMENTED
    int max_nsubsteps;
    int mmaximum_number_of_particles;
    std::vector<PFEM_Particle_Fluid> mparticles_vector; //Point<3>
    int mlast_elem_id;
    bool modd_timestep;
    bool mparticle_printing_tool_initialized;
    unsigned int mfilter_factor;
    unsigned int mlast_node_id;
    //ModelPart& mr_particle_model_part;

    // vector<int> mnumber_of_particles_in_elems;
    // vector<int> mnumber_of_particles_in_elems_aux; - perhaps delete
 
    typename BinBasedFastPointLocator<TDim>::Pointer mpSearchStructure;
    typename BinsObjectDynamic<Configure>::Pointer mpBinsObjectDynamic;
    typename BinsObjectDynamic<Configure>::Pointer mpTopographicBinsObjectDynamic;

    void CalculateNormal(Geometry<Node<3>> &pGeometry, array_1d<double, 3> &An);
};

template <>
void MoveParticleUtilityFullBinsPFEM2<2>::CalculateNormal(Geometry<Node<3>> &pGeometry, array_1d<double, 3> &An)
{
    array_1d<double, 2> v1;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();

    An[0] = -v1[1];
    An[1] = v1[0];
    An[2] = 0.0;

    //now checking orientation using the normal:
    const unsigned int NumNodes = 2;
    array_1d<double, 3> nodal_normal = ZeroVector(3);
    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        nodal_normal += pGeometry[iNode].FastGetSolutionStepValue(NORMAL);

    double dot_prod = nodal_normal[0] * An[0] + nodal_normal[1] * An[1];
    if (dot_prod < 0.0)
    {
        //std::cout << "inverting the normal" << std::endl;
        An *= -1.0; // inverting the direction of the normal!!!
    }
}

template <>
void MoveParticleUtilityFullBinsPFEM2<3>::CalculateNormal(Geometry<Node<3>> &pGeometry, array_1d<double, 3> &An)
{
    array_1d<double, 3> v1, v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An, v1, v2);
    An *= 0.5;

    //now checking orientation using the normal:
    const unsigned int NumNodes = 3;
    array_1d<double, 3> nodal_normal = ZeroVector(3);
    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        nodal_normal += pGeometry[iNode].FastGetSolutionStepValue(NORMAL);

    double dot_prod = nodal_normal[0] * An[0] + nodal_normal[1] * An[1] + nodal_normal[2] * An[2];
    if (dot_prod < 0.0)
    {
        //std::cout << "inverting the normal!!" << std::endl;
        An *= -1.0; // inverting the direction of the normal!!!
    }
}

} // namespace Kratos.

#endif // KRATOS_MOVE_PART_UTILITY_DIFF2_INCLUDED  defined
