//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//                   Pablo Becker
//

#if !defined(KRATOS_MOVE_SHALLOW_WATER_PARTICLE_UTILITY_H_INCLUDED)
#define  KRATOS_MOVE_SHALLOW_WATER_PARTICLE_UTILITY_H_INCLUDED

///@defgroup MoveShallowWaterParticleUtility
///@brief Utility to move particles on the eulerian mesh with an
/// explicit scheme. This is the basic tool of the pfem2 framework


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/checks.h"

#include "includes/dof.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "processes/node_erase_process.h"

#include "utilities/geometry_utilities.h"

#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "spatial_containers/spatial_containers.h"
#include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "utilities/spatial_containers_configure.h"

#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/point.h"

#include "shallow_water_application.h"
#include "shallow_water_particle.h"

#include "utilities/openmp_utils.h"

#include "time.h"

//#include "processes/process.h"

namespace Kratos
{

//this class is to be modified by the user to customize the interpolation process
template< unsigned int TDim>
class MoveShallowWaterParticleUtility
{
public:

    typedef SpatialContainersConfigure<TDim>                   Configure;
    typedef typename Configure::PointType                      PointType;
    typedef typename Configure::ContainerType                  ContainerType;
    typedef typename Configure::IteratorType                   IteratorType;
    typedef typename Configure::ResultContainerType            ResultContainerType;
    typedef typename Configure::ResultIteratorType             ResultIteratorType;
    typedef PointerVector< ShallowParticle, ShallowParticle*, std::vector<ShallowParticle*> > ParticlePointerVector;

    KRATOS_CLASS_POINTER_DEFINITION(MoveShallowWaterParticleUtility);

    //template<unsigned int TDim>
    MoveShallowWaterParticleUtility(ModelPart& rModelPart, Parameters rParameters) :
        mrModelPart(rModelPart),
        mScalarVar1(KratosComponents< Variable<double> >::Get( rParameters["convection_scalar_variable"].GetString() ) ),
        mVectorVar1(KratosComponents< Variable<array_1d<double,3> > >::Get( rParameters["convection_vector_variable"].GetString() ) )
    {
        KRATOS_TRY

        std::cout << "Initializing moveparticle utility for scalar transport" << std::endl;

        Parameters default_parameters( R"(
            {
                "convection_scalar_variable"  : "HEIGHT",
                "convection_vector_variable"  : "VELOCITY",
                "maximum_number_of_particles" : 16
            }  )" );

        // Now validate agains defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        m_scalar_var1_name = rParameters["convection_scalar_variable"].GetString();
        m_vector_var1_name = rParameters["convection_vector_variable"].GetString();
        mMaxNumberOfParticles = rParameters["maximum_number_of_particles"].GetDouble();

        Check();
        //storing water and air density and their inverses, just in case it is needed for the streamline integration
        //loop in elements to change their ID to their position in the array. Easier to get information later.
        //DO NOT PARALELIZE THIS! IT MUST BE SERIAL!!!!!!!!!!!!!!!!!!!!!!
        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();
        for(unsigned int  ii=0; ii<mrModelPart.Elements().size(); ii++)
        {
            ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
            ielem->SetId(ii+1);
        }
        mLastElemId= (mrModelPart.ElementsEnd()-1)->Id();
        int node_id=0;
        // we look for the smallest edge. could be used as a weighting function when going lagrangian->eulerian instead of traditional shape functions(method currently used)
        ModelPart::NodesContainerType::iterator inodebegin = mrModelPart.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Nodes().size(), node_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator pnode = inodebegin+ii;
                array_1d<double,3> position_node;
                double distance=0.0;
                position_node = pnode->Coordinates();
                WeakPointerVector< Node<3> >& rneigh = pnode->GetValue(NEIGHBOUR_NODES);
                //we loop all the nodes to check all the edges
                const double number_of_neighbours = static_cast<double>(rneigh.size());
                for( WeakPointerVector<Node<3> >::iterator inode = rneigh.begin(); inode!=rneigh.end(); inode++)
                {
                    array_1d<double,3> position_difference;
                    position_difference = inode->Coordinates() - position_node;
                    const double current_distance = norm_2( position_difference );
                    distance += current_distance / number_of_neighbours;
                }
                //and we save the largest edge.
                pnode->SetValue(MEAN_SIZE, distance);

                node_id=pnode->GetId();
            }
        }
        mLastNodeId=node_id;

        //we also calculate the element mean size in the same way, for the courant number
        //also we set the right size to the LHS column for the pressure enrichments, in order to recover correctly the enrichment pressure
        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Elements().size(), element_partition);

        //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
            {
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;

                double elem_size;
                array_1d<double,3> Edge(3,0.0);
                Edge = ielem->GetGeometry()[1].Coordinates() - ielem->GetGeometry()[0].Coordinates();
                elem_size = Edge[0]*Edge[0];
                for (unsigned int d = 1; d < TDim; d++)
                    elem_size += Edge[d]*Edge[d];

                for (unsigned int i = 2; i < (TDim+1); i++)
                    for(unsigned int j = 0; j < i; j++)
                    {
                        Edge = ielem->GetGeometry()[i].Coordinates() - ielem->GetGeometry()[j].Coordinates();
                        double Length = Edge[0]*Edge[0];
                        for (unsigned int d = 1; d < TDim; d++)
                            Length += Edge[d]*Edge[d];
                        if (Length < elem_size) elem_size = Length;
                    }
                elem_size = sqrt(elem_size);
                ielem->SetValue(MEAN_SIZE, elem_size);
            }
        }

        //matrix containing the position of the 4/15/45 particles that we will seed at the beggining
        BoundedMatrix<double, 5*(1+TDim), 3 > pos;
        BoundedMatrix<double, 5*(1+TDim), (1+TDim) > N;

        int particle_id=0;
        mNElems = mrModelPart.Elements().size();

        std::cout << "  about to resize vectors" << std::endl;

        //setting the right size to the vector containing the particles assigned to each element
        //particles vector. this vector contains ALL the particles in the simulation.
        mParticlesVector.resize(mNElems*mMaxNumberOfParticles);

        //and this vector contains the current number of particles that are in each element (currently zero)
        mNumOfParticlesInElems.resize(mNElems);
        mNumOfParticlesInElems=ZeroVector(mNElems);

        //when moving the particles, an auxiliary vector is necessary (to store the previous number)
        mNumOfParticlesInElemsAux.resize(mNElems);

        //each element will have a list of pointers to all the particles that are inside.
        //this vector contains the pointers to the vector of (particle) pointers of each element.
        mVectorOfParticlePointersVectors.resize(mNElems);
        //int artz;
        //std::cin >> artz;
        int i_int=0; //careful! it's not the id, but the position inside the array!
        std::cout << "  about to create particles" << std::endl;
        //now we seed: LOOP IN ELEMENTS
        //using loop index, DO NOT paralelize this! change lines : mparticles_in_elems_pointers((ii*mMaxNumberOfParticles)+mparticles_in_elems_integers(ii)) = pparticle; and the next one

        mOffset=0;
        //ShallowParticle& firstparticle = mParticlesVector[0];
        for(unsigned int ii=0; ii<mrModelPart.Elements().size(); ii++)
        {
            ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
            //(ielem->GetValue(BED_PARTICLE_POINTERS)) = ParticlePointerVector( mMaxNumberOfParticles*2, &firstparticle );
            //ParticlePointerVector&  particle_pointers =  (ielem->GetValue(BED_PARTICLE_POINTERS));
            //now we link the mpointers_to_particle_pointers_vectors to the corresponding element
            //mpointers_to_particle_pointers_vectors(ii) = &particle_pointers;
            //now we resize the vector of particle pointers. it is double sized because we move the particles from an initial position (first half) to a final position (second half).
            //for(int j=0; j<(mMaxNumberOfParticles*2); j++)
            //        particle_pointers.push_back(&firstparticle);
            mVectorOfParticlePointersVectors[ii] = ParticlePointerVector( mMaxNumberOfParticles*2 );
            ParticlePointerVector& particle_pointers = mVectorOfParticlePointersVectors[ii];
            //int & number_of_particles = ielem->GetValue(NUMBER_OF_BED_PARTICLES);
            int & number_of_particles = mNumOfParticlesInElems[ii];
            number_of_particles=0;

            Geometry< Node<3> >& geom = ielem->GetGeometry();
            //unsigned int elem_id = ielem->Id();
            ComputeGaussPointPositions_initial(geom, pos, N); //we also have the standard (4), and 45
            //now we seed the particles in the current element
            for (unsigned int j = 0; j < pos.size1(); j++)
            {
                ++particle_id;

                ShallowParticle& pparticle = mParticlesVector[particle_id-1];
                //~ pparticle.X()=pos(j,0);
                //~ pparticle.Y()=pos(j,1);
                //~ pparticle.Z()=pos(j,2);
                pparticle.Coordinates() = row(pos,j);

                pparticle.GetEraseFlag()=false;

                array_1d<float, 3 > & vector1 = pparticle.GetVector1();
                float & scalar1 = pparticle.GetScalar1();
                noalias(vector1) = ZeroVector(3);
                scalar1=0.0;

                for (unsigned int k = 0; k < (TDim+1); k++)
                {
                    scalar1          += N(j, k) * geom[k].FastGetSolutionStepValue(mScalarVar1);
                    noalias(vector1) += N(j, k) * geom[k].FastGetSolutionStepValue(mVectorVar1);
                }

                particle_pointers(j) = &pparticle;
                number_of_particles++ ;

            }
            ++i_int;
        }

        mNParticles=particle_id; //we save the last particle created as the total number of particles we have. For the moment this is true.
        std::cout << "  [Creating particles : " << mNParticles << " particles created]" << std::endl;

        mParticlePrintingToolInitialized=false;

        KRATOS_CATCH("")
    }


    ~MoveShallowWaterParticleUtility()
    {}


    void MountBin()
    {
        KRATOS_TRY

        //copy the elements to a new container, as the list will
        //be shuffled duringthe construction of the tree
        ContainerType& rElements           =  mrModelPart.ElementsArray();
        IteratorType it_begin              =  rElements.begin();
        IteratorType it_end                =  rElements.end();
        //const int number_of_elem         =   rElements.size();

        typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure>(it_begin, it_end  ) );
        paux.swap(mpBinsObjectDynamic);
        //BinsObjectDynamic<Configure>  mpBinsObjectDynamic(it_begin, it_end );

        std::cout << "  finished mounting Bins" << std::endl;

        KRATOS_CATCH("")
    }


    /// Calculates the mean velocity
    /** This function computes the mean velocity within an element and
     * stores it in MEAN_VEL_OVER_ELEM_SIZE variable.
     * This variable keeps the courant number aprox 0.1 in each substep
     *
     * @see MoveParticle
     * @see MoveParticleInverseWay
     */
    void CalculateVelOverElemSize()
    {
        KRATOS_TRY

        //ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();

        const double nodal_weight = 1.0/ (1.0 + double (TDim) );

        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();
        vector<unsigned int> element_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Elements().size(), element_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
            {
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
                Geometry<Node<3> >& geom = ielem->GetGeometry();

                array_1d<double, 3 >vector_mean_velocity=ZeroVector(3);

                for (unsigned int i=0; i != (TDim+1) ; i++)
                    vector_mean_velocity += geom[i].FastGetSolutionStepValue(VELOCITY);
                vector_mean_velocity *= nodal_weight;

                //~ const double mean_velocity = sqrt ( pow(vector_mean_velocity[0],2) + pow(vector_mean_velocity[1],2) + pow(vector_mean_velocity[2],2) );
                const double mean_velocity = norm_2( vector_mean_velocity );
                ielem->SetValue(MEAN_VEL_OVER_ELEM_SIZE, mean_velocity / ( ielem->GetValue(MEAN_SIZE) ) );
            }
        }
        KRATOS_CATCH("")
    }


    /// Reset the boundary conditions
    /** When a variable is fixed this function resets the nodal values
     * with the previous time step
     */
    void ResetBoundaryConditions()
    {
        KRATOS_TRY

        typedef VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        component_type vector_var_x = KratosComponents< component_type >::Get(m_vector_var1_name+std::string("_X"));
        component_type vector_var_y = KratosComponents< component_type >::Get(m_vector_var1_name+std::string("_Y"));
        component_type vector_var_z = KratosComponents< component_type >::Get(m_vector_var1_name+std::string("_Z"));

        ModelPart::NodesContainerType::iterator inodebegin = mrModelPart.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Nodes().size(), node_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;

                if (inode->IsFixed(mScalarVar1))
                {
                    inode->FastGetSolutionStepValue(mScalarVar1)=inode->GetSolutionStepValue(mScalarVar1,1);
                }
                if (inode->IsFixed(vector_var_x))
                {
                    inode->FastGetSolutionStepValue(vector_var_x)=inode->GetSolutionStepValue(vector_var_x,1);
                }
                if (inode->IsFixed(vector_var_y))
                {
                    inode->FastGetSolutionStepValue(vector_var_y)=inode->GetSolutionStepValue(vector_var_y,1);
                }
                if (inode->IsFixed(vector_var_z))
                {
                    inode->FastGetSolutionStepValue(vector_var_z)=inode->GetSolutionStepValue(vector_var_z,1);
                }
            }
        }

        KRATOS_CATCH("")
    }


    /// Auxiliar function to compute the "delta variables"
    /** Delta variables are the difference between two time steps.
     * It's value is used to update particles info
     *
     * @see CorrectParticlesWithoutMovingUsingDeltaVariables
     */
    void CalculateDeltaVariables()
    {
        KRATOS_TRY
        ModelPart::NodesContainerType::iterator inodebegin = mrModelPart.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Nodes().size(), node_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                inode->FastGetSolutionStepValue(DELTA_SCALAR1) = inode->FastGetSolutionStepValue(mScalarVar1) - inode->FastGetSolutionStepValue(PROJECTED_SCALAR1);  //PROJECTED_SCALAR1
                inode->FastGetSolutionStepValue(DELTA_VECTOR1) = inode->FastGetSolutionStepValue(mVectorVar1) - inode->FastGetSolutionStepValue(PROJECTED_VECTOR1);  //PROJECTED_VECTOR1
            }
        }
        KRATOS_CATCH("")
    }


    /// Auxiliar function
    /** This function copy a scalar variable value to the previous time step
     */
    void CopyScalarVarToPreviousTimeStep(const Variable<double>& OriginVariable,
                   ModelPart::NodesContainerType& rNodes)
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
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                inode->GetSolutionStepValue(OriginVariable,1) = inode->FastGetSolutionStepValue(OriginVariable);
            }
        }
        KRATOS_CATCH("")
    }


    /// Auxiliar function
    /** This function copy a vector variable value to the previous time step
     */
    void CopyVectorVarToPreviousTimeStep(const Variable<array_1d<double,3>>& OriginVariable,
                   ModelPart::NodesContainerType& rNodes)
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
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                noalias(inode->GetSolutionStepValue(OriginVariable,1)) = inode->FastGetSolutionStepValue(OriginVariable);
            }
        }
        KRATOS_CATCH("")
    }


    /// Move all the particles
    /** This function moves the particles across the streamlines
     * according to the velocity given by VELOCITY variable. The
     * movement is performed in nsubsteps, during a total time
     * of DELTA_TIME
     *
     * @see Moveparticle
     */
    void MoveParticles()
    {
        KRATOS_TRY

        ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();

        const int offset = mOffset; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
                                    //moveparticlesdiff reads from the pointers of one part (ie odd) and saves into the other part (ie even part)
                                    //since it is the only function in the whole procedure that does this, it must use alternatively one part and the other.

        bool even_timestep;
        if (offset!=0) even_timestep=false;
        else even_timestep=true;

        const int post_offset = mMaxNumberOfParticles * static_cast<int>(even_timestep); //and we also save the offset to know the location in which we will save the pointers after we've moved the particles

        double delta_t = CurrentProcessInfo[DELTA_TIME];

        array_1d<double,TDim+1> N;
        const unsigned int max_results = 10000;

        //double integration_distance= 2.0;

        mMaxSubSteps = 10;
        mMaxSubStepDt = delta_t / static_cast<double>(mMaxSubSteps);

        vector<unsigned int> element_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Elements().size(), element_partition);

        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();


        //before doing anything we must reset the vector of nodes contained by each element (particles that are inside each element.
        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
            {
                //ModelPart::ElementsContainerType::iterator old_element = ielembegin+ii;

                int & number_of_particles = mNumOfParticlesInElems[ii]; //old_element->GetValue(NUMBER_OF_BED_PARTICLES);

                mNumOfParticlesInElemsAux[ii] = number_of_particles;
                mNumOfParticlesInElems[ii] = 0;
                //we reset the local vectors for a faster access;
            }
        }
        std::cout << "convecting particles" << std::endl;
        //We move the particles across the fixed mesh and saving change data into them (using the function MoveParticle)

        #pragma omp barrier

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            ResultContainerType results(max_results);

            WeakPointerVector< Element > elements_in_trajectory;
            elements_in_trajectory.resize(20);

            for(unsigned int ielem = element_partition[kkk]; ielem<element_partition[kkk+1]; ielem++)
            {
                ModelPart::ElementsContainerType::iterator old_element = ielembegin+ielem;
                const int old_element_id = old_element->Id();

                ParticlePointerVector& old_element_particle_pointers = mVectorOfParticlePointersVectors(old_element_id-1);

                if ( (results.size()) != max_results )
                    results.resize(max_results);

                unsigned int number_of_elements_in_trajectory = 0; //excluding the origin one (current one, ielem)

                for(int ii=0; ii<(mNumOfParticlesInElemsAux(ielem)); ii++)
                {
                    ShallowParticle& pparticle = old_element_particle_pointers[offset+ii];

                    Element::Pointer pcurrent_element( *old_element.base() );
                    ResultIteratorType result_begin = results.begin();
                    bool & erase_flag=pparticle.GetEraseFlag();
                    if (erase_flag == false){
                        MoveParticle(pparticle,pcurrent_element,elements_in_trajectory,number_of_elements_in_trajectory,result_begin,max_results); //saqué N de los argumentos, no lo necesito ya q empieza SIEMPRE en un nodo y no me importa donde termina

                        const int current_element_id = pcurrent_element->Id();

                        int & number_of_particles_in_current_elem = mNumOfParticlesInElems(current_element_id-1);

                        if (number_of_particles_in_current_elem < mMaxNumberOfParticles && erase_flag == false)
                        {
                            ParticlePointerVector& current_element_particle_pointers = mVectorOfParticlePointersVectors(current_element_id-1);

                            #pragma omp critical
                            {
                                if (number_of_particles_in_current_elem < mMaxNumberOfParticles) // we cant go over this node, there's no room. otherwise we would be in the position of the first particle of the next element!!
                                {
                                    current_element_particle_pointers(post_offset+number_of_particles_in_current_elem) = &pparticle;
                                    number_of_particles_in_current_elem++ ;
                                    KRATOS_ERROR_IF( number_of_particles_in_current_elem > mMaxNumberOfParticles ) <<
                                        "In move shallow water particle utility: exceeded maximum number of particles" << std::endl;
                                    //~ if (number_of_particles_in_current_elem > mMaxNumberOfParticles)
                                        //~ KRATOS_WATCH("MAL");
                                }
                                else
                                {
                                    pparticle.GetEraseFlag()=true; //so we just delete it!
                                }
                            }
                        }
                        else
                        {
                            pparticle.GetEraseFlag()=true; //so we just delete it!
                        }
                    }
                }
            }
        }

        // After having changed everything we change the status of the mOddTimeStep flag:
        mOffset = post_offset;; //

        KRATOS_CATCH("")
    }


    /// Transfer particles information to the mesh nodes
    /** This function explicitly projects data from particles (lagrangian)
     * onto the eulerian mesh. Shape functions of the elements determine
     * the particle location within the element and its contribution to
     * each node as a weighting function.
     */
    void TransferLagrangianToEulerian() //explicit
    {
        KRATOS_TRY

        //ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        //const double delta_t =CurrentProcessInfo[DELTA_TIME];
        const double threshold = 1e-10 / (static_cast<double>(TDim)+1.0);

        std::cout << "projecting info to mesh" << std::endl;

        const int offset = mOffset;
        // the array of pointers for each element has twice the required size so that
        // we use a part in odd timesteps and the other in even ones.
        //(flag managed only by MoveParticles)

        // We must project data from the particles (lagrangian) onto the eulerian mesh
        //int nnodes = mrModelPart.Nodes().size();
        //array_1d<double,(n_nodes)> eulerian_nodes_sumweights;

        // We save data from previous time step of the eulerian mesh in case we must reuse it later
        // cos no particle was found around the nodes though we could've use a bigger buffer, to be changed later!
        // after having saved data, we reset them to zero, this way it's easier to add the contribution
        // of the surrounding particles.
        ModelPart::NodesContainerType::iterator inodebegin = mrModelPart.NodesBegin();
        vector<unsigned int> node_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Nodes().size(), node_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                inode->FastGetSolutionStepValue(PROJECTED_SCALAR1)=0.0;
                inode->FastGetSolutionStepValue(PROJECTED_VECTOR1)=ZeroVector(3);
                inode->FastGetSolutionStepValue(YP)=0.0;
            }
        }

        // Adding contribution, loop on elements, since each element has stored the particles found inside of it
        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Elements().size(), element_partition);

        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();
        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
            {
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;

                array_1d<double,3*(TDim+1)> nodes_positions;
                array_1d<double,3*(TDim+1)> nodes_added_vector1 = ZeroVector(3*(TDim+1));
                array_1d<double,(TDim+1)> nodes_added_scalar1 = ZeroVector((TDim+1));
                array_1d<double,(TDim+1)> nodes_added_weights = ZeroVector((TDim+1));
                //array_1d<double,(TDim+1)> weighting_inverse_divisor;

                Geometry<Node<3> >& geom = ielem->GetGeometry();

                for (int i=0 ; i!=(TDim+1) ; ++i)
                {
                    nodes_positions[i*3+0]=geom[i].X();
                    nodes_positions[i*3+1]=geom[i].Y();
                    nodes_positions[i*3+2]=geom[i].Z();
                    //weighting_inverse_divisor[i]=1.0/((geom[i].FastGetSolutionStepValue(MEAN_SIZE))*1.01);
                }

                int & number_of_particles_in_elem= mNumOfParticlesInElems[ii];
                ParticlePointerVector&  element_particle_pointers =  mVectorOfParticlePointersVectors[ii];

                for (int iii=0; iii<number_of_particles_in_elem ; iii++ )
                {
                    if (iii==mMaxNumberOfParticles) // It means we are out of our portion of the array, abort loop!
                        break;

                    ShallowParticle& pparticle = element_particle_pointers[offset+iii];

                    if (pparticle.GetEraseFlag()==false)
                    {
                        array_1d<double,3> & position = pparticle.Coordinates();
                        const float& particle_scalar1 = pparticle.GetScalar1();
                        const array_1d<float,3>& particle_vector1 = pparticle.GetVector1();

                        array_1d<double,TDim+1> N;
                        bool is_found = CalculatePosition(nodes_positions,position[0],position[1],position[2],N);
                        if (is_found==false) // Something went wrong. if it was close enough to the edge we simply send it inside the element.
                        {
                            KRATOS_INFO("MoveShallowWaterParticleUtility") << N << std::endl;
                            for (int j=0 ; j!=(TDim+1); j++)
                                if (N[j]<0.0 && N[j]> -1e-5)
                                    N[j]=1e-10;
                        }

                        for (int j=0 ; j!=(TDim+1); j++) //going through the 3/4 nodes of the element
                        {
                            // These lines for a weighting function based on the distance (or square distance) from the node insteadof the shape functions
                            //double sq_dist = 0;
                            //for (int k=0 ; k!=(TDim); k++) sq_dist += ((position[k] - nodes_positions[j*3+k])*(position[k] - nodes_positions[j*3+k]));
                            //double weight = (1.0 - (sqrt(sq_dist)*weighting_inverse_divisor[j] ) );

                            double weight=N(j)*N(j);
                            //weight=N(j)*N(j)*N(j);
                            if (weight<threshold) weight=1e-10;

                            nodes_added_weights[j] += weight;
                            nodes_added_scalar1[j] += weight*static_cast<double>(particle_scalar1);
                            for (int k=0 ; k!=(TDim); k++) //x,y,(z)
                            {
                                nodes_added_vector1[j*3+k] += weight * static_cast<double>(particle_vector1[k]);
                            }
                        }
                    }
                }

                for (int i=0 ; i!=(TDim+1) ; ++i) {
                    geom[i].SetLock();
                    geom[i].FastGetSolutionStepValue(PROJECTED_SCALAR1)   += nodes_added_scalar1[i];

                    geom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_X) += nodes_added_vector1[3*i+0];
                    geom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Y) += nodes_added_vector1[3*i+1];
                    geom[i].FastGetSolutionStepValue(PROJECTED_VECTOR1_Z) += nodes_added_vector1[3*i+2];

                    geom[i].FastGetSolutionStepValue(YP) += nodes_added_weights[i];
                    geom[i].UnSetLock();
                }
            }
        }

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=node_partition[kkk]; ii<node_partition[kkk+1]; ii++)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
                double sum_weights = inode->FastGetSolutionStepValue(YP);
                if (sum_weights>0.00001)
                {
                    double & scalar = inode->FastGetSolutionStepValue(PROJECTED_SCALAR1);
                    array_1d<double,3> & vector = inode->FastGetSolutionStepValue(PROJECTED_VECTOR1);
                    scalar /=sum_weights; // resetting the scalar1
                    vector /=sum_weights; // resetting the vector1
                }
                else // This should never happen because other ways to recover the information have been executed before, but leaving it just in case..
                {
                    inode->FastGetSolutionStepValue(PROJECTED_SCALAR1)=inode->FastGetSolutionStepValue(mScalarVar1,1); // Resetting the convected scalar
                    inode->FastGetSolutionStepValue(PROJECTED_VECTOR1)=inode->FastGetSolutionStepValue(mVectorVar1,1); // Resetting the convected vector
                }
            }
        }

        KRATOS_CATCH("")
    }


    /// Update all the particles without moving them
    /** This function updates all the particles variables using the
     * "delta variables" from the nodal database.
     *
     * @see CorrectParticleUsingDeltaVariables
     */
    void CorrectParticlesWithoutMovingUsingDeltaVariables()
    {
        KRATOS_TRY
        //std::cout << "updating particles" << std::endl;
        //ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();

        const int offset = mOffset; //the array of pointers for each element has twice the required size so that we use a part in odd timesteps and the other in even ones.
                                    //(flag managed only by MoveParticles)
        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();

        vector<unsigned int> element_partition;
        #ifdef _OPENMP
            int number_of_threads = omp_get_max_threads();
        #else
            int number_of_threads = 1;
        #endif
        OpenMPUtils::CreatePartition(number_of_threads, mrModelPart.Elements().size(), element_partition);

        #pragma omp parallel for
        for(int kkk=0; kkk<number_of_threads; kkk++)
        {
            for(unsigned int ii=element_partition[kkk]; ii<element_partition[kkk+1]; ii++)
            {
                //const int & elem_id = ielem->Id();
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
                Element::Pointer pelement(*ielem.base());
                Geometry<Node<3> >& geom = ielem->GetGeometry();

                //ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(BED_PARTICLE_POINTERS));
                //int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_BED_PARTICLES);
                int & number_of_particles_in_elem= mNumOfParticlesInElems[ii];
                ParticlePointerVector&  element_particle_pointers =  mVectorOfParticlePointersVectors[ii];

                for (int iii=0; iii<number_of_particles_in_elem ; iii++ )
                {
                    if (iii>mMaxNumberOfParticles) //it means we are out of our portion of the array, abort loop!
                        break;

                    ShallowParticle & pparticle = element_particle_pointers[offset+iii];

                    bool erase_flag= pparticle.GetEraseFlag();
                    if (erase_flag==false)
                    {
                        CorrectParticleUsingDeltaVariables(pparticle,pelement,geom); //'lite' version, we pass by reference the geometry, so much cheaper
                    }
                }
            }
        }
        KRATOS_CATCH("")
    }


    /// AddUniqueWeakPointer
    template< class TDataType >
    void AddUniqueWeakPointer
        (WeakPointerVector< TDataType >& v, const typename TDataType::WeakPointer candidate)
    {
        typename WeakPointerVector< TDataType >::iterator i = v.begin();
        typename WeakPointerVector< TDataType >::iterator endit = v.end();
        while ( i != endit && (i)->Id() != (candidate.lock())->Id())
        {
            i++;
        }
        if( i == endit )
        {
            v.push_back(candidate);
        }

    }


    /// Fill an element with particles
    /** This function is to be executed after moving particles and
     * before tranferring data from lagrangian particles to eulerian mesh
     * If an element finishes with less particles than "minimum number
     * of particles", then PreReseed adds particles inside it.
     * A minimal reseed is performed in order to not disturb the projection
     * from lagrangian to euelrian.
     *
     * @see MinimumNumberOfParticles
     *
     * @see MoveParticles
     * @see MoveParticleInverseWay: is called to get the particle values
     */
    void PreReseed(int MinimumNumberOfParticles)
    {
        KRATOS_TRY

        //ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        const int offset =mOffset;
        const int max_results = 1000;

        //tools for the paralelization
        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        vector<unsigned int> elem_partition;
        int number_of_rows = mrModelPart.Elements().size();
        elem_partition.resize(number_of_threads + 1);
        int elem_partition_size = number_of_rows / number_of_threads;
        elem_partition[0] = 0;
        elem_partition[number_of_threads] = number_of_rows;
        //KRATOS_WATCH(elem_partition_size);
        for (unsigned int i = 1; i < number_of_threads; i++)
        elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();

        #pragma omp parallel firstprivate(elem_partition)
        {
            ResultContainerType results(max_results);
            int k = OpenMPUtils::ThisThread();
            //ModelPart::ElementsContainerType::iterator it_begin = mrModelPart.ElementsBegin() +  elem_partition[k];
            //ModelPart::ElementsContainerType::iterator it_end = mrModelPart.ElementsBegin() + elem_partition[k+1] ;
            //ModelPart::NodesContainerType local_list=aux[k];
            //PointerVectorSet<ShallowParticle, IndexedObject> & list=aux[k];
            BoundedMatrix<double, (TDim+1), 3 > pos;
            BoundedMatrix<double, (TDim+1) , (TDim+1) > N;
            unsigned int freeparticle=0; //we start with the first position in the particles array

            //int local_id=1;
            for(unsigned int ii=elem_partition[k]; ii<elem_partition[k+1]; ii++)
            {
                //const int & elem_id = ielem->Id();
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;
                results.resize(max_results);
                //const int & elem_id = ielem->Id();
                //ParticlePointerVector&  element_particle_pointers =  (ielem->GetValue(BED_PARTICLE_POINTERS));
                //int & number_of_particles_in_elem=ielem->GetValue(NUMBER_OF_BED_PARTICLES);
                int & number_of_particles_in_elem = mNumOfParticlesInElems[ii];
                ParticlePointerVector&  element_particle_pointers =  mVectorOfParticlePointersVectors[ii];
                if (number_of_particles_in_elem < (MinimumNumberOfParticles)) // && (ielem->GetGeometry())[0].Y()<0.10 )
                {
                    Geometry< Node<3> >& geom = ielem->GetGeometry();
                    ComputeGaussPointPositionsForPreReseed(geom, pos, N);

                    for (unsigned int j = 0; j < (pos.size1()); j++) // I am dropping the last one, the one in the middle of the element
                    {
                        bool keep_looking = true;
                        while(keep_looking)
                        {
                            if (mParticlesVector[freeparticle].GetEraseFlag()==true)
                            {
                                #pragma omp critical
                                {
                                    if (mParticlesVector[freeparticle].GetEraseFlag()==true)
                                    {
                                        mParticlesVector[freeparticle].GetEraseFlag()=false;
                                        keep_looking=false;
                                    }
                                }
                                if (keep_looking==false)
                                    break;
                                else
                                    freeparticle++;
                            }
                            else
                                freeparticle++;
                        }

                        ShallowParticle pparticle(pos(j,0),pos(j,1),pos(j,2));

                        array_1d<double,TDim+1>aux2_N;
                        bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux2_N);
                        KRATOS_ERROR_IF_NOT( is_found ) <<
                            "In move shallow water particle utility: particle not found in domain" << std::endl;

                        pparticle.GetEraseFlag()=false;

                        ResultIteratorType result_begin = results.begin();
                        Element::Pointer pelement( *ielem.base() );
                        MoveParticleInverseWay(pparticle, pelement, result_begin, max_results);

                        //and we copy it to the array:
                        mParticlesVector[freeparticle] =  pparticle;

                        element_particle_pointers(offset+number_of_particles_in_elem) = &mParticlesVector[freeparticle];
                        pparticle.GetEraseFlag()=false;

                        number_of_particles_in_elem++;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    /// Fill an element with particles
    /** This function is to be executed after the mesh stage solver is
     * called and the particles are updated.
     * If an element contains less particles than "minimum number of
     * particles", then PostReseed adds particles inside it.
     * A full reseed is performed and the particle gets it's convected
     * variables directly from the eulerian mesh
     *
     * @param MinimumNumberOfParticles
     *
     * @see PreReseed
     */
    void PostReseed(int MinimumNumberOfParticles) //pooyan's way
    {
        KRATOS_TRY

        //ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        const int offset = mOffset;

        //TOOLS FOR THE PARALELIZATION
        //int last_id= (mr_linea_model_part.NodesEnd()-1)->Id();
        unsigned int number_of_threads = OpenMPUtils::GetNumThreads();
        vector<unsigned int> elem_partition;
        int number_of_rows=mrModelPart.Elements().size();
        //KRATOS_THROW_ERROR(std::logic_error, "Add  ----NODAL_H---- variable!!!!!! ERROR", "");
        elem_partition.resize(number_of_threads + 1);
        int elem_partition_size = number_of_rows / number_of_threads;
        elem_partition[0] = 0;
        elem_partition[number_of_threads] = number_of_rows;

        for (unsigned int i = 1; i < number_of_threads; i++)
        elem_partition[i] = elem_partition[i - 1] + elem_partition_size;
        //typedef Node < 3 > PointType;
        //std::vector<ModelPart::NodesContainerType> aux;// aux;
        //aux.resize(number_of_threads);

        //ModelPart::NodesContainerType::iterator it_begin_particle_model_part = mr_linea_model_part.NodesBegin();
        //ModelPart::NodesContainerType::iterator it_end_particle_model_part = mr_linea_model_part.NodesEnd();
        ModelPart::ElementsContainerType::iterator ielembegin = mrModelPart.ElementsBegin();

        #pragma omp parallel firstprivate(elem_partition) // firstprivate(results)//we will add the nodes in different parts of aux and later assemple everything toghether, remaming particles ids to get consecutive ids
        {
            unsigned int reused_particles=0;

            unsigned int freeparticle = 0; //we start by the first position;

            int k = OpenMPUtils::ThisThread();

            BoundedMatrix<double, (3+2*TDim), 3 > pos; //7 particles (2D) or 9 particles (3D)
            BoundedMatrix<double, (3+2*TDim), (TDim+1) > N;

            double mesh_scalar1;
            array_1d<double,3> mesh_vector1;

            array_1d<int, (3+2*TDim) > positions;

            unsigned int number_of_reseeded_particles;

            for(unsigned int ii=elem_partition[k]; ii<elem_partition[k+1]; ii++)
            {
                //const int & elem_id = ielem->Id();
                ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii;

                int & number_of_particles_in_elem = mNumOfParticlesInElems[ii];
                ParticlePointerVector&  element_particle_pointers =  mVectorOfParticlePointersVectors[ii];

                Geometry< Node<3> >& geom = ielem->GetGeometry();
                if ( number_of_particles_in_elem < (MinimumNumberOfParticles) ) // && (geom[0].Y()<0.10) ) || (number_of_water_particles_in_elem>2 && number_of_particles_in_elem<(MinimumNumberOfParticles) ) )
                {
                    //bool reseed_more=false;
                    number_of_reseeded_particles = 0;

                    //reseed_more=true;
                    number_of_reseeded_particles = 3 + 2*TDim;
                    ComputeGaussPointPositionsForPostReseed(geom, pos, N);

                    for (unsigned int j = 0; j < number_of_reseeded_particles; j++)
                    {
                        // Now we have to find an empty space (a particle that was about to be deleted) in the
                        // particles model part. once found. there will be our renewed particle:
                        bool keep_looking = true;
                        while(keep_looking)
                        {
                            if (mParticlesVector[freeparticle].GetEraseFlag()==true)
                            {
                                #pragma omp critical
                                {
                                    if (mParticlesVector[freeparticle].GetEraseFlag()==true)
                                    {
                                        mParticlesVector[freeparticle].GetEraseFlag()=false;
                                        keep_looking=false;
                                    }
                                }
                                if (keep_looking==false)
                                    break;

                                else
                                    freeparticle++;
                            }
                            else
                                freeparticle++;
                        }

                        ShallowParticle pparticle(pos(j,0),pos(j,1),pos(j,2));

                        array_1d<double,TDim+1>aux_N;
                        bool is_found = CalculatePosition(geom,pos(j,0),pos(j,1),pos(j,2),aux_N);
                        KRATOS_ERROR_IF_NOT( is_found ) <<
                            "In move shallow water particle utility: particle not found in domain" << std::endl;

                        mesh_scalar1 = 0.0;
                        mesh_vector1 = ZeroVector(3);

                        for (unsigned int l = 0; l < (TDim+1); l++)
                        {
                            mesh_scalar1 +=  N(j,l) * geom[l].FastGetSolutionStepValue(mScalarVar1);
                            noalias(mesh_vector1) += N(j, l) * geom[l].FastGetSolutionStepValue(mVectorVar1);
                        }
                        pparticle.GetScalar1()=mesh_scalar1;
                        pparticle.GetVector1()=mesh_vector1;
                        pparticle.GetEraseFlag()=false;

                        mParticlesVector[freeparticle]=pparticle;
                        element_particle_pointers(offset+number_of_particles_in_elem) = &mParticlesVector[freeparticle];
                        number_of_particles_in_elem++;

                        KRATOS_ERROR_IF( keep_looking ) <<
                            "In move shallow water particle utility: Finished the list and couldnt find a free cell for the new particle!" << std::endl;
                        reused_particles++;
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }


    /// Fill a model part with particles
    /** This function prints the particles to a model part
     *
     * @param rLagrangianModelPart: empty model part to print particles
     * @param FilterFactor: the function will print one particle of every "filter factor"
     */
    void ExecuteParticlesPrintingTool( ModelPart& rLagrangianModelPart, unsigned int FilterFactor )
    {
        KRATOS_TRY
        // We will only print one out of every "filter factor" particles of the total particle list

        if (mParticlePrintingToolInitialized == false)
        {
            KRATOS_ERROR_IF( rLagrangianModelPart.NodesBegin() - rLagrangianModelPart.NodesEnd() > 0 ) <<
                "In move shallow water particle utility: an empty model part is required for the particles printing tool" << std::endl;

            rLagrangianModelPart.AddNodalSolutionStepVariable(mScalarVar1);
            rLagrangianModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);

            for (unsigned int i = 0; i != ((mMaxNumberOfParticles*mNElems)/FilterFactor) + FilterFactor; i++)
            {
                Node < 3 > ::Pointer pnode = rLagrangianModelPart.CreateNewNode( i+mLastNodeId+1 , 0.0, 0.0, 0.0);  //recordar que es el nueevo model part!!
                //pnode->SetBufferSize(mrModelPart.NodesBegin()->GetBufferSize());
                pnode->SetBufferSize(1);
            }
            mParticlePrintingToolInitialized=true;
        }

        // Resetting data of the unused particles
        const double inactive_particle_position = -10.0;
        array_1d<double,3>inactive_particle_position_vector;
        inactive_particle_position_vector(0)=inactive_particle_position;
        inactive_particle_position_vector(1)=inactive_particle_position;
        inactive_particle_position_vector(2)=inactive_particle_position;
        ModelPart::NodesContainerType::iterator inodebegin = rLagrangianModelPart.NodesBegin();
        for(unsigned int ii = 0; ii < rLagrangianModelPart.Nodes().size(); ii++)
        {
            ModelPart::NodesContainerType::iterator inode = inodebegin+ii;
            inode->FastGetSolutionStepValue(mScalarVar1)  = 0.0;
            inode->FastGetSolutionStepValue(DISPLACEMENT) = inactive_particle_position_vector;
        }

        int counter = 0;
        //ModelPart::NodesContainerType::iterator it_begin = rLagrangianModelPart.NodesBegin();
        for (int i = 0; i != mMaxNumberOfParticles*mNElems; i++)
        {
            ShallowParticle& pparticle = mParticlesVector[i];
            if(pparticle.GetEraseFlag() == false && i%FilterFactor == 0)
            {
                ModelPart::NodesContainerType::iterator inode = inodebegin + counter; //copying info from the particle to the (printing) node.
                inode->FastGetSolutionStepValue(mScalarVar1)  = pparticle.GetScalar1();
                inode->FastGetSolutionStepValue(DISPLACEMENT) = pparticle.Coordinates();
                counter++;
            }
        }

        KRATOS_CATCH("")
    }


protected:


private:


    /// Move a particle
    /** this function moves a particle according to the velocity given
     * by VELOCITY variable. The movement is performed in nsubsteps,
     * during a total time of DELTA_TIME
     *
     * @param pParticle
     * @param pElement
     * @param rElementsInTrajectory
     * @param rNumberOfElementsInTrajectory
     * @param ResultBegin
     * @param MaxNumberOfResults
     *
     * @see MoveParticles
     */
    void MoveParticle(ShallowParticle & pParticle,
                      Element::Pointer & pElement,
                      WeakPointerVector< Element >& rElementsInTrajectory,
                      unsigned int & rNumberOfElementsInTrajectory,
                      ResultIteratorType ResultBegin,
                      const unsigned int MaxNumberOfResults)
    {
        ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        unsigned int nsubsteps;
        double substep_dt;

        bool keep_integrating = false;
        bool is_found;

        array_1d<double,3> vel;
        array_1d<double,3> vel_without_other_phase_nodes=ZeroVector(3);
        array_1d<double,3> position;
        array_1d<double,3> mid_position;
        array_1d<double,TDim+1> N;

        //we start with the first position, then it will enter the loop.
        position = pParticle.Coordinates(); //initial coordinates

        double only_integral  = 0.0 ;

        is_found = FindNodeOnMesh(position, N, pElement, ResultBegin, MaxNumberOfResults); //good, now we know where this point is:
        if(is_found == true)
        {
            keep_integrating=true;
            Geometry< Node<3> >& geom = pElement->GetGeometry();//the element we're in
            vel=ZeroVector(3);

            for(unsigned int j=0; j<(TDim+1); j++)
            {
                noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j];
            }

            //calculating substep to get +- courant(substep) = 0.1
            nsubsteps = 10.0 * (delta_t * pElement->GetValue(MEAN_VEL_OVER_ELEM_SIZE));
            if (nsubsteps<1)
                nsubsteps=1;
            substep_dt = delta_t / double(nsubsteps);

            only_integral = 1.0;// weight;//*double(nsubsteps);
            position += vel*substep_dt;//weight;

            // DONE THE FIRST LOCATION OF THE PARTICLE, NOW WE PROCEED TO STREAMLINE INTEGRATION USING THE MESH VELOCITY
            unsigned int check_from_element_number = 0;

            for(unsigned int i=0; i<(nsubsteps-1); i++)// this is for the substeps n+1. in the first one we already knew the position of the particle.
            {
                if (keep_integrating == true)
                {
                    is_found = FindNodeOnMesh(position, N, pElement, rElementsInTrajectory, rNumberOfElementsInTrajectory, check_from_element_number, ResultBegin, MaxNumberOfResults); //good, now we know where this point is:
                    if(is_found == true)
                    {
                        Geometry< Node<3> >& geom = pElement->GetGeometry();//the element we're in

                        vel = ZeroVector(3);
                        for(unsigned int j=0; j<(TDim+1); j++)
                        {
                            noalias(vel) += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j];
                        }

                        only_integral += 1.0; //values saved for the current time step

                        position+=vel*substep_dt;//weight;

                    }
                    else
                    {
                        keep_integrating=false;
                        break;
                    }
                }
                else
                    break;
            }
        }

        if (keep_integrating == false) (pParticle.GetEraseFlag()=true);
        else is_found = FindNodeOnMesh(position, N ,pElement,ResultBegin,MaxNumberOfResults); //we must save the pointer of the last element that we're in (inside the pointervector pElement)

        if (is_found == false) ( pParticle.GetEraseFlag()=true);

         pParticle.Coordinates() = position;
    }


    /// This function updates a particle
    /** This function updates a particle variables using the "delta
     * variables" from the nodal database.
     *
     * @param pParticle
     * @param pElement
     * @param rGeom
     *
     * @see CorrectParticlesWithoutMovingUsingDeltaVariables
     */
    void CorrectParticleUsingDeltaVariables(ShallowParticle & pParticle,
                                            Element::Pointer & pElement,
                                            Geometry< Node<3> >& rGeom)
    {
        array_1d<double,TDim+1> N;

        //we start with the first position, then it will enter the loop.
        array_1d<double,3> coords = pParticle.Coordinates();
        float & particle_scalar1 = pParticle.GetScalar1();
        array_1d<float,3> & particle_vector1 = pParticle.GetVector1();

        //double distance=0.0;
        double delta_scalar1 = 0.0;
        array_1d<double,3> delta_vector1 = ZeroVector(3);

        bool is_found = CalculatePosition(rGeom,coords[0],coords[1],coords[2],N);
        if(is_found == false)
        {
            KRATOS_INFO("MoveShallowWaterParticleUtility") << N << std::endl;
            for (int j=0 ; j!=(TDim+1); j++)
                                if (N[j]<0.0 )
                                    N[j]=1e-10;
        }

        for(unsigned int j=0; j<(TDim+1); j++)
        {
            delta_scalar1          += rGeom[j].FastGetSolutionStepValue(DELTA_SCALAR1)*N[j];
            noalias(delta_vector1) += rGeom[j].FastGetSolutionStepValue(DELTA_VECTOR1)*N[j];
        }
        particle_scalar1 = particle_scalar1 + delta_scalar1;
        particle_vector1 = particle_vector1 + delta_vector1;
    }


    /// Move a particle in the inverse way
    /** this function moves a particle according to the -velocity given
     * by VELOCITY variable. The movement is performed by a backward
     * integration in nsubsteps, during a total time of DELTA_TIME
     * Before the particle goes out of the element, gets the value
     * of the eulerian mesh and stores it
     *
     * @param pParticle
     * @param pElement
     * @param ResultBegin
     * @param MaxNumberOfResults
     *
     * @see PreReseed
     */
    void MoveParticleInverseWay(ShallowParticle & pParticle,
                                  Element::Pointer & pElement, //NOT A REFERENCE!! WE SHALL NOT OVERWRITE THE ELEMENT IT BELONGS TO!
                                  ResultIteratorType ResultBegin,
                                  const unsigned int MaxNumberOfResults)
    {
        ProcessInfo& CurrentProcessInfo = mrModelPart.GetProcessInfo();
        double delta_t = CurrentProcessInfo[DELTA_TIME];
        unsigned int nsubsteps;
        double substep_dt;

        bool keep_integrating = false;
        bool is_found;

        double scalar1 = 0.0;
        array_1d<double,3> vector1;
        array_1d<double,3> vel;
        array_1d<double,3> position;
        array_1d<double,3> mid_position;
        array_1d<double,TDim+1> N;

        //we start with the first position, then it will enter the loop.
        position = pParticle.Coordinates(); // + (pParticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates

        double only_integral  = 0.0 ;

        is_found = FindNodeOnMesh(position, N, pElement, ResultBegin, MaxNumberOfResults); //good, now we know where this point is:
        if(is_found == true)
        {
            keep_integrating = true;
            Geometry< Node<3> >& geom = pElement->GetGeometry(); //the element we're in

            scalar1 = 0.0;
            vector1 = ZeroVector(3);
            vel     = ZeroVector(3);

            for(unsigned int j=0; j<(TDim+1); j++)
            {
                scalar1          += geom[j].FastGetSolutionStepValue(mScalarVar1)*N[j];
                noalias(vector1) += geom[j].FastGetSolutionStepValue(mVectorVar1)*N[j];
                noalias(vel)     += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j];
            }
            //calculating substep to get +- courant(substep) = 1/4
            nsubsteps = 10.0 * (delta_t * pElement->GetValue(MEAN_VEL_OVER_ELEM_SIZE));
            if (nsubsteps<1)
                nsubsteps=1;
            substep_dt = delta_t / double(nsubsteps);

            only_integral = 1.0; // weight;//*double(nsubsteps);
            position -= vel*substep_dt; //weight;

            for(unsigned int i=0; i<(nsubsteps-1); i++) // this is for the substeps n+1. in the first one we already knew the position of the particle.
            {
                if (keep_integrating == true)
                {
                    is_found = FindNodeOnMesh(position, N, pElement, ResultBegin, MaxNumberOfResults); //good, now we know where this point is:
                    if (is_found == true)
                    {
                        Geometry< Node<3> >& geom = pElement->GetGeometry();//the element we're in

                        scalar1 = 0.0;
                        vector1 = ZeroVector(3);
                        vel     = ZeroVector(3);

                        for(unsigned int j=0; j<(TDim+1); j++)
                        {
                            scalar1          += geom[j].FastGetSolutionStepValue(mScalarVar1)*N(j);
                            noalias(vector1) += geom[j].FastGetSolutionStepValue(mVectorVar1)*N[j];
                            noalias(vel)     += geom[j].FastGetSolutionStepValue(VELOCITY)*N[j];
                        }

                        only_integral += 1.0; //weight ; //values saved for the current time step
                        position -= vel*substep_dt; //weight;
                    }
                    else keep_integrating = false;
                }
            }

            pParticle.GetScalar1() = scalar1;
            pParticle.GetVector1() = vector1;
        }
    }


    /// Find the element into which a given node is located
    /** This function should find the element into which a given node
     * is located and return a pointer to the element and the vector
     * containing the shape functions that define the positions within
     * the element.
     * If false is returned the element is not found
     *
     * @param position of the node
     * @param N: return shape functions that define the positions within the elem
     * @param pElement: return a pointer to the element
     * @param ResultBegin
     * @param MaxNumberOfResults
     * @return FindNodeOnMesh if the element is found of not
     *
     * @see CalculatePosition
     */
    bool FindNodeOnMesh( const array_1d<double,3>& rPosition,
                         array_1d<double,TDim+1>& N,
                         Element::Pointer & pElement,
                         ResultIteratorType ResultBegin,
                         const unsigned int MaxNumberOfResults)
    {
        typedef std::size_t SizeType;

        array_1d<double,TDim+1> aux_N;
        //before using the bin to search for possible elements we check first the last element in which the particle was.
        Geometry<Node<3> >& geom_default = pElement->GetGeometry(); //(*(i))->GetGeometry();
        bool is_found_1 = CalculatePosition(geom_default,rPosition[0],rPosition[1],rPosition[2],N);
        if (is_found_1) //that was easy!
        {
            return true;
        }

        // To begin with we check the neighbour elements; it is a bit more expensive
        WeakPointerVector< Element >& neighb_elems = pElement->GetValue(NEIGHBOUR_ELEMENTS);
        for (unsigned int i=0;i!=(neighb_elems.size());i++)
        {
            Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
            bool is_found_2 = CalculatePosition(geom,rPosition[0],rPosition[1],rPosition[2],N);
            if (is_found_2)
            {
                pElement=Element::Pointer(((neighb_elems(i))));
                return true;
            }
        }

        // If checking all the neighbour elements did not work, we have to use the bins
        // ask to the container for the list of candidate elements
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(rPosition, ResultBegin, MaxNumberOfResults );

        if (results_found>0)
        {
            //loop over the candidate elements and check if the particle falls within
            for(SizeType i = 0; i< results_found; i++)
            {
                Geometry<Node<3> >& geom = (*(ResultBegin+i))->GetGeometry();

                //find local position
                bool is_found_3 = CalculatePosition(geom,rPosition[0],rPosition[1],rPosition[2],N);

                if (is_found_3)
                {
                    pElement=Element::Pointer((*(ResultBegin+i)));
                    return true;
                }
            }
        }

        //if nothing worked, then:
        //not found case
        return false;
    }


    /// Find the element into which a given node is located
    /** This function should find the element into which a given node
     * is located and return a pointer to the element and the vector
     * containing the shape functions that define the positions within
     * the element.
     * If false is returned the element is not found
     * This version includes predefined elements following a trajectory
     *
     * @param rPosition of the node
     * @param N Output shape functions that define the positions within the elem
     * @param pElement Output a pointer to the element
     * @param rElementsInTrajectory
     * @param rNumberOfElementsInTrajectory Output
     * @param CheckFromElementNumber
     * @param ResultBegin
     * @param MaxNumberOfResults
     * @return FindNodeOnMesh if the element is found of not
     *
     * @see CalculatePosition
     */
    bool FindNodeOnMesh( const array_1d<double,3>& rPosition,
                         array_1d<double,TDim+1>& N,
                         Element::Pointer & pElement,
                         WeakPointerVector< Element >& rElementsInTrajectory,
                         unsigned int & rNumberOfElementsInTrajectory,
                         unsigned int & rCheckFromElementNumber,
                         ResultIteratorType ResultBegin,
                         const unsigned int MaxNumberOfResults)
    {
        typedef std::size_t SizeType;

        //~ const array_1d<double,3>& coords = rPosition;
         array_1d<double,TDim+1> aux_N;
        //before using the bin to search for possible elements we check first the last element in which the particle was.
        Geometry<Node<3> >& geom_default = pElement->GetGeometry(); //(*(i))->GetGeometry();
        bool is_found_1 = CalculatePosition(geom_default,rPosition[0],rPosition[1],rPosition[2],N);
        if(is_found_1 == true)
        {
            return true; //that was easy!
        }

        // If it was not found in the first element, we can proceed to check in the following elements (in the trajectory defined by previous particles that started from the same element.
        for (unsigned int i=(rCheckFromElementNumber);i!=rNumberOfElementsInTrajectory;i++)
        {
            Geometry<Node<3> >& geom = rElementsInTrajectory[i].GetGeometry();
            bool is_found_2 = CalculatePosition(geom,rPosition[0],rPosition[1],rPosition[2],aux_N);
            if (is_found_2)
            {
                pElement = Element::Pointer(((rElementsInTrajectory(i))));
                N = aux_N;
                rCheckFromElementNumber = i+1 ; //now i element matches pElement, so to avoid cheching twice the same element we send the counter to the following element.
                return true;
            }
        }

        // Now we check the neighbour elements:
        WeakPointerVector< Element >& neighb_elems = pElement->GetValue(NEIGHBOUR_ELEMENTS);
        for (unsigned int i=0;i!=(neighb_elems.size());i++)
        {
            Geometry<Node<3> >& geom = neighb_elems[i].GetGeometry();
            bool is_found_2 = CalculatePosition(geom,rPosition[0],rPosition[1],rPosition[2],N);
            if (is_found_2)
            {
                pElement=Element::Pointer(((neighb_elems(i))));
                if (rNumberOfElementsInTrajectory<20)
                {
                    rElementsInTrajectory(rNumberOfElementsInTrajectory) = pElement;
                    rNumberOfElementsInTrajectory++;
                    rCheckFromElementNumber = rNumberOfElementsInTrajectory;  //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the rElementsInTrajectory list. we are the particle that is adding elements to the list
                }
                return true;
            }
        }

        // If checking all the neighbour elements did not work, we have to use the bins
        // ask to the container for the list of candidate elements
        SizeType results_found = mpBinsObjectDynamic->SearchObjectsInCell(rPosition, ResultBegin, MaxNumberOfResults );

        if(results_found>0)
        {
            //loop over the candidate elements and check if the particle falls within
            for(SizeType i = 0; i< results_found; i++)
            {
                Geometry<Node<3> >& geom = (*(ResultBegin+i))->GetGeometry();

                //find local position
                bool is_found = CalculatePosition(geom,rPosition[0],rPosition[1],rPosition[2],N);

                if (is_found)
                {
                    pElement=Element::Pointer((*(ResultBegin+i)));
                    if (rNumberOfElementsInTrajectory<20)
                    {
                        rElementsInTrajectory(rNumberOfElementsInTrajectory) = pElement;
                        rNumberOfElementsInTrajectory++;
                        rCheckFromElementNumber = rNumberOfElementsInTrajectory;  //we do it after doing the ++ to the counter, so we woudlnt enter the loop that searches in the rElementsInTrajectory list. we are the particle that is adding elements to the list
                    }
                    return true;
                }
            }
        }

        //not found case
        return false;
    }


    /// Calculate the position of a given particle inside an element
    /** This function calculates the position of a given particle inside
     * an element and returns the shape functions that define it position
     * within the element and returns false if the particle is otuside
     * the element
     *
     * @param rGeom: the element (a triangle)
     * @param xc: the postition of the particle
     * @param yc: the postition of the particle
     * @param zc: the postition of the particle
     * @param N: the shape functions to define the particle position
     *
     * @return CalculatePosition
     */
    inline bool CalculatePosition( const Geometry<Node < 3 > >&rGeom,
                                   const double xc,
                                   const double yc,
                                   const double zc,
                                   array_1d<double,3> & N )
    {
        double x0 = rGeom[0].X();
        double y0 = rGeom[0].Y();
        double x1 = rGeom[1].X();
        double y1 = rGeom[1].Y();
        double x2 = rGeom[2].X();
        double y2 = rGeom[2].Y();

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        KRATOS_ERROR_IF( area == 0.0 ) << "In move shallow water particle utility: element with zero area found" << std::endl;
        double inv_area = 1.0 / area;

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0)
            //if the xc yc is inside the triangle return true
            return true;

        return false;
    }


    /// Calculate the position of a given particle inside an element
    /** This function calculates the position of a given particle inside
     * an element and returns the shape functions that define it position
     * within the element and returns false if the particle is otuside
     * the element
     *
     * @param rNodesPositions of the element (a triangle)
     * @param xc: the postition of the particle
     * @param yc: the postition of the particle
     * @param zc: the postition of the particle
     * @param N: the shape functions to define the particle position
     *
     * @return CalculatePosition
     */
    inline bool CalculatePosition( const array_1d<double,3*(TDim+1)>& rNodesPositions,
                                   const double xc,
                                   const double yc,
                                   const double zc,
                                   array_1d<double,3> & N )
    {
        const double& x0 = rNodesPositions[0];
        const double& y0 = rNodesPositions[1];
        const double& x1 = rNodesPositions[3];
        const double& y1 = rNodesPositions[4];
        const double& x2 = rNodesPositions[6];
        const double& y2 = rNodesPositions[7];

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        KRATOS_ERROR_IF( area == 0.0 ) << "In move shallow water particle utility: element with zero area found" << std::endl;
        double inv_area = 1.0 / area;

        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0)
            //if the xc yc is inside the triangle return true
            return true;

        return false;
    }


    /// Calculate the position of a given particle inside an element
    /** This function calculates the position of a given particle inside
     * an element and returns the shape functions that define it position
     * within the element and returns false if the particle is otuside
     * the element
     *
     * @param rGeom: the element (a tetrahedron)
     * @param xc: the postition of the particle
     * @param yc: the postition of the particle
     * @param zc: the postition of the particle
     * @param N: the shape functions to define the particle position
     *
     * @return CalculatePosition
     */
    inline bool CalculatePosition( const Geometry<Node < 3 > >&rGeom,
                                   const double xc,
                                   const double yc,
                                   const double zc,
                                   array_1d<double, 4 > & N )
    {
        double x0 = rGeom[0].X();
        double y0 = rGeom[0].Y();
        double z0 = rGeom[0].Z();
        double x1 = rGeom[1].X();
        double y1 = rGeom[1].Y();
        double z1 = rGeom[1].Z();
        double x2 = rGeom[2].X();
        double y2 = rGeom[2].Y();
        double z2 = rGeom[2].Z();
        double x3 = rGeom[3].X();
        double y3 = rGeom[3].Y();
        double z3 = rGeom[3].Z();

        double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
        KRATOS_ERROR_IF( vol == 0.0 ) << "In move shallow water particle utility: element with zero vol found" << std::endl;
        double inv_vol = 1.0 / vol;

        N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
            //if the xc yc zc is inside the tetrahedron return true
            return true;

        return false;
    }


    /// Calculate the position of a given particle inside an element
    /** This function calculates the position of a given particle inside
     * an element and returns the shape functions that define it position
     * within the element and returns false if the particle is otuside
     * the element
     *
     * @param rNodesPositions of the element (a tetrahedron)
     * @param xc: the postition of the particle
     * @param yc: the postition of the particle
     * @param zc: the postition of the particle
     * @param N: the shape functions to define the particle position
     *
     * @return CalculatePosition
     */
    inline bool CalculatePosition( const array_1d<double,3*(TDim+1)>& rNodesPositions,
                                   const double xc,
                                   const double yc,
                                   const double zc,
                                   array_1d<double, 4 > & N )
    {
        const double& x0 = rNodesPositions[0];
        const double& y0 = rNodesPositions[1];
        const double& z0 = rNodesPositions[2];
        const double& x1 = rNodesPositions[3];
        const double& y1 = rNodesPositions[4];
        const double& z1 = rNodesPositions[5];
        const double& x2 = rNodesPositions[6];
        const double& y2 = rNodesPositions[7];
        const double& z2 = rNodesPositions[8];
        const double& x3 = rNodesPositions[9];
        const double& y3 = rNodesPositions[10];
        const double& z3 = rNodesPositions[11];

        double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);
        KRATOS_ERROR_IF( vol == 0.0 ) << "In move shallow water particle utility: element with zero vol found" << std::endl;
        double inv_vol = 1.0 / vol;

        N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
        N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
        N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;

        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
            //if the xc yc zc is inside the tetrahedron return true
            return true;

        return false;
    }


    /// Calculate the volume
    /** This function computes the area of a triangle
     */
    inline double CalculateVol( const double x0, const double y0,
                                const double x1, const double y1,
                                const double x2, const double y2 )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }


    /// Calculate the volume
    /** This function computes the volume of a tetrahedron
     */
    inline double CalculateVol( const double x0, const double y0, const double z0,
                                const double x1, const double y1, const double z1,
                                const double x2, const double y2, const double z2,
                                const double x3, const double y3, const double z3 )
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


    /// Compute the Gauss points
    /**
     */
    void ComputeGaussPointPositions_4( Geometry< Node < 3 > >& geom,
                                       BoundedMatrix<double, 7, 3 > & pos,
                                       BoundedMatrix<double, 7, 3 > & N )
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


    /// Compute the Gauss points
    /** For a triangle
     *
     * @see PostReseed
     */
    void ComputeGaussPointPositionsForPostReseed( Geometry< Node < 3 > >& geom,
                                                  BoundedMatrix<double, 7, 3 > & pos,
                                                  BoundedMatrix<double, 7, 3 > & N ) //2d
    {
        double one_third = 1.0 / 3.0;
        double one_eight = 0.12; //1.0 / 6.0;
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


    /// Compute the Gauss points
    /** For a tetrahedron
     *
     * @see PostReseed
     */
    void ComputeGaussPointPositionsForPostReseed( Geometry< Node < 3 > >& geom,
                                                  BoundedMatrix<double, 9, 3 > & pos,
                                                  BoundedMatrix<double, 9, 4 > & N ) //3D
    {
        double one_quarter = 0.25;
        double small_fraction = 0.1; //1.0 / 6.0;
        double big_fraction = 0.7; //2.0 * one_third;
        double mid_fraction = 0.3; //2.0 * one_third;

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

        pos=ZeroMatrix(9,3);
        for (unsigned int i=0; i!=4; i++) //going through the 4 nodes
        {
            array_1d<double, 3 > & coordinates = geom[i].Coordinates();
            for (unsigned int j=0; j!=9; j++) //going through the 9 particles
            {
                for (unsigned int k=0; k!=3; k++) //x,y,z
                    pos(j,k) += N(j,i) * coordinates[k];
            }
        }
    }


    /// Compute the Gauss points
    /** For a triangle
     *
     * @see PreReseed
     */
    void ComputeGaussPointPositionsForPreReseed( Geometry< Node < 3 > >& geom,
                                                 BoundedMatrix<double, 3, 3 > & pos,
                                                 BoundedMatrix<double, 3, 3 > & N ) //2D
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


    /// Compute the Gauss points
    /** For a tetrahedron
     *
     * @see PreReseed
     */
    void ComputeGaussPointPositionsForPreReseed( Geometry< Node < 3 > >& geom,
                                                 BoundedMatrix<double, 4, 3 > & pos,
                                                 BoundedMatrix<double, 4, 4 > & N ) //3D
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

        pos=ZeroMatrix(4,3);
        for (unsigned int i=0; i!=4; i++) //going through the 4 nodes
        {
            array_1d<double, 3 > & coordinates = geom[i].Coordinates();
            for (unsigned int j=0; j!=4; j++) //going through the 4 particles
            {
                for (unsigned int k=0; k!=3; k++) //x,y,z
                    pos(j,k) += N(j,i) * coordinates[k];
            }
        }

    }


    /// Compute the Gauss points
    /**
     */
    void ComputeGaussPointPositions_45( Geometry< Node < 3 > >& geom,
                                        BoundedMatrix<double, 45, 3 > & pos,
                                        BoundedMatrix<double, 45, 3 > & N )
    {
        unsigned int counter=0;
        for (unsigned int i=0; i!=9;i++)
        {
            for (unsigned int j=0; j!=(9-i);j++)
            {
                N(counter,0)=0.05+double(i)*0.1;
                N(counter,1)=0.05+double(j)*0.1;
                N(counter,2)=1.0 - ( N(counter,1)+ N(counter,0) ) ;
                pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X();
                pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y();
                pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z();

                counter++;
            }
        }
    }


    /// Compute the Gauss points
    /**
     */
    void ComputeGaussPointPositions_initial( Geometry< Node < 3 > >& geom,
                                             BoundedMatrix<double, 15, 3 > & pos,
                                             BoundedMatrix<double, 15, 3 > & N ) //2D
    {
        unsigned int counter=0;
        for (unsigned int i=0; i!=5;i++)
        {
            for (unsigned int j=0; j!=(5-i);j++)
            {
                N(counter,0)=0.05+double(i)*0.2;
                N(counter,1)=0.05+double(j)*0.2;
                N(counter,2)=1.0 - ( N(counter,1)+ N(counter,0) ) ;
                pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X();
                pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y();
                pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z();

                counter++;
            }
        }
    }


    /// Compute the Gauss points
    /**
     */
    void ComputeGaussPointPositions_initial( Geometry< Node < 3 > >& geom,
                                             BoundedMatrix<double, 20, 3 > & pos,
                                             BoundedMatrix<double, 20, 4 > & N ) //3D
    {
        double fraction_increment;
        unsigned int counter=0;
        for (unsigned int i=0; i!=4;i++) //going to build a particle "pyramid"(tetrahedra) by layers. the first layer will be made by a triangle of 4 base X 4 height. since it is a triangle, it means it will have 10 particles
        {
            for (unsigned int j=0; j!=(4-i);j++)
            {
                for (unsigned int k=0; k!=(4-i-j);k++)
                {
                    N(counter,0)= 0.27 * ( 0.175 + double(i) ) ; //this is our "surface" in which we will build each layer, so we must construct a triangle using what's left of the shape functions total (a total of 1)

                    //total = 1.0 - N(counter,0);
                    fraction_increment = 0.27; //

                    N(counter,1)=fraction_increment * (0.175 + double(j));
                    N(counter,2)=fraction_increment * (0.175 + double(k));
                    N(counter,3)=1.0 - ( N(counter,0)+ N(counter,1) + N(counter,2) ) ;
                    pos(counter, 0) = N(counter,0) * geom[0].X() + N(counter,1) * geom[1].X() + N(counter,2) * geom[2].X() + N(counter,3) * geom[3].X();
                    pos(counter, 1) = N(counter,0) * geom[0].Y() + N(counter,1) * geom[1].Y() + N(counter,2) * geom[2].Y() + N(counter,3) * geom[3].Y();
                    pos(counter, 2) = N(counter,0) * geom[0].Z() + N(counter,1) * geom[1].Z() + N(counter,2) * geom[2].Z() + N(counter,3) * geom[3].Z();

                    counter++;
                }

            }
        }
    }


    /// check function
    virtual int Check()
    {
        KRATOS_TRY

        Node<3>& rnode = *mrModelPart.NodesBegin();

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mVectorVar1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(mScalarVar1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_VECTOR1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DELTA_SCALAR1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MEAN_SIZE, rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(YP, rnode)

        return 0;

        KRATOS_CATCH("")
    }


    /// Member variables
    ModelPart& mrModelPart;
    int mNParticles;
    int mNElems;
    int mOffset;
    int mMaxSubSteps;
    double mMaxSubStepDt;
    int mMaxNumberOfParticles;
    std::vector< ShallowParticle > mParticlesVector;
    int mLastElemId;
    bool mOddTimeStep;
    bool mParticlePrintingToolInitialized;
    unsigned int mLastNodeId;

    vector<int> mNumOfParticlesInElems;
    vector<int> mNumOfParticlesInElemsAux;
    vector<ParticlePointerVector> mVectorOfParticlePointersVectors;

    typename BinsObjectDynamic<Configure>::Pointer mpBinsObjectDynamic;

    Variable<double> mScalarVar1;
    Variable<array_1d<double,3>> mVectorVar1;
    std::string m_scalar_var1_name;
    std::string m_vector_var1_name;

}; // class MoveShallowWaterParticleUtility

}  // namespace Kratos.

#endif // KRATOS_MOVE_SHALLOW_WATER_PARTICLE_UTILITY_H_INCLUDED  defined
