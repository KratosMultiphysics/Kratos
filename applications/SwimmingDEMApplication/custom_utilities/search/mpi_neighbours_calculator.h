//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2011-6-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change
//something so that you can pass from a kdtree to a bin data search structure;

#if !defined(KRATOS_MPI_NEIGHBOURS_CALCULATOR )
#define  KRATOS_MPI_NEIGHBOURS_CALCULATOR

//M: we are using static bins for objects...

#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "containers/weak_pointer_vector.h"
#include "containers/pointer_vector.h"
#include "containers/pointer_vector_set.h"

#include "custom_utilities/search/mpi_discrete_particle_configure.h"
#include "custom_utilities/neighbours_calculator.h"

#include "../applications/mpi_search_application/custom_utilities/bins_dynamic_objects_mpi.h"
#include "../applications/mpi_search_application/custom_utilities/morton_parallel_partitioner.h"

namespace Kratos {

    template< //pot no compilar en windows aquest tipus d'assignacio per template.
    class TParticle
    >

    class Mpi_Neighbours_Calculator: public Neighbours_Calculator<TParticle> {
    public:
        typedef MpiDiscreteParticleConfigure < 3 > ConfigureType;

        typedef TParticle Particle; // es el objecte
        typedef typename Particle::Pointer ParticlePointer; // es punter al objecte
        typedef ModelPart::ElementsContainerType::ContainerType ParticleVector; // es un vector d'objectes.
        typedef ParticleVector::iterator ParticleIterator; // es un iterador d'objectes

        typedef ModelPart::ElementsContainerType ParticlePointerVector; // es un vector de iteradors
        typedef ParticlePointerVector::iterator ParticlePointerIterator; // es un iterador de punters

        typedef ConfigureType::PointType PointType;
        typedef ConfigureType::DistanceIteratorType DistanceIteratorType;
        typedef ConfigureType::ContainerType ContainerType;
        typedef ConfigureType::PointerType PointerType;
        typedef ConfigureType::IteratorType IteratorType; // iterador de punteros.
        typedef ConfigureType::ResultContainerType ResultContainerType;
        typedef ConfigureType::ResultPointerType ResultPointerType;
        typedef ConfigureType::ResultIteratorType ResultIteratorType;
        typedef ConfigureType::ContactPairType ContactPairType;
        typedef ConfigureType::ContainerContactType ContainerContactType;
        typedef ConfigureType::IteratorContactType IteratorContactType;
        typedef ConfigureType::PointerContactType PointerContactType;
        typedef ConfigureType::PointerTypeIterator PointerTypeIterator;

        typedef WeakPointerVector<Element> ParticleWeakVector;
        typedef typename ParticleWeakVector::iterator ParticleWeakIterator;
        typedef ParticleWeakVector::ptr_iterator ParticleWeakIteratorType_ptr;

        typedef std::vector<double> DistanceVector;
        typename DistanceVector::iterator DistanceIterator;

        typedef std::vector<array_1d<double, 3 > > TangDisplacementsVectorType;
        typedef TangDisplacementsVectorType::iterator TangDisplacementsIteratorType;

        typedef BinsObjectDynamicMpi <ConfigureType> Bins;
        typedef MortonParallelPartitioner <ConfigureType> Part;

        /// Pointer definition of Neighbour_calculator
        //KRATOS_CLASS_POINTER_DEFINITION(Neighbours_Calculator);  R: necesitu? C: No ho se :S

        virtual ~Mpi_Neighbours_Calculator() {
        };
          
        static void Parallel_partitioning(ModelPart& r_model_part, bool extension_option, int CalculateBoundary)
        {
            KRATOS_TRY
                      
            ContainerType& pLocalElements = r_model_part.GetCommunicator().LocalMesh().ElementsArray();
//             ContainerType& pGhostElements = r_model_part.GetCommunicator().GhostMesh().ElementsArray();

            ProcessInfo& rCurrentProcessInfo = r_model_part.GetProcessInfo();
            
            double radius_extend = 0.0;
            if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
            
            static double MaxNodeRadius = 0.0f;
            if(MaxNodeRadius == 0.0f) //TODO
                for (IteratorType particle_pointer_it = pLocalElements.begin(); particle_pointer_it != pLocalElements.end(); ++particle_pointer_it)
                {
                    double NodeRaidus = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()[0].GetSolutionStepValue(RADIUS); //TODO: MA: WATCH OUT! The use of SEARCH_RADIUS_EXTENSION might have changed
                    MaxNodeRadius = NodeRaidus > MaxNodeRadius ? NodeRaidus : MaxNodeRadius;
                }
            
            static Part partitioner;
            partitioner.LloydsBasedPartitioner(r_model_part,MaxNodeRadius,CalculateBoundary);
            
            KRATOS_CATCH("")
        }
        
        virtual void Add_To_Modelpart(ModelPart& r_model_part, ResultIteratorType neighbour_it)
        {
            Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();
            
//             ContainerType& pLocalElements = r_model_part.GetCommunicator().LocalMesh().ElementsArray();
            ContainerType& pGhostElements = r_model_part.GetCommunicator().GhostMesh().ElementsArray();
            
            int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();
            int destination = -1;
            
            bool IsInGhostMesh = false;
//             bool IsInLocalMesh = false;
          
            for(int i = 0; i < NumberOfRanks; i++)
                if((*neighbour_it)->GetGeometry()[0].GetSolutionStepValue(PARTITION_INDEX) == communicator_ranks[i])
                    destination = i;
                            
            if(destination > -1)
            {   
                for(IteratorType element_it = pGhostElements.begin(); !IsInGhostMesh && element_it != pGhostElements.end(); ++element_it)
                    if((*element_it)->GetGeometry()[0].Id() == (*neighbour_it)->GetGeometry()[0].Id())
                        IsInGhostMesh = true;
                                
//                 for(IteratorType element_it = pLocalElements.begin(); !IsInLocalMesh && element_it != pLocalElements.end(); ++element_it)
//                     if((*element_it)->GetGeometry()[0].Id() == (*neighbour_it)->GetGeometry()[0].Id())
//                         IsInLocalMesh = true;
                        
                if(!IsInGhostMesh /*&& !IsInLocalMesh*/)
                {
                    r_model_part.GetCommunicator().GhostMesh().Elements().push_back((*neighbour_it));
                    r_model_part.GetCommunicator().GhostMesh().Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                }
                
                IsInGhostMesh = false;
//                 IsInLocalMesh = false;
              
                ContainerType& pMyGhostElements = r_model_part.GetCommunicator().GhostMesh(destination).ElementsArray();
//                 ContainerType& pMyLocalElements = r_model_part.GetCommunicator().LocalMesh(destination).ElementsArray();
          
                for(IteratorType element_it = pMyGhostElements.begin(); !IsInGhostMesh && element_it != pMyGhostElements.end(); ++element_it)
                    if((*element_it)->GetGeometry()[0].Id() == (*neighbour_it)->GetGeometry()[0].Id())
                        IsInGhostMesh = true;
                           
//                 for(IteratorType element_it = pMyLocalElements.begin(); !IsInLocalMesh && element_it != pMyLocalElements.end(); ++element_it)
//                     if((*element_it)->GetGeometry()[0].Id() == (*particle_pointer_it)->GetGeometry()[0].Id())
//                         IsInLocalMesh = true;
                
                if(!IsInGhostMesh)
                {
                    r_model_part.GetCommunicator().GhostMesh(destination).Elements().push_back((*neighbour_it));
                    r_model_part.GetCommunicator().GhostMesh(destination).Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                }
                
//                if(!IsInLocalMesh)
//                {
//                    r_model_part.GetCommunicator().LocalMesh(destination).Elements().push_back((*particle_pointer_it));
//                    r_model_part.GetCommunicator().LocalMesh(destination).Nodes().push_back((*particle_pointer_it)->GetGeometry()(0));
//                }
            }
        }
                
        virtual void Clean_Modelpart(ModelPart& r_model_part)
        {
            Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();
            
            unsigned int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();
          
            ModelPart::ElementsContainerType ETempGhost[NumberOfRanks];
            ModelPart::ElementsContainerType ETempLocal[NumberOfRanks];
            ModelPart::NodesContainerType NTempGhost[NumberOfRanks];
            ModelPart::NodesContainerType NTempLocal[NumberOfRanks];
                    
            //Clean the ghost(i) and local(i) meshes
            for(unsigned int i = 0; i < NumberOfRanks; i++)
            {
                ETempGhost[i].swap(r_model_part.GetCommunicator().GhostMesh(i).Elements());
                ETempLocal[i].swap(r_model_part.GetCommunicator().LocalMesh(i).Elements());
                NTempGhost[i].swap(r_model_part.GetCommunicator().GhostMesh(i).Nodes());
                NTempLocal[i].swap(r_model_part.GetCommunicator().LocalMesh(i).Nodes());
            }
            
            //Celan the ghost mesh
            ModelPart::ElementsContainerType ETempGhostGlobal;
            ModelPart::NodesContainerType NTempGhostGlobal;
            
            ETempGhostGlobal.swap(r_model_part.GetCommunicator().GhostMesh().Elements());
            NTempGhostGlobal.swap(r_model_part.GetCommunicator().GhostMesh().Nodes());
        }
        
        virtual void Sort_Modelpart(ModelPart& r_model_part)
        {
            for (unsigned int i = 0; i < r_model_part.GetCommunicator().LocalMeshes().size(); i++)
                r_model_part.GetCommunicator().LocalMesh(i).Nodes().Unique();
            
            for (unsigned int i = 0; i < r_model_part.GetCommunicator().GhostMeshes().size(); i++)
                r_model_part.GetCommunicator().GhostMesh(i).Nodes().Unique();
        }
        
        virtual ContainerType& Get_Elements(ModelPart& r_model_part)
        {
            return r_model_part.GetCommunicator().LocalMesh().ElementsArray();
        }
        
        virtual void SearchNeighbours(ModelPart& r_model_part,
                                      ContainerType& pIteratorElements,
                                      int NumberOfElements,
                                      int MaximumNumberOfResults,
                                      std::vector<std::size_t> &NumberOfResults, 
                                      std::vector<std::vector<PointerType> > &Results,
                                      std::vector<std::vector<double> > &ResultsDistances,
                                      std::vector<double> &Radius
        )
        {
            Bins particle_bin(pIteratorElements.begin(), pIteratorElements.end());
            particle_bin.SearchObjectsMpi(r_model_part,pIteratorElements.begin(),NumberOfElements,Radius,Results,ResultsDistances,NumberOfResults,MaximumNumberOfResults,r_model_part.pGetCommunicator());
        }

    protected:


    private:

        inline void Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
            unsigned int buffer_size = node_it->GetBufferSize();
            for (unsigned int step = 0; step < buffer_size; step++) {
                //getting the data of the solution step
                double* step_data = (node_it)->SolutionStepData().Data(step);
                //copying this data in the position of the vector we are interested in
                for (int j = 0; j < step_data_size; j++) {
                    step_data[j] = 0.0;

                }

            }
        }

        inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
            /*
            array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
            noalias(Aux_var) = ZeroVector(3);
             * */
        }

        inline void ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {
            /*
            double& Aux_var = (particle_it->GetPointerToCenterNode()).FastGetSolutionStepValue(rVariable, 0);
            Aux_var = 0.0;
             * */
        }

        Mpi_Neighbours_Calculator & operator=(Mpi_Neighbours_Calculator const& rOther);

    }; // Class Neighbours_calculator


} // namespace Kratos.

#endif // KRATOS_MPI_NEIGHBOURS_CALCULATOR  defined 



