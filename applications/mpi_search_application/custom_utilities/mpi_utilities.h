//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Carlos A. Roig $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1 $
//
//

#if !defined(KRATOS_MPI_UTILITIES)
#define  KRATOS_MPI_UTILITIES


// System includes

// Project includes
#include "utilities/timer.h"

/* System includes */
#include <limits>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/line_3d_2.h"

/* Search */
#include "custom_utilities/lloyd_parallel_partitioner.h"
#include "processes/graph_coloring_process.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// Short class definition.
  /** Detail class definition.
  */

  class MpiUtilities
  {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MpiUtilities
      KRATOS_CLASS_POINTER_DEFINITION(MpiUtilities);

      typedef SpatialSearch                          SearchType;
      typedef MpiDiscreteParticleConfigure<3>        Configure;

      typedef SearchType::ElementsContainerType      ElementsContainerType;
      typedef SearchType::NodesContainerType         NodesContainerType;
      typedef SearchType::ConditionsContainerType    ConditionsContainerType;

      typedef GraphColoringProcess::GraphType        GraphType;


      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MpiUtilities(){}

      MpiUtilities(ModelPart& model_part,
                   ModelPart& contacts_model_part
      )
      {
      }

      /// Destructor.
      virtual ~MpiUtilities(){}

      /**
       *    Migrate Scheme:
       *
       *                            X(i) stands for:
       *                    Object X has PARTITION_INDEX = i
       *
       *    Initial ModelParts                  Final ModelParts
       *    in each process                     in each process
       *
       *    Process-0      Process-1            Process-0      Process-1
       *    +----------+   +----------+         +----------+   +----------+
       *    |  A(0)    |   |  D(0)    |         |  A(0)    |   |  B(1)    |
       *    |  B(1)    |   |  E(1)    |         |  D(0)    |   |  E(1)    |
       *    |  C(2)    |   |  F(2)    |         |  J(0)    |   |  G(1)    |
       *    |          |   |          |         |          |   |  K(1)    |
       *    +----------+   +----------+         +----------+   +----------+
       *
       *    Process-2      Process-3            Process-2      Process-3
       *    +----------+   +----------+         +----------+   +----------+
       *    |  G(1)    |   |  J(0)    |         |  C(2)    |   |  I(3)    |
       *    |  H(2)    |   |  K(1)    |         |  F(2)    |   |  L(3)    |
       *    |  I(3)    |   |  L(2)    |         |  H(2)    |   |          |
       *    |          |   |  M(3)    |         |  L(2)    |   |          |
       *    +----------+   +----------+         +----------+   +----------+
       *
      **/

      /**
       * MigrateNodes for an specific meshgroup.
      **/
      void MigrateMeshNodes(ModelPart& rModelPart, NodesContainerType::ContainerType& pNodes, std::vector<NodesContainerType>& SendObjects, std::vector<NodesContainerType> RecvObjects, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          //Fill the send buffer with elements to be transfered
          for (NodesContainerType::ContainerType::iterator i_node = pNodes.begin(); i_node != pNodes.end(); ++i_node)
          {
              int PartitionIndex = (*i_node)->GetSolutionStepValue(PARTITION_INDEX);

              if(PartitionIndex != mpi_rank)
              {
                  SendObjects[PartitionIndex].push_back(*i_node);
              }
          }

          rModelPart.GetCommunicator().TransferObjects(SendObjects,RecvObjects);
          BuildNewNodesPartitions(rModelPart,RecvObjects,groupId);

          //Clear the buffers
          SendObjects.clear();
          RecvObjects.clear();
      }

      /**
       * Transfer all nodes in a given ModelPart to the partition indicated by the PARTITION_INDEX variable of
       * every node in that ModelPart
       * @param Modelpart: Input ModelPart
       **/
      void MigrateNodes(ModelPart& rModelPart)
      {
          KRATOS_TRY

          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          std::vector<NodesContainerType> SendObjects(mpi_size);
          std::vector<NodesContainerType> RecvObjects(mpi_size);

          for(int i = 0; i < mpi_size; i++)
          {
              SendObjects[i].reserve(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes());
          }

          NodesContainerType::ContainerType& pNodes = rModelPart.GetCommunicator().LocalMesh().NodesArray();

          MigrateMeshNodes(rModelPart,pNodes,SendObjects,RecvObjects,0);
          FinalizeNewPartition(rModelPart);

          KRATOS_CATCH("")
      }

      /**
       * MigrateElements for an specific meshgroup.
       **/
      void MigrateMeshElements(ModelPart& rModelPart, std::vector<ElementsContainerType>& SendObjects, std::vector<ElementsContainerType> RecvObjects, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          ElementsContainerType::ContainerType& pElements = rModelPart.GetMesh(groupId).ElementsArray();

          //Fill the send buffer with elements to be transfered
          for (ElementsContainerType::ContainerType::iterator i_element = pElements.begin(); i_element != pElements.end(); ++i_element)
          {
              int PartitionIndex = (*i_element)->GetValue(PARTITION_INDEX);

              if(PartitionIndex != mpi_rank)
              {
                  SendObjects[PartitionIndex].push_back(*i_element);
              }
          }

          rModelPart.GetCommunicator().TransferObjects(SendObjects,RecvObjects);
          BuildNewElementsPartitions(rModelPart,RecvObjects,groupId);

          //Clear the buffers
          SendObjects.clear();
          RecvObjects.clear();
      }

      void MigrateMeshElementsId(ModelPart& rModelPart, std::vector<std::vector<int> >& SendObjectsId, std::vector<std::vector<int> >& RecvObjectsId, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          NodesContainerType::ContainerType& pNodes = rModelPart.GetMesh(groupId).NodesArray();

          // Fill the send buffer with elements to be transfered
          for (NodesContainerType::ContainerType::iterator i_node = pNodes.begin(); i_node != pNodes.end(); ++i_node)
          {
              int PartitionIndex = (*i_node)->GetSolutionStepValue(PARTITION_INDEX);

              if(PartitionIndex != mpi_rank)
              {
                  SendObjectsId[PartitionIndex].push_back((*i_node)->Id());
              }
          }

          std::stringstream * serializer_buffer;
          std::vector<std::string> buffer(mpi_size);

          int * msgSendSize = new int[mpi_size];
          int * msgRecvSize = new int[mpi_size];

          for(int i = 0; i < mpi_size; i++)
          {
              msgSendSize[i] = 0;
              msgRecvSize[i] = 0;
          }

          for(int i = 0; i < mpi_size; i++)
          {
              if(mpi_rank != i)
              {
                  Kratos::Serializer particleSerializer;
                  particleSerializer.save("nodes",SendObjectsId[i]);

                  serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                  buffer[i] = std::string(serializer_buffer->str());
                  msgSendSize[i] = buffer[i].size()+1;
              }
          }

          MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);

          int NumberOfCommunicationEvents = 0;
          int NumberOfCommunicationEventsIndex = 0;

          char ** message = new char * [mpi_size];
          char ** mpi_send_buffer = new char *[mpi_size];

          for(int j = 0; j < mpi_size; j++)
          {
              if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents++;
              if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents++;
          }

          MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
          MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

          //Set up all receive and send events
          for(int i = 0; i < mpi_size; i++)
          {
              if((i != mpi_rank) && msgRecvSize[i])
              {
                  message[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

                  MPI_Irecv(message[i],msgRecvSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }

              if((i != mpi_rank) && msgSendSize[i])
              {
                  mpi_send_buffer[i] = (char *)malloc(sizeof(char) * msgSendSize[i]);
                  memcpy(mpi_send_buffer[i],buffer[i].c_str(),msgSendSize[i]);

                  MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
              }
          }

          //wait untill all communications finish
          MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);
          MPI_Barrier(MPI_COMM_WORLD);

          for(int i = 0; i < mpi_size; i++)
          {
              if (i != mpi_rank && msgRecvSize[i])
              {
                  Kratos::Serializer particleSerializer;
                  std::stringstream * serializer_buffer;
                  serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
                  serializer_buffer->write(message[i], msgRecvSize[i]);
                  particleSerializer.load("nodes",RecvObjectsId[i]);
              }

              MPI_Barrier(MPI_COMM_WORLD);
          }

          BuildNewMeshPartitions(rModelPart,RecvObjectsId,groupId);

          // Clear the buffers
          for(int i = 0; i < mpi_size; i++)
          {
              SendObjectsId[i].clear();
              RecvObjectsId[i].clear();
          }

          delete [] reqs;
          delete [] stats;

          delete [] message;
          delete [] mpi_send_buffer;

          delete [] msgSendSize;
          delete [] msgRecvSize;
      }

      /**
      * Transfer all elemens and nodes in a given ModelPart to the partition indicated by the PARTITION_INDEX variable of
      * every element in that ModelPart
      * @param rModelPart: Input ModelPart
      **/
      void MigrateElements(ModelPart& rModelPart)
      {
          // KRATOS_TRY

          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          std::vector<ElementsContainerType> SendObjects(mpi_size);
          std::vector<ElementsContainerType> RecvObjects(mpi_size);

          std::vector<std::vector<int> >  SendObjectsId(mpi_size, std::vector<int>(0));
          std::vector<std::vector<int> >  RecvObjectsId(mpi_size, std::vector<int>(0));

          for(int i = 0; i < mpi_size; i++)
          {
              SendObjects[i].reserve(rModelPart.GetCommunicator().LocalMesh().NumberOfElements());
              SendObjectsId[i].reserve(rModelPart.GetCommunicator().LocalMesh().NumberOfElements());
          }

          // Main Mesh
          MigrateMeshElements(rModelPart,SendObjects,RecvObjects,0);
          FinalizeNewPartition(rModelPart);

          // Rest of the meshes
          for(unsigned int i = 1; i < rModelPart.NumberOfMeshes(); i++)
          {
              MigrateMeshElementsId(rModelPart,SendObjectsId,RecvObjectsId,i);
          }

          // KRATOS_CATCH("")
      }

      void ParallelPartitioning(ModelPart& rModelPart, bool extension_option, int CalculateBoundary)
      {
          KRATOS_TRY

          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          typedef GraphColoringProcess::GraphType GraphType;
          GraphType domains_graph(mpi_size, mpi_size, 0);
          GraphType global_domains_graph(mpi_size, mpi_size, 0);
          GraphType domains_colored_graph;

          ElementsContainerType::ContainerType& pLocalElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray();

          //ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

          //double search_increment = 0.0;
          //if (extension_option) search_increment = rCurrentProcessInfo[SEARCH_RADIUS_INCREMENT];

          static double MaxNodeRadius = 0.0f;
          if(MaxNodeRadius == 0.0f) {//TODO
            for (ElementsContainerType::ContainerType::iterator particle_pointer_it = pLocalElements.begin(); particle_pointer_it != pLocalElements.end(); ++particle_pointer_it) {
              double NodeRadius = Configure::GetObjectRadius(*particle_pointer_it);
              //double NodeRadius = search_increment + (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
              MaxNodeRadius = NodeRadius > MaxNodeRadius ? NodeRadius : MaxNodeRadius;
            }
          }

          int colors_number;

          LloydParallelPartitioner<Configure> partitioner(pLocalElements.begin(), pLocalElements.end());

          if(mpi_rank == 0) {
            KRATOS_WATCH(domains_graph)
          }

          partitioner.SerialPartition();
          MigrateElements(rModelPart);
          partitioner.UpdateDomainGraph(pLocalElements.begin(), pLocalElements.end(), domains_graph);

          MPI_Allreduce(&domains_graph(0,0), &global_domains_graph(0, 0), mpi_size*mpi_size, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          GraphColoringProcess(mpi_size, global_domains_graph, domains_colored_graph, colors_number).Execute();

          //allocate space needed in the communicator
          rModelPart.GetCommunicator().SetNumberOfColors(colors_number);
          rModelPart.GetCommunicator().NeighbourIndices().resize(colors_number);

          // Resize the neighbour index vector with the new value
          Communicator::NeighbourIndicesContainerType& neighbours_indices = rModelPart.GetCommunicator().NeighbourIndices();
          if(neighbours_indices.size() != static_cast<unsigned int>(colors_number)) {
            neighbours_indices.resize(colors_number,false);
          }

          // Reset the list of the neighbour numbers with the ones calculated in global_domains_graph
          for(int i = 0; i  < colors_number; i++) {
            neighbours_indices[i] = domains_colored_graph(mpi_rank,i);
          }

          // Adding local, ghost and interface meshes to ModelPart if is necessary
          int number_of_meshes =  ModelPart::Kratos_Ownership_Size + colors_number; // (all + local + ghost) + (colors_number for interfaces)
          if(rModelPart.GetMeshes().size() < static_cast<unsigned int>(number_of_meshes)) {
            for(int i = rModelPart.GetMeshes().size() ; i < number_of_meshes ; i++) {
              rModelPart.GetMeshes().push_back(boost::make_shared<ModelPart::MeshType>());
            }
          }

          KRATOS_CATCH("")
      }

      void CalculateModelNewIds(ModelPart& mModelPart, int offset)
      {
          CalculateElementsNewId(mModelPart,offset);
          CalculateNodesNewId(mModelPart,offset);
          CalculateConditionsNewId(mModelPart,offset);
      }

      void CalculateElementsNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY

          int element_size = mModelPart.GetCommunicator().LocalMesh().Elements().size();
          int iteratorId = -1;

          Configure::ReduceIds(element_size,offset,iteratorId);

          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;

          for (ElementsContainerType::ContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Elements().ptr_end(); ++it)
              (*it)->SetId(iteratorId++);

          KRATOS_CATCH("")
      }

      void CalculateConditionsNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY

          int conditions_size = mModelPart.GetCommunicator().LocalMesh().Conditions().size();
          int iteratorId = -1;

          Configure::ReduceIds(conditions_size,offset,iteratorId);

          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;

          for (ConditionsContainerType::ContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_end(); ++it)
              (*it)->SetId(iteratorId++);

          KRATOS_CATCH("")
      }

      void CalculateNodesNewId(ModelPart& mModelPart, int offset)
      {
          KRATOS_TRY

          int nodes_size = mModelPart.GetCommunicator().LocalMesh().Nodes().size();
          int iteratorId = -1;

          Configure::ReduceIds(nodes_size,offset,iteratorId);

          if(iteratorId == -1)
              std::cout << "Invalid starting Id" << std::endl;

          for (NodesContainerType::ContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_end(); ++it)
          {
              (*it)->SetId(iteratorId);
              iteratorId++;
          }

          KRATOS_CATCH("")
      }

      /**
       * Cleans up the partition in order to be able to recieve new elements from other partitions.
       * This function must not be called before calculating new partition id
       * @param rModelPart: Input ModelPart
       */
      bool FinalizeNewPartition(ModelPart& mModelPart)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          // Elements
          // Global
          ElementsContainerType::ContainerType& ElementsGlobal  = mModelPart.ElementsArray();
          ElementsContainerType::ContainerType  temp_Elements;

          temp_Elements.reserve(ElementsGlobal.size());
          temp_Elements.swap(ElementsGlobal);

          for (ElementsContainerType::ContainerType::iterator i_element = temp_Elements.begin();
              i_element != temp_Elements.end(); ++i_element)
          {
              if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.Elements().push_back((*i_element));
              }
          }

          // Local
          ElementsContainerType::ContainerType& ElementsLocal   = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
          ElementsContainerType::ContainerType  temp_ElementsLocal;

          temp_ElementsLocal.reserve(ElementsLocal.size());
          temp_ElementsLocal.swap(ElementsLocal);

          for (ElementsContainerType::ContainerType::iterator i_element = temp_ElementsLocal.begin();
              i_element != temp_ElementsLocal.end(); ++i_element)
          {
              if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.GetCommunicator().LocalMesh().Elements().push_back((*i_element));
              }
          }

          // MeshGroups
          for( unsigned int meshId = 1; meshId < mModelPart.NumberOfMeshes(); meshId++)
          {
              ElementsContainerType::ContainerType& ElementsMesh  = mModelPart.GetMesh(meshId).ElementsArray();
              ElementsContainerType::ContainerType  temp_ElementsMesh;

              temp_ElementsMesh.reserve(ElementsMesh.size());
              temp_ElementsMesh.swap(ElementsMesh);

              for (ElementsContainerType::ContainerType::iterator i_element = temp_ElementsMesh.begin();
                  i_element != temp_ElementsMesh.end(); ++i_element)
              {
                  if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
                  {
                      mModelPart.GetMesh(meshId).Elements().push_back((*i_element));
                  }
              }
          }

          // Conditions
          // Global
          ConditionsContainerType::ContainerType& ConditionsGlobal  = mModelPart.ConditionsArray();
          ConditionsContainerType::ContainerType  temp_Condition;

          temp_Condition.reserve(ConditionsGlobal.size());
          temp_Condition.swap(ConditionsGlobal);

          for (ConditionsContainerType::ContainerType::iterator i_condition = temp_Condition.begin();
              i_condition != temp_Condition.end(); ++i_condition)
          {
              if((*i_condition)->GetValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.Conditions().push_back((*i_condition));
              }
          }

          // Local
          ConditionsContainerType::ContainerType& ConditionsLocal   = mModelPart.GetCommunicator().LocalMesh().ConditionsArray();
          ConditionsContainerType::ContainerType  temp_ConditionLocal;

          temp_ConditionLocal.reserve(ConditionsLocal.size());
          temp_ConditionLocal.swap(ConditionsLocal);

          for (ConditionsContainerType::ContainerType::iterator i_condition = temp_ConditionLocal.begin();
              i_condition != temp_ConditionLocal.end(); ++i_condition)
          {
              if((*i_condition)->GetValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.GetCommunicator().LocalMesh().Conditions().push_back((*i_condition));
              }
          }

          // MeshGroups
          for( unsigned int meshId = 1; meshId < mModelPart.NumberOfMeshes(); meshId++)
          {
              ConditionsContainerType::ContainerType& ConditionsMesh  = mModelPart.GetMesh(meshId).ConditionsArray();
              ConditionsContainerType::ContainerType  temp_ConditionMesh;

              temp_ConditionMesh.reserve(ConditionsMesh.size());
              temp_ConditionMesh.swap(ConditionsMesh);

              for (ConditionsContainerType::ContainerType::iterator i_condition = temp_ConditionMesh.begin();
                  i_condition != temp_ConditionMesh.end(); ++i_condition)
              {
                  if((*i_condition)->GetValue(PARTITION_INDEX) == mpi_rank)
                  {
                      mModelPart.GetMesh(meshId).Conditions().push_back((*i_condition));
                  }
              }
          }

          // Nodes
          // Global
          NodesContainerType::ContainerType&    NodesGlobal     = mModelPart.NodesArray();
          NodesContainerType::ContainerType     temp_Nodes;

          temp_Nodes.reserve(NodesGlobal.size());
          temp_Nodes.swap(NodesGlobal);

          for (NodesContainerType::ContainerType::iterator i_node = temp_Nodes.begin();
              i_node != temp_Nodes.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.Nodes().push_back((*i_node));
              }
          }

          // Local
          NodesContainerType::ContainerType&    NodesLocal      = mModelPart.GetCommunicator().LocalMesh().NodesArray();
          NodesContainerType::ContainerType     temp_NodesLocal;

          temp_NodesLocal.reserve(NodesLocal.size());
          temp_NodesLocal.swap(NodesLocal);

          for (NodesContainerType::ContainerType::iterator i_node = temp_NodesLocal.begin();
              i_node != temp_NodesLocal.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {
                  mModelPart.GetCommunicator().LocalMesh().Nodes().push_back((*i_node));
              }
          }

//           // MeshGroups
//           for( unsigned int meshId = 1; meshId < mModelPart.NumberOfMeshes(); meshId++)
//           {
//               NodesContainerType::ContainerType& NodesMesh  = mModelPart.GetMesh(meshId).NodesArray();
//               NodesContainerType::ContainerType  temp_NodesMesh;
//
//               temp_NodesMesh.reserve(NodesMesh.size());
//               temp_NodesMesh.swap(NodesMesh);
//
//               for (NodesContainerType::ContainerType::iterator i_node = temp_NodesMesh.begin();
//                   i_node != temp_NodesMesh.end(); ++i_node)
//               {
//                   if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
//                   {
//                       mModelPart.GetMesh(meshId).Nodes().push_back((*i_node));
//                   }
//               }
//           }

          // Sort both the elements and nodes of the modelpart. Otherwise the results will be unpredictable
          mModelPart.Elements().Unique();
          mModelPart.Nodes().Unique();

          // Sort both the elements and nodes of the modelpart. Otherwise the results will be unpredictable
          mModelPart.GetCommunicator().LocalMesh().Elements().Unique();
          mModelPart.GetCommunicator().LocalMesh().Nodes().Unique();

          for (unsigned int i = 0; i < mModelPart.GetCommunicator().LocalMeshes().size(); i++)
          {
              mModelPart.GetCommunicator().LocalMesh(i).Elements().Unique();
              mModelPart.GetCommunicator().LocalMesh(i).Nodes().Unique();
          }

          for (unsigned int i = 0; i < mModelPart.GetCommunicator().GhostMeshes().size(); i++)
          {
              mModelPart.GetCommunicator().GhostMesh(i).Elements().Unique();
              mModelPart.GetCommunicator().GhostMesh(i).Nodes().Unique();
          }

          return false;
      }

      /**
       * Fill a paritition with a given set of elements
       * @param rModelPart: Input ModelPart
       * @param RecvObjects: Additional elements to be added in this paritition
       **/
      bool BuildNewElementsPartitions(ModelPart& mModelPart, std::vector<ElementsContainerType> &RecvObjects, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          // Add new elements and nodes
          for(int i = 0; i < mpi_size; i++)
          {
              for (ElementsContainerType::iterator i_element = RecvObjects[i].begin();
              i_element != RecvObjects[i].end(); ++i_element)
              {

                  if(i_element->GetValue(PARTITION_INDEX) == mpi_rank)
                  {
                      mModelPart.Elements().push_back(*i_element.base());

                      if(groupId == 0)
                      {
                          // If it is the first transfer:
                          //  - add the element to the 0 mesh
                          //  - add the element to local mesh
                          //  - add the element nodes to the modelpart
                          mModelPart.GetCommunicator().LocalMesh().Elements().push_back(*i_element.base());

                          for (unsigned int i = 0; i < i_element->GetGeometry().PointsNumber(); i++)
                          {
                              ModelPart::NodeType::Pointer pNode = i_element->GetGeometry().pGetPoint(i);

                              pNode->GetSolutionStepValue(PARTITION_INDEX) = mpi_rank;

                              mModelPart.Nodes().push_back(pNode);
                              mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(pNode);
                          }
                      }
                      else  // WARNING: I'm assuming that the @ for a given element is kept between transfers, wich is most likely to be incorrect...
                      {
                          // If the group its different form 0, its not the first transfer, hence the element is already in the mesh
                          //  - add the element to the GroupId mesh
                          mModelPart.GetMesh(groupId).Elements().push_back(*i_element.base());
                      }
                  }
              }
          }

          return true;
      }

            /**
       * Fill a paritition with a given set of elements
       * @param rModelPart: Input ModelPart
       * @param RecvObjects: Additional elements to be added in this paritition
       **/
      bool BuildNewMeshPartitions(ModelPart& mModelPart, std::vector<std::vector<int> > &RecvObjects, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          std::cout << "Be Nodesize: " << mModelPart.GetMesh(groupId).Nodes().size() << std::endl;
          std::cout << "Lo Nodesize: " << mModelPart.GetCommunicator().LocalMesh().Nodes().size() << std::endl;

          // Add new elements and nodes
          for(int i = 0; i < mpi_size; i++)
          {
              for( unsigned int j = 0; j < RecvObjects[i].size(); j++)
              {
                  mModelPart.GetMesh(groupId).Nodes().push_back(mModelPart.pGetNode(RecvObjects[i][j]));
              }
          }

          std::cout << "Af Nodesize: " << mModelPart.GetMesh(groupId).Nodes().size() << std::endl;

          return true;
      }

      /**
       * Build a paritition with a given set of nodes
       * @param rModelPart: Input ModelPart
       * @param RecvObjects: Additional nodes to be added in this paritition
       **/
      bool BuildNewNodesPartitions(ModelPart& mModelPart, std::vector<NodesContainerType> &RecvObjects, const int groupId)
      {
          int mpi_rank;
          int mpi_size;

          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

          // Add new  and nodes
          for(int i = 0; i < mpi_size; i++)
          {
              for (NodesContainerType::iterator i_node = RecvObjects[i].begin();
                  i_node != RecvObjects[i].end(); ++i_node)
              {

                  if(i_node->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
                  {
                      mModelPart.Nodes().push_back(*i_node.base());

                      if(groupId == 0)
                      {
                          // If it is the first transfer:
                          //  - add the node to the 0 mesh
                          //  - add the node to local mesh
                          //  - add the node nodes to the modelpart
                          mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(*i_node.base());
                      }
                      else  // WARNING: I'm assuming that the @ for a given node is kept between transfers, wich is most likely to be incorrect...
                      {
                          // If the group its different form 0, its not the first transfer, hence the node is already in the mesh
                          //  - add the node to the GroupId mesh
                          mModelPart.GetMesh(groupId).Nodes().push_back(*i_node.base());
                      }
                  }
              }
          }

          // Sort both the nodes of the modelpart. Otherwise the results will be unpredictable
          mModelPart.GetCommunicator().LocalMesh().Nodes().Unique();

          for (unsigned int i = 0; i < mModelPart.GetCommunicator().LocalMeshes().size(); i++)
          {
              mModelPart.GetCommunicator().LocalMesh(i).Nodes().Unique();
          }

          for (unsigned int i = 0; i < mModelPart.GetCommunicator().GhostMeshes().size(); i++)
          {
              mModelPart.GetCommunicator().GhostMesh(i).Nodes().Unique();
          }

          return true;
      }

    protected:

    private:

  }; // Class MpiUtilities

 /*
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                    MpiUtilities& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                    const MpiUtilities& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    */
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
