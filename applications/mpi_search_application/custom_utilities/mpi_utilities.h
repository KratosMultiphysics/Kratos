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
      
      typedef SpatialSearch                                         SearchType;
      typedef MpiDiscreteParticleConfigure<3>                       Configure;

      typedef SearchType::ElementsContainerType::ContainerType      ElementsContainerType;
      typedef SearchType::NodesContainerType::ContainerType         NodesContainerType;
      typedef SearchType::ConditionsContainerType::ContainerType    ConditionsContainerType;

      
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
      
            
      void TransferModelNodes(ModelPart& rModelPart)
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
          
          NodesContainerType& pNodes = rModelPart.GetCommunicator().LocalMesh().NodesArray();
          
          //Fill the buffer with elements to be transfered
          for (NodesContainerType::iterator i_node = pNodes.begin(); i_node != pNodes.end(); ++i_node)
          {
              int PartitionIndex = (*i_node)->GetSolutionStepValue(PARTITION_INDEX);
              if(PartitionIndex != mpi_rank)
              {
                  SendObjects[PartitionIndex].push_back(*i_node);
              }
          }
          
          rModelPart.GetCommunicator().TransferObjects(SendObjects,RecvObjects);
          BuildNewNodesPartitions(rModelPart,RecvObjects);
          
          KRATOS_CATCH("")
      }
      
      void TransferModelElements(ModelPart& rModelPart)
      {
          KRATOS_TRY
          
          int mpi_rank;
          int mpi_size;
          
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
          
          std::vector<ElementsContainerType> SendObjects(mpi_size);
          std::vector<ElementsContainerType> RecvObjects(mpi_size);
          
          for(int i = 0; i < mpi_size; i++)
          {
              SendObjects[i].reserve(rModelPart.GetCommunicator().LocalMesh().NumberOfElements());
          }
          
          ElementsContainerType& pElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray();
          
          //Fill the buffer with elements to be transfered
          for (ElementsContainerType::iterator i_element = pElements.begin(); i_element != pElements.end(); ++i_element)
          {
              int PartitionIndex = (*i_element)->GetValue(PARTITION_INDEX);
              
              if(PartitionIndex != mpi_rank)
              {
                  SendObjects[PartitionIndex].push_back(*i_element);
              }
          }
          
          rModelPart.GetCommunicator().TransferObjects(SendObjects,RecvObjects);
          BuildNewElementsPartitions(rModelPart,RecvObjects);
          
          KRATOS_CATCH("")
      }
    
      void ParallelPartitioning(ModelPart& rModelPart, bool extension_option, int CalculateBoundry)
      {
          KRATOS_TRY
                    
          ElementsContainerType& pLocalElements = rModelPart.GetCommunicator().LocalMesh().ElementsArray();

          ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
          
          double radius_extend = 0.0;
          if (extension_option) radius_extend = rCurrentProcessInfo[SEARCH_RADIUS_EXTENSION];
          
          static double MaxNodeRadius = 0.0f;
          if(MaxNodeRadius == 0.0f) //TODO
              for (ElementsContainerType::iterator particle_pointer_it = pLocalElements.begin(); particle_pointer_it != pLocalElements.end(); ++particle_pointer_it)
              {
                  double NodeRaidus = (1.0 + radius_extend) * (*particle_pointer_it)->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                  MaxNodeRadius = NodeRaidus > MaxNodeRadius ? NodeRaidus : MaxNodeRadius;
              }
          
          LloydParallelPartitioner<Configure> partitioner;
          partitioner.LloydsBasedParitioner(rModelPart,MaxNodeRadius,CalculateBoundry);
          TransferModelElements(rModelPart);
          partitioner.CalculatePartitionInterface(rModelPart);
          
          Timer::SetOuputFile("TimesPartitioner");
          Timer::PrintTimingInformation();
          
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
          
          for (ElementsContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Elements().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Elements().ptr_end(); ++it)
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
          
          std::cout << "STARTING NODE ID: " << iteratorId << std::endl;
          std::cout << "ENDING   NODE ID: " << iteratorId + nodes_size << std::endl;

          for (NodesContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Nodes().ptr_end(); ++it)
          {
              (*it)->SetId(iteratorId);
              iteratorId++;
          }
          
          std::cout << "New ids assigned" << std::endl;
          
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
          
          for (ConditionsContainerType::iterator it = mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_begin(); it != mModelPart.GetCommunicator().LocalMesh().Conditions().ptr_end(); ++it)
              (*it)->SetId(iteratorId++);
          
          KRATOS_CATCH("")
      }
      
      bool BuildNewElementsPartitions(ModelPart& mModelPart, std::vector<ElementsContainerType> &RecvObjects)
      { 
          int mpi_rank;
          int mpi_size;
          
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
          
          ElementsContainerType& ElementsLocal  = mModelPart.GetCommunicator().LocalMesh().ElementsArray();
          ElementsContainerType& ElementsGlobal = mModelPart.ElementsArray();
          
          ElementsContainerType temp_particles_container_local;
          ElementsContainerType temp_particles_container_global;
          
          NodesContainerType& NodesLocal  = mModelPart.GetCommunicator().LocalMesh().NodesArray();
          NodesContainerType& NodesGlobal = mModelPart.NodesArray();
          
          NodesContainerType temp_nodes_container_local;
          NodesContainerType temp_nodes_container_global;
          
          temp_particles_container_local.reserve(ElementsLocal.size());
          temp_particles_container_global.reserve(ElementsGlobal.size());
          
          temp_nodes_container_local.reserve(NodesLocal.size());
          temp_nodes_container_global.reserve(NodesGlobal.size());

          temp_particles_container_local.swap(ElementsLocal);
          temp_particles_container_global.swap(ElementsGlobal);
          
          temp_nodes_container_local.swap(NodesLocal);
          temp_nodes_container_global.swap(NodesGlobal);
          
          // Keep Global elements and nodes from our domain
          for (ElementsContainerType::iterator i_element = temp_particles_container_global.begin();
              i_element != temp_particles_container_global.end(); ++i_element)
          {
              if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.Elements().push_back((*i_element));
              }
          }
          
          for (NodesContainerType::iterator i_node = temp_nodes_container_global.begin();
              i_node != temp_nodes_container_global.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.Nodes().push_back((*i_node));
              }
          }

          // Keep Local elements and nodes from our domain
          for (ElementsContainerType::iterator i_element = temp_particles_container_local.begin();
              i_element != temp_particles_container_local.end(); ++i_element)
          {
              if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.GetCommunicator().LocalMesh().Elements().push_back((*i_element));
              }
          }
          
          for (NodesContainerType::iterator i_node = temp_nodes_container_local.begin();
              i_node != temp_nodes_container_local.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.GetCommunicator().LocalMesh().Nodes().push_back((*i_node));
              }
          }

          // Add new elements and nodes
          for(int i = 0; i < mpi_size; i++)
          {   
              for (ElementsContainerType::iterator i_element = RecvObjects[i].begin();
              i_element != RecvObjects[i].end(); ++i_element)
              {
                  if((*i_element)->GetValue(PARTITION_INDEX) == mpi_rank)
                  {                 
                      mModelPart.Elements().push_back(*i_element);
                      mModelPart.GetCommunicator().LocalMesh().Elements().push_back(*i_element);
                      
                      for (unsigned int i = 0; i < (*i_element)->GetGeometry().PointsNumber(); i++)
                      {
                          ModelPart::NodeType::Pointer pNode = (*i_element)->GetGeometry().pGetPoint(i);
                          
                          pNode->GetSolutionStepValue(PARTITION_INDEX) = mpi_rank;
                          
                          mModelPart.Nodes().push_back(pNode);
                          mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(pNode);
                      }
                  }
              }
          }
          
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
          
          return true;
      }
      
      bool BuildNewNodesPartitions(ModelPart& mModelPart, std::vector<NodesContainerType> &RecvObjects)
      {
          int mpi_rank;
          int mpi_size;
          
          MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
          MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
          
          NodesContainerType& NodesLocal  = mModelPart.GetCommunicator().LocalMesh().NodesArray();
          NodesContainerType& NodesGlobal = mModelPart.NodesArray();
    
          NodesContainerType temp_nodes_container_local;
          NodesContainerType temp_nodes_container_global;
          
          temp_nodes_container_local.reserve(NodesLocal.size());
          temp_nodes_container_global.reserve(NodesGlobal.size());

          temp_nodes_container_local.swap(NodesLocal);
          temp_nodes_container_global.swap(NodesGlobal);
          
          // Keep Global and nodes from our domain
          for (NodesContainerType::iterator i_node = temp_nodes_container_global.begin();
              i_node != temp_nodes_container_global.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.Nodes().push_back((*i_node));
              }
          }

          // Keep Local nodes from our domain
          for (NodesContainerType::iterator i_node = temp_nodes_container_local.begin();
              i_node != temp_nodes_container_local.end(); ++i_node)
          {
              if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
              {                 
                  mModelPart.GetCommunicator().LocalMesh().Nodes().push_back((*i_node));
              }
          }

          // Add new elements and nodes
          for(int i = 0; i < mpi_size; i++)
          {
              for (NodesContainerType::iterator i_node = RecvObjects[i].begin();
                  i_node != RecvObjects[i].end(); ++i_node)
              {
                  if((*i_node)->GetSolutionStepValue(PARTITION_INDEX) == mpi_rank)
                  {                 
                      mModelPart.Nodes().push_back(*i_node);
                      mModelPart.GetCommunicator().LocalMesh().Nodes().push_back(*i_node);
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




