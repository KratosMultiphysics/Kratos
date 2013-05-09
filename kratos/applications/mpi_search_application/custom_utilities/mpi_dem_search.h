//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_MPI_DEM_SEARCH_H_INCLUDED )
#define  KRATOS_MPI_DEM_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "bins_dynamic_objects_mpi.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/dem_search.h"

namespace Kratos
{

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
class MPI_DEMSearch : public DEMSearch<MPI_DEMSearch>
{
    public:
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of MPI_DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(MPI_DEMSearch);
      
      typedef MpiDiscreteParticleConfigure<3>   Configure;
      typedef BinsObjectDynamicMpi<Configure>   BinsType;
      
      typedef Configure::IteratorType           IteratorType;
    
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      MPI_DEMSearch(){}

      /// Destructor.
      ~MPI_DEMSearch(){}

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
      void SearchElementsInRadiusExclusiveImplementation (
          ModelPart& rModelPart,
          ElementsContainerType rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          Clean_Modelpart(rModelPart);
        
          // Get the data
          ElementsContainerType::ContainerType& elements_array      = rElements.GetContainer();
          ElementsContainerType::ContainerType& elements_ModelPart  = rModelPart.GetCommunicator().LocalMesh().ElementsArray();
         
          int NumberOfSearchElements = elements_array.size();
          int NumberOfModelPElements = elements_ModelPart.size();
                   
          // Generate the bins
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          // Perform the search
          std::vector<std::size_t> NumberOfResults(NumberOfModelPElements);
          
          for (IteratorType particle_pointer_it = elements_array.begin();
                particle_pointer_it != elements_array.end(); ++particle_pointer_it)
          {                   
              rResults[particle_pointer_it-elements_array.begin()].resize(NumberOfModelPElements);
              rResultsDistance[particle_pointer_it-elements_array.begin()].resize(NumberOfModelPElements);
          }
          
          bins.SearchObjectsMpi(rModelPart,elements_array,NumberOfSearchElements,Radius,rResults,rResultsDistance,NumberOfResults,NumberOfModelPElements,rModelPart.pGetCommunicator());
       
          // Update the modelpart interface and keep the coherence between domains
          int ResultCounter = 0;
          
          for (IteratorType particle_pointer_it = elements_array.begin();
               particle_pointer_it != elements_array.end(); ++particle_pointer_it, ++ResultCounter)
          {                   
              unsigned int neighbour_counter = 0;

              for (ResultIteratorType neighbour_it = rResults[ResultCounter].begin(); neighbour_counter < NumberOfResults[ResultCounter]; ++neighbour_it, ++neighbour_counter)
              {         
                  Add_To_Modelpart(rModelPart,neighbour_it);
              }
          }
          
          // Finally sort model for correct sync
          Sort_Modelpart(rModelPart);
          
          KRATOS_CATCH(" ")
      }
      
      void SearchNodesInRadiusExclusiveImplementation (
          ModelPart& rModelPart,
          NodesContainerType rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
      }
        
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "MPIDemSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MPIDemSearch";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      virtual void Clean_Modelpart(ModelPart& r_model_part)
      {
          KRATOS_TRY
          
          Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();

          unsigned int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();

          ModelPart::ElementsContainerType    ETempGhost[NumberOfRanks];
          ModelPart::ElementsContainerType    ETempLocal[NumberOfRanks];
          ModelPart::NodesContainerType       NTempGhost[NumberOfRanks];
          ModelPart::NodesContainerType       NTempLocal[NumberOfRanks];

          //Clean the ghost(i) and local(i) meshes

          for(unsigned int i = 0; i < NumberOfRanks; i++)
          {
              ETempGhost[i].swap(r_model_part.GetCommunicator().GhostMesh(i).Elements());
              ETempLocal[i].swap(r_model_part.GetCommunicator().LocalMesh(i).Elements());
              NTempGhost[i].swap(r_model_part.GetCommunicator().GhostMesh(i).Nodes());
              NTempLocal[i].swap(r_model_part.GetCommunicator().LocalMesh(i).Nodes());
          }

          //Celan the ghost mesh

          ModelPart::ElementsContainerType  ETempGhostGlobal;
          ModelPart::NodesContainerType     NTempGhostGlobal;

          ETempGhostGlobal.swap(r_model_part.GetCommunicator().GhostMesh().Elements());
          NTempGhostGlobal.swap(r_model_part.GetCommunicator().GhostMesh().Nodes());
          
          KRATOS_CATCH(" ")
      }
        
      //TODO: Enable Local nodes again and remove them from the search function
      void Add_To_Modelpart(ModelPart& r_model_part, ResultIteratorType neighbour_it)
      {
          KRATOS_TRY
          
          #pragma omp critical
          {
              Communicator::NeighbourIndicesContainerType communicator_ranks = r_model_part.GetCommunicator().NeighbourIndices();
              
              ElementsContainerType::ContainerType& pGhostElements = r_model_part.GetCommunicator().GhostMesh().ElementsArray();
              
              int NumberOfRanks = r_model_part.GetCommunicator().GetNumberOfColors();
              int destination = -1;
              
              bool IsInGhostMesh = false;
              //bool IsInLocalMesh = false;
            
              for(int i = 0; i < NumberOfRanks; i++)
                  if((*neighbour_it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) == communicator_ranks[i])
                      destination = i;
                              
              if(destination > -1)
              {   
                  for(IteratorType element_it = pGhostElements.begin(); !IsInGhostMesh && element_it != pGhostElements.end(); ++element_it)
                      if((*element_it)->GetGeometry()(0)->Id() == (*neighbour_it)->GetGeometry()(0)->Id())
                          IsInGhostMesh = true;
                  
                  /*
                  for(IteratorType element_it = pLocalElements.begin(); !IsInLocalMesh && element_it != pLocalElements.end(); ++element_it)
                      if((*element_it)->GetGeometry()(0)->Id() == (*neighbour_it)->GetGeometry()(0)->Id())
                          IsInLocalMesh = true;
                */
                          
                  if(!IsInGhostMesh /*&& !IsInLocalMesh*/)
                  {
                      r_model_part.GetCommunicator().GhostMesh().Elements().push_back((*neighbour_it));
                      r_model_part.GetCommunicator().GhostMesh().Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                  }
                  
                  IsInGhostMesh = false;
                  //IsInLocalMesh = false;
                
                  ElementsContainerType::ContainerType& pMyGhostElements = r_model_part.GetCommunicator().GhostMesh(destination).ElementsArray();
                  //ContainerType& pMyLocalElements = r_model_part.GetCommunicator().LocalMesh(destination).ElementsArray();
            
                  for(IteratorType element_it = pMyGhostElements.begin(); !IsInGhostMesh && element_it != pMyGhostElements.end(); ++element_it)
                      if((*element_it)->GetGeometry()(0)->Id() == (*neighbour_it)->GetGeometry()(0)->Id())
                          IsInGhostMesh = true;
                  
                  /*
                  for(IteratorType element_it = pMyLocalElements.begin(); !IsInLocalMesh && element_it != pMyLocalElements.end(); ++element_it)
                      if((*element_it)->GetGeometry()(0)->Id() == (*particle_pointer_it)->GetGeometry()(0)->Id())
                          IsInLocalMesh = true;
                  */
                  
                  if(!IsInGhostMesh)
                  {   
                      r_model_part.GetCommunicator().GhostMesh(destination).Elements().push_back((*neighbour_it));
                      r_model_part.GetCommunicator().GhostMesh(destination).Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                  }
                  
                  /*
                  if(!IsInLocalMesh)
                  {
                      r_model_part.GetCommunicator().LocalMesh(destination).Elements().push_back((*particle_pointer_it));
                      r_model_part.GetCommunicator().LocalMesh(destination).Nodes().push_back((*particle_pointer_it)->GetGeometry()(0));
                  }
                  */
              }
          }
          
          KRATOS_CATCH(" ")
      }
             
      void Sort_Modelpart(ModelPart& r_model_part)
      {
          KRATOS_TRY
          
          for (unsigned int i = 0; i < r_model_part.GetCommunicator().LocalMeshes().size(); i++)
              r_model_part.GetCommunicator().LocalMesh(i).Nodes().Unique();
          
          for (unsigned int i = 0; i < r_model_part.GetCommunicator().GhostMeshes().size(); i++)
              r_model_part.GetCommunicator().GhostMesh(i).Nodes().Unique();
          
          KRATOS_CATCH(" ")
      }
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      MPI_DEMSearch& operator=(MPI_DEMSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      MPI_DEMSearch(MPI_DEMSearch const& rOther)
      {
          *this = rOther;
      }

        
      ///@}    
        
    }; // Class DEMSearch

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
                    MPI_DEMSearch& rThis){return rIStream;}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const MPI_DEMSearch& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
    
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_MPI_DEM_SEARCH_H_INCLUDED  defined 


