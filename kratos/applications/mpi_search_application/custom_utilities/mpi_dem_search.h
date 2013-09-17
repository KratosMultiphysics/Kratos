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
      MPI_DEMSearch(Communicator& comm) : mCommunicator(comm){
      }

      /// Destructor.
      ~MPI_DEMSearch(){}

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
      void SearchElementsInRadiusExclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          //TODO: Temporal compile fix. Commincator should be a member class from now on.
          ModelPart rModelPart;
          Clean_Modelpart(rModelPart);

          // Get the data
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
         
          int NumberOfSearchElements = elements_array.size();
          int NumberOfModelPElements = elements_ModelPart.size();
                   
          // Generate the bins
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          // Perform the search
          std::vector<std::size_t> NumberOfResults(NumberOfModelPElements);
          
          for (IteratorType particle_pointer_it = elements_array.begin();
                particle_pointer_it != elements_array.end(); ++particle_pointer_it)
          {                   
              rResults[particle_pointer_it-elements_array.begin()].resize(20);
              rResultsDistance[particle_pointer_it-elements_array.begin()].resize(20);
          }
          
          bins.SearchObjectsMpi(elements_array,NumberOfSearchElements,Radius,rResults,rResultsDistance,NumberOfResults,20,mCommunicator);

          for(int i = 0; i < NumberOfSearchElements; i++)
          {
              rResults[i].resize(NumberOfResults[i]);
              rResultsDistance[i].resize(NumberOfResults[i]);
          }

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
      
     void SearchElementsInRadiusInclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
      
      }

      void SearchElementsInRadiusExclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults )
      {     

      }
      
      void SearchElementsInRadiusInclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults )
      {     
    
      }
        
        
      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     

      }
      
      void SearchNodesInRadiusInclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     

      }
      
      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults )
      {     

      }
      
      void SearchNodesInRadiusInclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults )
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
            
      Communicator& mCommunicator;
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
      virtual void Clean_Modelpart(ModelPart& r_model_part)
      {
          KRATOS_TRY
          
          Communicator::NeighbourIndicesContainerType communicator_ranks = mCommunicator.NeighbourIndices();

          unsigned int NumberOfRanks = mCommunicator.GetNumberOfColors();

          ModelPart::ElementsContainerType    ETempGhost[NumberOfRanks];
          ModelPart::ElementsContainerType    ETempLocal[NumberOfRanks];
          ModelPart::NodesContainerType       NTempGhost[NumberOfRanks];
          ModelPart::NodesContainerType       NTempLocal[NumberOfRanks];

          //Clean the ghost(i) and local(i) meshes

          for(unsigned int i = 0; i < NumberOfRanks; i++)
          {
              ETempGhost[i].swap(mCommunicator.GhostMesh(i).Elements());
              ETempLocal[i].swap(mCommunicator.LocalMesh(i).Elements());
              NTempGhost[i].swap(mCommunicator.GhostMesh(i).Nodes());
              NTempLocal[i].swap(mCommunicator.LocalMesh(i).Nodes());
          }

          //Celan the ghost mesh

          ModelPart::ElementsContainerType  ETempGhostGlobal;
          ModelPart::NodesContainerType     NTempGhostGlobal;

          ETempGhostGlobal.swap(mCommunicator.GhostMesh().Elements());
          NTempGhostGlobal.swap(mCommunicator.GhostMesh().Nodes());
          
          KRATOS_CATCH(" ")
      }
        
      //TODO: Enable Local nodes again and remove them from the search function
      void Add_To_Modelpart(ModelPart& r_model_part, ResultIteratorType neighbour_it)
      {
          KRATOS_TRY
          
          #pragma omp critical
          {
              Communicator::NeighbourIndicesContainerType communicator_ranks = mCommunicator.NeighbourIndices();
              
              ElementsContainerType::ContainerType& pGhostElements = mCommunicator.GhostMesh().ElementsArray();
              
              int NumberOfRanks = mCommunicator.GetNumberOfColors();
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
                      mCommunicator.GhostMesh().Elements().push_back((*neighbour_it));
                      mCommunicator.GhostMesh().Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                  }
                  
                  IsInGhostMesh = false;
                  //IsInLocalMesh = false;
                
                  ElementsContainerType::ContainerType& pMyGhostElements = mCommunicator.GhostMesh(destination).ElementsArray();
                  //ContainerType& pMyLocalElements = mCommunicator.LocalMesh(destination).ElementsArray();
            
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
                      mCommunicator.GhostMesh(destination).Elements().push_back((*neighbour_it));
                      mCommunicator.GhostMesh(destination).Nodes().push_back((*neighbour_it)->GetGeometry()(0));
                  }
                  
                  /*
                  if(!IsInLocalMesh)
                  {
                      mCommunicator.LocalMesh(destination).Elements().push_back((*particle_pointer_it));
                      mCommunicator.LocalMesh(destination).Nodes().push_back((*particle_pointer_it)->GetGeometry()(0));
                  }
                  */
              }
          }
          
          KRATOS_CATCH(" ")
      }
             
      void Sort_Modelpart(ModelPart& r_model_part)
      {
          KRATOS_TRY
          
          for (unsigned int i = 0; i < mCommunicator.LocalMeshes().size(); i++)
              mCommunicator.LocalMesh(i).Nodes().Unique();
          
          for (unsigned int i = 0; i < mCommunicator.GhostMeshes().size(); i++)
              mCommunicator.GhostMesh(i).Nodes().Unique();
          
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
      MPI_DEMSearch(MPI_DEMSearch const& rOther) : mCommunicator(rOther.mCommunicator)
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


