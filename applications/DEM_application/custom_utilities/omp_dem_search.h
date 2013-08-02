//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_OMP_DEM_SEARCH_H_INCLUDED )
#define  KRATOS_OMP_DEM_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/dem_search.h"
#include "discrete_particle_configure.h"
#include "node_configure.h"
#include "utilities/openmp_utils.h"

// External includes

/* Timer defines */
#include "utilities/timer.h"
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

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

class OMP_DEMSearch : public DEMSearch<OMP_DEMSearch>
{
    public:
      ///@name Type Definitions
      ///@{
        
      enum SearchType { 
        ELEM2ELEM = 1,
        ELEM2NODE = 2,
        ELEM2COND = 4,
        NODE2ELEM = 8,
        NODE2NODE = 16,
        NODE2COND = 32,
        COND2ELEM = 64,
        COND2NODE = 128,
        COND2COND = 256
      };
      
      enum SearchMode {
        INCLUSIVE = 1,
        EXCLUSIVE = 2
      };
    
      /// Pointer definition of OMP_DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(OMP_DEMSearch);
      
      //Configure Types
      typedef DiscreteParticleConfigure<3>              ElementConfigureType;   //Element
      typedef NodeConfigure<3>                          NodeConfigureType;      //Node
      
      typedef BinsObjectDynamic<ElementConfigureType>   BinsType;
      typedef BinsObjectDynamic<NodeConfigureType>      NodeBinsType;
     
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      OMP_DEMSearch(){}

      /// Destructor.
      ~OMP_DEMSearch(){}
      

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
          
          int MaxNumberOfElements = rStructureElements.size();
          
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
        
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
                    
          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              DistanceType                  localResultsDistances(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadius(elements_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);      
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchElementsInRadiusInclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType& Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          int MaxNumberOfElements = rStructureElements.size();
          
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
        
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              DistanceType                  localResultsDistances(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusInner(elements_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);      
              }
          }
          
          KRATOS_CATCH("")        
      }

      void SearchElementsInRadiusExclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults )
      {     
          KRATOS_TRY
          KRATOS_TIMER_START("SN-LVLN1")
          
          int MaxNumberOfElements = rStructureElements.size();
          
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
        
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer = localResults.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadius(elements_array[i],Radius[i],ResultsPointer,MaxNumberOfElements);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);     
              }
          }
          
          KRATOS_TIMER_STOP("SN-LVLN1")
          KRATOS_CATCH("")
      }
      
      void SearchElementsInRadiusInclusiveImplementation (
          ElementsContainerType const& rStructureElements,
          ElementsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults )
      {     
          KRATOS_TRY
          
          int MaxNumberOfElements = rStructureElements.size();
          
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
        
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer = localResults.begin();
                        
                  NumberOfResults = bins.SearchObjectsInRadiusInner(elements_array[i],Radius[i],ResultsPointer,MaxNumberOfElements);
                          
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);    
              } 
          }
          
          KRATOS_CATCH("")        
      }

      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          int MaxNumberOfNodes = rStructureNodes.size();
          
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());
          NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
              
          NodeBinsType bins(nodes_ModelPart.begin(), nodes_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(MaxNumberOfNodes);
              DistanceType              localResultsDistances(MaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < nodes_array.size(); i++)
              {
                  ResultNodesContainerType::iterator    ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadius(nodes_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfNodes);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);      
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchNodesInRadiusInclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          int MaxNumberOfNodes = rStructureNodes.size();
      
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());
          NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
          
          NodeBinsType bins(nodes_ModelPart.begin(), nodes_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(MaxNumberOfNodes);
              DistanceType              localResultsDistances(MaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < nodes_array.size(); i++)
              {
                  ResultNodesContainerType::iterator    ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusInner(nodes_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfNodes);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);      
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults )
      {     
          KRATOS_TRY
          
          int MaxNumberOfNodes = rStructureNodes.size();
          
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());
          NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
        
          NodeBinsType bins(nodes_ModelPart.begin(), nodes_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(MaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < nodes_array.size(); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadius(nodes_array[i],Radius[i],ResultsPointer,MaxNumberOfNodes);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);     
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchNodesInRadiusInclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults )
      {     
          KRATOS_TRY
          
          int MaxNumberOfNodes = rStructureNodes.size();
          
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());
          NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
        
          NodeBinsType bins(nodes_ModelPart.begin(), nodes_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(MaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;
              
              #pragma omp parallel for
              for(std::size_t i = 0; i < nodes_array.size(); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
                        
                  NumberOfResults = bins.SearchObjectsInRadiusInner(nodes_array[i],Radius[i],ResultsPointer,MaxNumberOfNodes);
                          
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);    
              }
          }
          
          KRATOS_CATCH("")
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
          buffer << "OpenMPDemSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "OpenMPDemSearch";}

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
      OMP_DEMSearch& operator=(OMP_DEMSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      OMP_DEMSearch(OMP_DEMSearch const& rOther)
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
//   inline std::istream& operator >> (std::istream& rIStream, 
//                     DEMSearch& rThis){return rIStream;}
// 
//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
//                     const DEMSearch& rThis)
//   {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
// 
//     return rOStream;
//   }
  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DEM_SEARCH_H_INCLUDED  defined 


