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
#include "spatial_containers/dem_search.h"
#include "utilities/openmp_utils.h"

// Configures
#include "discrete_particle_configure.h"
#include "geometrical_object_configure.h"
#include "node_configure.h"

// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"
#include "custom_search/bins_dynamic_objects_periodic.h"

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
    
      /// Pointer definition of OMP_DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(OMP_DEMSearch);
      
      typedef PointType*                                    PtrPointType;
      typedef std::vector<PtrPointType>*                    PointVector;
      typedef std::vector<PtrPointType>::iterator           PointIterator;
      
      typedef double*                                       DistanceVector;
      typedef double*                                       DistanceIterator;
      
      //Configure Types
      typedef DiscreteParticleConfigure<3>                  ElementConfigureType;       //Element
      typedef NodeConfigure<3>                              NodeConfigureType;          //Node
      typedef GeometricalConfigure<3>                       GeometricalConfigureType;   //Generic Geometry
      
      //Bin Types
      typedef BinsObjectDynamic<ElementConfigureType>               BinsType;
      typedef BinsObjectDynamicPeriodic<ElementConfigureType>       BinsTypePeriodic;
      typedef std::unique_ptr<BinsType>                             BinsUniquePointerType;
      typedef BinsObjectDynamic<NodeConfigureType>                  NodeBinsType;
      typedef BinsObjectDynamic<GeometricalConfigureType>           GeometricalBinsType;
      
      //GeoimetricalObject
      typedef PointerVectorSet<GeometricalObject, IndexedObject>     GeometricalObjectType;
     
      
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.

      OMP_DEMSearch(const double domain_min_x = 0.0, const double domain_min_y = 0.0, const double domain_min_z = 0.0,
                    const double domain_max_x = -1.0, const double domain_max_y = -1.0, const double domain_max_z = -1.0)
      {
            mDomainPeriodicity = (domain_min_x <= domain_max_x) ? true : false;
      }

      /// Destructor.
      ~OMP_DEMSearch(){
      }
      

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
//           KRATOS_TRY
//           
//           int MaxNumberOfElements = rStructureElements.size();
//           
//           ElementsContainerType::ContainerType& elements_bins = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
//           ElementsContainerType::ContainerType& elements_sear = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
//           
//           GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
//           GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;
// 
//           BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());
//           SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
// 
//           for(ElementsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
//               BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);          
//           
//           for(ElementsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
//               SearElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
//           
//           GeometricalBinsType bins(BinsElementPointerToGeometricalObjecPointerTemporalVector.begin(), BinsElementPointerToGeometricalObjecPointerTemporalVector.end());
//           
//           #pragma omp parallel
//           {
//               GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
//               DistanceType                          localResultsDistances(MaxNumberOfElements);
//               std::size_t                           NumberOfResults = 0;
//               
//               #pragma omp for
//               for(std::size_t i = 0; i < elements_sear.size(); i++)
//               {
//                   GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
//                   DistanceType::iterator                                                        ResultsDistancesPointer = localResultsDistances.begin();
// 
//                   NumberOfResults = bins.SearchObjectsInRadiusExclusive(SearElementPointerToGeometricalObjecPointerTemporalVector[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
// 
//                   rResults[i].reserve(NumberOfResults);
//                   
//                   for(GeometricalObjectType::ContainerType::iterator it = localResults.begin(); it != localResults.begin() + NumberOfResults; it++)
//                   {
//                       Element::Pointer elem = Kratos::dynamic_pointer_cast<Element>(*it);
//                       rResults[i].push_back(elem);
//                       rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
//                   }
//               }
//           }
//           
//           KRATOS_CATCH("")

          KRATOS_TRY

          int MaxNumberOfElements = rStructureElements.size();
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
          BinsUniquePointerType p_bins = GetBins(elements_ModelPart);

          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              DistanceType                  localResultsDistances(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for schedule(dynamic, 100) //schedule(guided)
              for(int i = 0; i < static_cast<int>(elements_array.size()); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                  
                  SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*elements_array[i]);
                  const double radius = p_particle->GetSearchRadius();

                  NumberOfResults = p_bins->SearchObjectsInRadiusExclusive(elements_array[i],radius,ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);      
              }
          }         
          //MAJOR TODO: creating and destroying (when leaving the function) this BINS is not parallel and takes a significant time if we search at every time step. Can we re-use a bins and avoid allocation and deallocation?? MA
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
          BinsUniquePointerType p_bins = GetBins(elements_ModelPart);

          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              DistanceType                  localResultsDistances(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_array.size()); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                  
                  SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*elements_array[i]);
                  const double radius = p_particle->GetSearchRadius();
                
                  NumberOfResults = p_bins->SearchObjectsInRadius(elements_array[i],radius,ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
                  
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

          int MaxNumberOfElements = rStructureElements.size();
          
          ElementsContainerType::ContainerType& elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
          ElementsContainerType::ContainerType& elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
          BinsUniquePointerType p_bins = GetBins(elements_ModelPart);

          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_array.size()); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer = localResults.begin();
                  
                  SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*elements_array[i]);
                  const double radius = p_particle->GetSearchRadius();

                  NumberOfResults = p_bins->SearchObjectsInRadiusExclusive(elements_array[i],radius,ResultsPointer,MaxNumberOfElements);

                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
              } 
          }
          
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
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_array.size()); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer = localResults.begin();
                  
                  SphericParticle* p_particle = dynamic_cast<SphericParticle*>(&*elements_array[i]);
                  const double radius = p_particle->GetSearchRadius();
                        
                  NumberOfResults = bins.SearchObjectsInRadius(elements_array[i],radius,ResultsPointer,MaxNumberOfElements);
                          
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
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator    ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(nodes_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfNodes);
                  
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
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
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
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(nodes_array[i],Radius[i],ResultsPointer,MaxNumberOfNodes);
                  
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
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
                        
                  NumberOfResults = bins.SearchObjectsInRadius(nodes_array[i],Radius[i],ResultsPointer,MaxNumberOfNodes);
                          
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);    
              }
          }
          
          KRATOS_CATCH("")
      }

      void SearchGeometricalInRadiusExclusiveImplementation (
          ElementsContainerType   const& rStructureElements,
          ConditionsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultConditionsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY

          int MaxNumberOfElements = rStructureElements.size();

          ElementsContainerType::ContainerType& elements_bins   = const_cast<ElementsContainerType::ContainerType&>  (rStructureElements.GetContainer());
          ConditionsContainerType::ContainerType& elements_sear = const_cast<ConditionsContainerType::ContainerType&>(rElements.GetContainer());

          GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
          GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;

          SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
          BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());
              
          for(ElementsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
              BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          for(ConditionsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
              SearElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          GeometricalBinsType bins(BinsElementPointerToGeometricalObjecPointerTemporalVector.begin(), BinsElementPointerToGeometricalObjecPointerTemporalVector.end());
          
          #pragma omp parallel
          {
              GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
              DistanceType                          localResultsDistances(MaxNumberOfElements);
              std::size_t                           NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_sear.size()); i++)
              {
                  GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
                  DistanceType::iterator                                                        ResultsDistancesPointer = localResultsDistances.begin();

                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(SearElementPointerToGeometricalObjecPointerTemporalVector[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

                  rResults[i].reserve(NumberOfResults);
                  
                  for(GeometricalObjectType::ContainerType::iterator it = localResults.begin(); it != localResults.begin() + NumberOfResults; it++)
                  {
                      Condition::Pointer elem = Kratos::dynamic_pointer_cast<Condition>(*it);
                      rResults[i].push_back(elem);
                      rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
                  }
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchGeometricalInRadiusInclusiveImplementation (
          ElementsContainerType   const& rStructureElements,
          ConditionsContainerType const& rElements,
          const RadiusArrayType& Radius, 
          VectorResultConditionsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY

          int MaxNumberOfElements = rStructureElements.size();

          ElementsContainerType::ContainerType& elements_bins   = const_cast<ElementsContainerType::ContainerType&>  (rStructureElements.GetContainer());
          ConditionsContainerType::ContainerType& elements_sear = const_cast<ConditionsContainerType::ContainerType&>(rElements.GetContainer());

          GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
          GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;

          SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
          BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());
              
          for(ElementsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
              BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          for(ConditionsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
              SearElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          GeometricalBinsType bins(BinsElementPointerToGeometricalObjecPointerTemporalVector.begin(), BinsElementPointerToGeometricalObjecPointerTemporalVector.end());
          
          #pragma omp parallel
          {
              GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
              DistanceType                          localResultsDistances(MaxNumberOfElements);
              std::size_t                           NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_sear.size()); i++)
              {
                  GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
                  DistanceType::iterator                                                        ResultsDistancesPointer = localResultsDistances.begin();

                  NumberOfResults = bins.SearchObjectsInRadius(SearElementPointerToGeometricalObjecPointerTemporalVector[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

                  rResults[i].reserve(NumberOfResults);
                  
                  for(GeometricalObjectType::ContainerType::iterator it = localResults.begin(); it != localResults.begin() + NumberOfResults; it++)
                  {
                      Condition::Pointer elem = Kratos::dynamic_pointer_cast<Condition>(*it);
                      rResults[i].push_back(elem);
                      rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
                  }
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchGeometricalInRadiusExclusiveImplementation (
          ConditionsContainerType const& rStructureElements,
          ElementsContainerType   const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {
          KRATOS_TRY

          int MaxNumberOfElements = rStructureElements.size();

          ConditionsContainerType::ContainerType& elements_bins = const_cast<ConditionsContainerType::ContainerType&>(rStructureElements.GetContainer());
          ElementsContainerType::ContainerType& elements_sear   = const_cast<ElementsContainerType::ContainerType&>  (rElements.GetContainer());

          GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
          GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;

          SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
          BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());

          for(ElementsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
              SearElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
              
          for(ConditionsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
              BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          GeometricalBinsType bins(BinsElementPointerToGeometricalObjecPointerTemporalVector.begin(), BinsElementPointerToGeometricalObjecPointerTemporalVector.end());
          
          #pragma omp parallel
          {
              GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
              DistanceType                          localResultsDistances(MaxNumberOfElements);
              std::size_t                           NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_sear.size()); i++)
              {
                  GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
                  DistanceType::iterator                                                        ResultsDistancesPointer = localResultsDistances.begin();

                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(SearElementPointerToGeometricalObjecPointerTemporalVector[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

                  rResults[i].reserve(NumberOfResults);
                  
                  for(GeometricalObjectType::ContainerType::iterator it = localResults.begin(); it != localResults.begin() + NumberOfResults; it++)
                  {
                      Element::Pointer elem = Kratos::dynamic_pointer_cast<Element>(*it);
                      rResults[i].push_back(elem);
                      rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
                  }
              }
          }
          
          KRATOS_CATCH("")
      }
      
      void SearchGeometricalInRadiusInclusiveImplementation (
          ConditionsContainerType const& rStructureElements,
          ElementsContainerType   const& rElements,
          const RadiusArrayType& Radius, 
          VectorResultElementsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          int MaxNumberOfElements = rStructureElements.size();

          ConditionsContainerType::ContainerType& elements_bins = const_cast<ConditionsContainerType::ContainerType&>(rStructureElements.GetContainer());
          ElementsContainerType::ContainerType& elements_sear   = const_cast<ElementsContainerType::ContainerType&>  (rElements.GetContainer());

          GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
          GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;

          SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
          BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());

          for(ElementsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
              SearElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
              
          for(ConditionsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
              BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          GeometricalBinsType bins(BinsElementPointerToGeometricalObjecPointerTemporalVector.begin(), BinsElementPointerToGeometricalObjecPointerTemporalVector.end());
          
          #pragma omp parallel
          {
              GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
              DistanceType                          localResultsDistances(MaxNumberOfElements);
              std::size_t                           NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_sear.size()); i++)
              {
                  GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
                  DistanceType::iterator                                                        ResultsDistancesPointer = localResultsDistances.begin();

                  NumberOfResults = bins.SearchObjectsInRadius(SearElementPointerToGeometricalObjecPointerTemporalVector[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

                  rResults[i].reserve(NumberOfResults);
                  
                  for(GeometricalObjectType::ContainerType::iterator it = localResults.begin(); it != localResults.begin() + NumberOfResults; it++)
                  {
                      Element::Pointer elem = Kratos::dynamic_pointer_cast<Element>(*it);
                      rResults[i].push_back(elem);
                      rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
                  }
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
      virtual std::string Info() const override
      {
          std::stringstream buffer;
          buffer << "OpenMPDemSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "OpenMPDemSearch";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override {}
      
            
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
      ///
        BinsUniquePointerType GetBins(ElementsContainerType::ContainerType& r_model_part_container)
        {
            if (mDomainPeriodicity){
                return std::unique_ptr<BinsType>(new BinsTypePeriodic(r_model_part_container.begin(), r_model_part_container.end(), this->mDomainMin, this->mDomainMax));
            }

            else {
                return std::unique_ptr<BinsType>(new BinsType(r_model_part_container.begin(), r_model_part_container.end()));
            }
        }
        
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


