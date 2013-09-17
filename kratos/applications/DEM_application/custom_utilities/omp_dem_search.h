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
#include "node_configure.h"
#include "geometry_configure.h"

// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"

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
      
      typedef RadiusPoint<Dimension>                    PointType;
      typedef PointType*                                PtrPointType;
      typedef std::vector<PtrPointType>*                PointVector;
      typedef std::vector<PtrPointType>::iterator       PointIterator;
      
      typedef double*                                   DistanceVector;
      typedef double*                                   DistanceIterator;
      
      //Configure Types
      typedef DiscreteParticleConfigure<3>              ElementConfigureType;   //Element
      typedef NodeConfigure<3>                          NodeConfigureType;      //Node
      typedef GeometryConfigure<3>                      GeometryConfigureType;  //Geometry
      
      //Bin Types
      typedef BinsObjectDynamic<ElementConfigureType>   BinsType;
      typedef BinsObjectDynamic<NodeConfigureType>      NodeBinsType;
      typedef BinsObjectDynamic<GeometryConfigureType>  GemoetryBinsType; 
      
      typedef BinsDynamic<3,PointType,PointVector,PtrPointType,PointIterator,DistanceIterator,PointDistance2<Dimension,PointType> >    PointBins;
     
      
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      OMP_DEMSearch(){
//           searchPoints = new std::vector<PtrPointType>(0);
      }

      /// Destructor.
      ~OMP_DEMSearch(){
//           delete searchPoints;
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
              
              #pragma omp for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(elements_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
                  
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
              
              #pragma omp for
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
      
      //ELEMENT SEARCH BASED
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
        
          BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());
          
          #pragma omp parallel
          {
              ResultElementsContainerType   localResults(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for
              for(std::size_t i = 0; i < elements_array.size(); i++)
              {
                  ResultElementsContainerType::iterator ResultsPointer = localResults.begin();
                        
                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(elements_array[i],Radius[i],ResultsPointer,MaxNumberOfElements);
  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);    
              } 
          }
          
          KRATOS_CATCH("")        
      }

//       //GEOMETRY SEARCH BASED
//       void SearchElementsInRadiusExclusiveImplementation (
//           ElementsContainerType const& rStructureElements,
//           ElementsContainerType const& rElements,
//           const RadiusArrayType & Radius, 
//           VectorResultElementsContainerType& rResults )
//       {     
//           KRATOS_TRY
//           
//           int MaxNumberOfElements = rStructureElements.size();
// 
//           searchPoints->resize(MaxNumberOfElements);
//           
//           ElementsContainerType::ContainerType& elements_array = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
//           
//           PointIterator PointsBegin = searchPoints->begin();
//           PointIterator PointsEnd   = searchPoints->end();
//            
//           static int a = 1;
//           if(a)
//           {
//               for(ElementsContainerType::ContainerType::iterator i_elem = elements_array.begin(); i_elem != elements_array.end(); i_elem++)
//               {
//                   std::size_t index = i_elem - elements_array.begin();
//                   
//                   (*searchPoints)[index] = new PointType();
//               }
//               a = 0;
//           }
//             
//           for(ElementsContainerType::ContainerType::iterator i_elem = elements_array.begin(); i_elem != elements_array.end(); i_elem++)
//           {
//               std::size_t index = i_elem - elements_array.begin();
// 
//               (*searchPoints)[index]->Initialize((*i_elem));
//           }
// 
//           PointBins bins(PointsBegin, PointsEnd);
//           
//           #pragma omp parallel
//           {
//               PointVector localResults = new std::vector<PtrPointType>(MaxNumberOfElements);
//               double distances[MaxNumberOfElements];
//               std::size_t                   NumberOfResults = 0;
//               
//               #pragma omp for
//               for(std::size_t i = 0; i < searchPoints->size(); i++)
//               {
//                   std::vector<PtrPointType>::iterator ResultsPointer = localResults->begin();
//                 
//                   NumberOfResults = bins.SearchInRadius((*(*searchPoints)[i]),Radius[i],ResultsPointer,distances,MaxNumberOfElements);
// 
//                   rResults[i].resize(NumberOfResults-1);
//                   int k = 0;
//                   
//                   for(std::size_t j = 0; j < NumberOfResults; j++)
//                       if((*localResults)[j]->pNaseElem->Id() != (*searchPoints)[i]->pNaseElem->Id())
//                           rResults[i][k++] = boost::make_shared<ElementType>(*(*localResults)[j]->pNaseElem);
//               }
//           }
// 
//           KRATOS_CATCH("")
//       }

//       //POINT SEARCH BASED
//       void SearchElementsInRadiusExclusiveImplementation (
//           ElementsContainerType const& rStructureElements,
//           ElementsContainerType const& rElements,
//           const RadiusArrayType & Radius, 
//           VectorResultElementsContainerType& rResults )
//       {     
//           KRATOS_TRY
//           
//           int MaxNumberOfElements = rStructureElements.size();
//           
//           ElementsContainerType::ContainerType& elements_array = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());
//           
// //           KRATOS_TIMER_START("OMP_DEM_ALLOCATE_MEMORY")
//           if(searchPoints->size() < MaxNumberOfElements) 
//           {
//               //Keep allocating data chunks until the right size is achieved
//               int old_size = searchPoints->size();
//               int cur_size = searchPoints->size();
//               int fnl_size = MaxNumberOfElements;
//               int chk_size = 1024;                                          //Charlie: To be defined by the user?
//               
//               std::cout << "omp_dem_search.h: Not enought space in the result buffer. Allocating memory" << std::endl;
//               
//               while(cur_size < fnl_size) cur_size += chk_size;      //Charlie: It can be calculated rather than do it in a while...
//               
//               std::cout << "\told buffer size: " << old_size << std::endl;
//               std::cout << "\tnew buffer size: " << cur_size << std::endl;
//               
//               searchPoints->resize(MaxNumberOfElements);
//               for(ElementsContainerType::ContainerType::iterator i_elem = elements_array.begin() + old_size; i_elem != elements_array.end(); i_elem++)
//               {
//                   std::size_t index = i_elem - elements_array.begin();
//                   
//                   (*searchPoints)[index] = new PointType();
//               }
//           }
// //           KRATOS_TIMER_STOP("OMP_DEM_ALLOCATE_MEMORY")
//                    
//           PointIterator PointsBegin = searchPoints->begin();
//           PointIterator PointsEnd   = searchPoints->end();
//             
// //           KRATOS_TIMER_START("OMP_DEM_INITIALIZE")
//           for(ElementsContainerType::ContainerType::iterator i_elem = elements_array.begin(); i_elem != elements_array.end(); i_elem++)
//           {
//               std::size_t index = i_elem - elements_array.begin();
// 
//               (*searchPoints)[index]->Initialize((*i_elem));
//           }
// //           KRATOS_TIMER_STOP("OMP_DEM_INITIALIZE")
// 
// //           KRATOS_TIMER_START("OMP_DEM_GENERATE_BINS")
//           PointBins bins(PointsBegin, PointsEnd);
// //           KRATOS_TIMER_STOP("OMP_DEM_GENERATE_BINS")
//           
//           #pragma omp parallel
//           {
//               PointVector localResults = new std::vector<PtrPointType>(MaxNumberOfElements);
//               double distances[MaxNumberOfElements];
//               std::size_t                   NumberOfResults = 0;
//               
//               #pragma omp for
//               for(std::size_t i = 0; i < searchPoints->size(); i++)
//               {
//                   std::vector<PtrPointType>::iterator ResultsPointer = localResults->begin();
//                 
// //                   KRATOS_TIMER_START("OMP_DEM_SEARCH")
//                   NumberOfResults = bins.SearchInRadius((*(*searchPoints)[i]),Radius[i],ResultsPointer,distances,MaxNumberOfElements);
// //                   KRATOS_TIMER_STOP("OMP_DEM_SEARCH")
// 
//                   rResults[i].resize(NumberOfResults-1);
//                   int k = 0;
//                   
// //                   KRATOS_TIMER_START("OMP_DEM_CPY_RES")
//                   for(std::size_t j = 0; j < NumberOfResults; j++)
//                       if((*localResults)[j]->pNaseElem->Id() != (*searchPoints)[i]->pNaseElem->Id())
//                           rResults[i][k++] = boost::make_shared<ElementType>(*(*localResults)[j]->pNaseElem);
// //                   KRATOS_TIMER_STOP("OMP_DEM_CPY_RES")
//               }
//               
//               delete localResults;
//           }
// 
//           KRATOS_CATCH("")
//       }
      
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
                        
                  NumberOfResults = bins.SearchObjectsInRadius(elements_array[i],Radius[i],ResultsPointer,MaxNumberOfElements);
                          
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


