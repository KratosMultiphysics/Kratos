//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//

#if !defined(KRATOS_SPATIAL_SEARCH_H_INCLUDED )
#define  KRATOS_SPATIAL_SEARCH_H_INCLUDED

// system includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <vector>

// include kratos definitions
#include "includes/define.h"

// kratos includes
#include "includes/element.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"

// kratos utils
#include "utilities/spatial_containers_configure.h"

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

class SpatialSearch 
{
    public:
      
      
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of SpatialSearch
      KRATOS_CLASS_POINTER_DEFINITION(SpatialSearch);

      enum { Dimension = 3,
             MAX_LEVEL = 16,
             MIN_LEVEL = 2
      };
      
      /// Common Defines
      typedef Point                          PointType;
     
      typedef ModelPart::ElementsContainerType                  ElementsContainerType;
      typedef ModelPart::ElementType                            ElementType;
      typedef ModelPart::ElementType::Pointer                   ElementPointerType;
      typedef ElementsContainerType::ContainerType              ResultElementsContainerType;
      typedef std::vector<ResultElementsContainerType>          VectorResultElementsContainerType;

      typedef ModelPart::NodesContainerType                     NodesContainerType;
      typedef ModelPart::NodeType                               NodeType;
      typedef ModelPart::NodeType::Pointer                      NodePointerType;
      typedef NodesContainerType::ContainerType                 ResultNodesContainerType;
      typedef std::vector<ResultNodesContainerType>             VectorResultNodesContainerType;
      
      typedef ModelPart::ConditionsContainerType                ConditionsContainerType;
      typedef ModelPart::ConditionType                          ConditionType;
      typedef ModelPart::ConditionType::Pointer                 ConditionPointerType;
      typedef ConditionsContainerType::ContainerType            ResultConditionsContainerType;
      typedef std::vector<ResultConditionsContainerType>        VectorResultConditionsContainerType;

      typedef std::vector<double>                               RadiusArrayType;
      typedef std::vector<double>                               DistanceType;
      typedef std::vector<DistanceType>                         VectorDistanceType;
      
      typedef ElementsContainerType::ContainerType::iterator    ResultIteratorType;
      
      
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      SpatialSearch(){}

      /// Destructor.
      virtual ~SpatialSearch(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
      //************************************************************************
      // Elemental Exclusive search with distance calculation                            
      //************************************************************************
    
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                InputElements, 
                                                Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every element in "StructureElements" excluding itself
       * @param StructureElements   List of elements modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusExclusive(StructureElements,
                                                StructureElements,
                                                Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Elemental Inclusive search with distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every element in "rModelpart" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every element in "InputElements" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                InputElements, 
                                                Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every element in "StructureElements" including itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusInclusive(StructureElements, 
                                                StructureElements, 
                                                Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every element in "InputElements" including itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Elemental Exclusive search without distance calculation                            
      //************************************************************************
    
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                Radius,rResults);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          this->SearchElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                InputElements, 
                                                Radius,rResults);
      }
      
      /**
       * Search neighbours for every element in "StructureElements" excluding itself
       * @param StructureElements   List of nodes against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          this->SearchElementsInRadiusExclusive(StructureElements,
                                                StructureElements,
                                                Radius,rResults);
          
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Elemental Inclusive search without distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every element in "rModelpart" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {     
          this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                Radius,rResults);
      }
      
      /**
       * Search neighbours for every element in "InputElements" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {     
          this->SearchElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                InputElements, 
                                                Radius,rResults);
      }
      
      /**
       * Search neighbours for every element in "StructureElements" including itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {     
          this->SearchElementsInRadiusInclusive(StructureElements, 
                                                StructureElements, 
                                                Radius,rResults);
      }
      
      /**
       * Search neighbours for every element in "InputElements" including itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       */
      virtual void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Nodal Exclusive search with distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every node in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), 
                                             rModelPart.GetCommunicator().LocalMesh().Nodes(), 
                                             Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every node in "InputNodes" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputNodes          List of nodes to be searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusExclusive (
          ModelPart& rModelPart,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          this->SearchNodesInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(), 
                                             InputNodes,
                                             Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every node in "StructureNodes" excluding itself
       * @param StructureNodes      Lis of nodes against which the neighbours are searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusExclusive(StructureNodes,
                                             StructureNodes,
                                             Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every node in "InputNodes" excluding itself
       * @param rModelPart          List of nodes against which the neighbours are searched
       * @param InputNodes          List of nodes to be searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Nodal Inclusive search with distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every node in "rModelpart" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                             rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                             Radius,rResults,rResultsDistance);
      }

     /**
      * Search neighbours for every node in "InputNodes" including itself
      * @param rModelPart          Input modelpart against which the neighbours are searched
      * @param InputNodes          List of nodes to be searched
      * @param Radius              List of search radius for every node
      * @param rResults            Array of results for each node
      * @param rResultDistance     Array of distances for each result of each node
      */
      virtual void SearchNodesInRadiusInclusive (
          ModelPart& rModelPart,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          this->SearchNodesInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Nodes(),
                                             InputNodes,
                                             Radius,rResults,rResultsDistance);
      }
      
      /**
       * Search neighbours for every node in "StructureNodes" including itself
       * @param StructureNodes      List of nodes against which the neighbours are searched
       * @param Radius              List of search radius for every node
       * @param rResults            Array of results for each node
       * @param rResultDistance     Array of distances for each result of each node
       */
      virtual void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusInclusive(StructureNodes, 
                                             StructureNodes, 
                                             Radius,rResults,rResultsDistance);
      }

     /**
      * Search neighbours for every node in "InputNodes" including itself
      * @param StructureNodes      List of nodes against which the neighbours are searched
      * @param InputNodes          List of nodes to be searched
      * @param Radius              List of search radius for every node
      * @param rResults            Array of results for each node
      * @param rResultDistance     Array of distances for each result of each node
      */
      virtual void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Condition Exclusive search with distance calculation                            
      //************************************************************************

      /**
       * Search neighbours for every condition in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every condition
       * @param rResults            Array of results for each condition
       * @param rResultDistance     Array of distances for each result of each condition
       */      
      virtual void SearchConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchConditionsInRadiusExclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                  Radius,rResults,rResultsDistance);
      }

     /**
      * Search neighbours for every condition in "InputConditions" excluding itself
      * @param rModelPart          Input modelpart against which the neighbours are searched
      * @param InputConditions     List of conditions to be searched
      * @param Radius              List of search radius for every condition
      * @param rResults            Array of results for each condition
      * @param rResultDistance     Array of distances for each result of each condition
      */
      virtual void SearchConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Condition Inclusive search with distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every condition in "rModelpart" including itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every condition
       * @param rResults            Array of results for each condition
       * @param rResultDistance     Array of distances for each result of each condition
       */
      virtual void SearchConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchConditionsInRadiusInclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                  Radius,rResults,rResultsDistance);
      }

     /**
      * Search neighbours for every condition in "InputConditions" including itself
      * @param rModelPart          Input modelpart against which the neighbours are searched
      * @param InputConditions     List of conditions to be searched
      * @param Radius              List of search radius for every condition
      * @param rResults            Array of results for each condition
      * @param rResultDistance     Array of distances for each result of each condition
      */
      virtual void SearchConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      //************************************************************************
      //************************************************************************
      
      //************************************************************************
      // Elemental vs Condition Exclusive search with distance calculation                            
      //************************************************************************
    
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchConditionsOverElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                              rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputConditions     List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusExclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchConditionsOverElementsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                              InputConditions, 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputConditions     List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Elemental vs Condition Inclusive search with distance calculation                            
      //************************************************************************
    
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchConditionsOverElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(),
                                                              rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputConditions     List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusInclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchConditionsOverElementsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                              InputConditions, 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputConditions     List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchConditionsOverElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Condition vs Elemental Exclusive search with distance calculation                            
      //************************************************************************
    
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsOverConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                              rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsOverConditionsInRadiusExclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                              InputElements, 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusExclusive (
          ConditionsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      // Condition vs Elemental Inclusive search with distance calculation                            
      //************************************************************************
      
      /**
       * Search neighbours for every element in "rModelpart" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsOverConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(),
                                                              rModelPart.GetCommunicator().LocalMesh().Elements(), 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param rModelPart          Input modelpart against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsOverConditionsInRadiusInclusive(rModelPart.GetCommunicator().LocalMesh().Conditions(), 
                                                              InputElements, 
                                                              Radius,rResults,rResultsDistance);
      }
     
      /**
       * Search neighbours for every element in "Inputelements" excluding itself
       * @param StructureElements   List of elements against which the neighbours are searched
       * @param InputElements       List of elements to be searched
       * @param Radius              List of search radius for every element
       * @param rResults            Array of results for each element
       * @param rResultDistance     Array of distances for each result of each element
       */
      virtual void SearchElementsOverConditionsInRadiusInclusive (
          ConditionsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_THROW_ERROR(std::runtime_error,"Direct call of an abstract method","")
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
          buffer << "SpatialSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "SpatialSearch";}

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
      SpatialSearch& operator=(SpatialSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      SpatialSearch(SpatialSearch const& rOther)
      {
          *this = rOther;
      }

        
      ///@}    
        
    }; // Class SpatialSearch

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
    
  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
                    const SpatialSearch& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_SPATIAL_SEARCH_H_INCLUDED  defined 


