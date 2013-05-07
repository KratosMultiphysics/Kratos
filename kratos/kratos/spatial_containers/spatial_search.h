//
//   Project Name:        Kratos
//   Last Modified by:    $Author: croig $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
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
      typedef Point<Dimension, double>                          PointType;
     
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
      //************************************************************************
     
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusExclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().ElementsArray(), 
                                                Radius,rResults,rResultsDistance);
          
      }
     
      virtual void SearchElementsInRadiusExclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          this->SearchElementsInRadiusInclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().ElementsArray(), 
                                                Radius,rResults,rResultsDistance);
      }
      
      virtual void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      virtual void SearchNodesInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusExclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().NodesArray(), 
                                             Radius,rResults,rResultsDistance);
      }
      
      virtual void SearchNodesInRadiusExclusive (
          ModelPart& rModelPart,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      virtual void SearchNodesInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchNodesInRadiusInclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().NodesArray(), 
                                             Radius,rResults,rResultsDistance);
      }
      
      virtual void SearchNodesInRadiusInclusive (
          ModelPart& rModelPart,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      virtual void SearchConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchConditionsInRadiusExclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().ConditionsArray(), 
                                                  Radius,rResults,rResultsDistance);
      }
      
      virtual void SearchConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
      
      virtual void SearchConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance ) 
      {
          this->SearchConditionsInRadiusInclusive(rModelPart, rModelPart.GetCommunicator().LocalMesh().ConditionsArray(), 
                                                  Radius,rResults,rResultsDistance);
      }
      
      virtual void SearchConditionsInRadiusInclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& InputConditions,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          /* Abstract */
          KRATOS_ERROR(std::runtime_error,"Direct call of an abstract method","")
      }
      
      //************************************************************************
      //************************************************************************
        
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


