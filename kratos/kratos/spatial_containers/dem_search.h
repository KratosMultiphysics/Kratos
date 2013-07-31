//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_DEM_SEARCH_H_INCLUDED )
#define  KRATOS_DEM_SEARCH_H_INCLUDED

// include kratos definitions
#include "includes/define.h"

// System includes
#include <string>
#include <iostream> 

// External includes
#include "spatial_containers/spatial_search.h"

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
template< class TDerived >
class DEMSearch : public SpatialSearch
{
    public:
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(DEMSearch);
    
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      DEMSearch(){}

      /// Destructor.
      virtual ~DEMSearch(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
      void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          static_cast<TDerived*>(this)->SearchElementsInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }
      
      void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          static_cast<TDerived*>(this)->SearchElementsInRadiusInclusiveImplementation(StructureElements,InputElements,Radius,rResults,rResultsDistance);
      }
      
      void SearchElementsInRadiusExclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          static_cast<TDerived*>(this)->SearchElementsInRadiusExclusiveImplementation(StructureElements,InputElements,Radius,rResults);
      }
      
      void SearchElementsInRadiusInclusive (
          ElementsContainerType const& StructureElements,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultElementsContainerType& rResults )
      {     
          static_cast<TDerived*>(this)->SearchElementsInRadiusInclusiveImplementation(StructureElements,InputElements,Radius,rResults);
      }
      
      void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }
      
      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusInclusiveImplementation(StructureNodes,InputNodes,Radius,rResults,rResultsDistance);
      }
      
      void SearchNodesInRadiusExclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(StructureNodes,InputNodes,Radius,rResults);
      }
      
      void SearchNodesInRadiusInclusive (
          NodesContainerType const& StructureNodes,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusInclusiveImplementation(StructureNodes,InputNodes,Radius,rResults);
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
          buffer << "DemSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DemSearch";}

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
      DEMSearch& operator=(DEMSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      DEMSearch(DEMSearch const& rOther)
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
        
 
//   /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
//                     DEMSearch& rThis){return rIStream;}
// 
//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
//                     const DEMSearch& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);
// 
//       return rOStream;
//     }
    
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DEM_SEARCH_H_INCLUDED  defined 


