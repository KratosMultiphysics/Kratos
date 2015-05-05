//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_POINT_SEARCH_H_INCLUDED )
#define  KRATOS_POINT_SEARCH_H_INCLUDED

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
class PointSearch : public SpatialSearch
{
    public:
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(PointSearch);
    
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      PointSearch(){}

      /// Destructor.
      virtual ~PointSearch(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
        
      void SearchNodesInRadiusExclusive (
          ModelPart& rModelPart,
          NodesContainerType const& InputNodes,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
          static_cast<TDerived*>(this)->SearchNodesInRadiusExclusiveImplementation(rModelPart,InputNodes,Radius,rResults,rResultsDistance);
      }
      
      void SearchElementsInRadiusInclusive (
          ModelPart& rModelPart,
          ElementsContainerType const& InputElements,
          const RadiusArrayType & Radius,
          VectorResultNodesContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_THROW_ERROR(std::runtime_error,"Elemental search is not defined for point based structures","");
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
          buffer << "PointSearch" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "PointSearch";}

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
      PointSearch& operator=(PointSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      PointSearch(PointSearch const& rOther)
      {
          *this = rOther;
      }

        
      ///@}    
        
    }; // Class PointSearch

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
//                     PointSearch& rThis){return rIStream;}
// 
//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
//                     const PointSearch& rThis)
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

#endif // KRATOS_POINT_SEARCH_H_INCLUDED  defined 


