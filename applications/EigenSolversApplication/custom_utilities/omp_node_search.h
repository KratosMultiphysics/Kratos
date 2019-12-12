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
#include "utilities/openmp_utils.h"
#include "spatial_containers/bins_dynamic_objects.h"

// Configures
#include "node_configure.h"

// Search
#include "spatial_containers/point_search.h"
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

class OMP_NodeSearch 
{
    public:
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of OMP_NodeSearch
      KRATOS_CLASS_POINTER_DEFINITION(OMP_NodeSearch);
      
      //Configure Types
      typedef NodeConfigure<3>                              NodeConfigureType;          //Node
      typedef ModelPart::NodesContainerType                 NodesContainerType;
      //Bin Types
      typedef BinsObjectDynamic<NodeConfigureType>      NodeBinsType;
     
      typedef NodesContainerType::ContainerType             ResultNodesContainerType;
      typedef std::vector<ResultNodesContainerType>         VectorResultNodesContainerType;

      typedef std::vector<double>                           RadiusArrayType;

      typedef std::vector<double>                           DistanceType;
      typedef std::vector<DistanceType>                     VectorDistanceType;
      
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.

      OMP_NodeSearch(const double domain_min_x = 0.0, const double domain_min_y = 0.0, const double domain_min_z = 0.0,
                    const double domain_max_x = -1.0, const double domain_max_y = -1.0, const double domain_max_z = -1.0)
      {
          mIsInitialized = false;
      }

      /// Destructor.
      ~OMP_NodeSearch(){
      }
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

      void InitializeSearch(NodesContainerType const& rStructureNodes)
      {
        KRATOS_TRY;
        if(!mIsInitialized) {
            NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
            mBins = new NodeBinsType(nodes_ModelPart.begin(), nodes_ModelPart.end());

            mIsInitialized = true;
        }
        KRATOS_CATCH("");
      }
      /* 
      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          NodesContainerType const& rNodes,
          const RadiusArrayType & Radius, 
          VectorResultNodesContainerType& rResults )
      {     
          KRATOS_TRY
          int MaxNumberOfNodes = rStructureNodes.size();
          
          NodesContainerType::ContainerType& nodes_array = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());

          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(MaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;
       
              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
                
                  NumberOfResults = mBins->SearchObjectsInRadiusExclusive(nodes_array[i],Radius[i],ResultsPointer,MaxNumberOfNodes);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);     
              }
          }
          
          KRATOS_CATCH("")
      }*/

      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          int const Id,
          double const Radius, 
          ResultNodesContainerType& rResults )
      {     
            KRATOS_TRY
            int MaxNumberOfNodes = rStructureNodes.size();
            
            NodesContainerType::ContainerType& nodes_array = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());

            ResultNodesContainerType  localResults(MaxNumberOfNodes);
            std::size_t               NumberOfResults = 0;
            
            ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();
        
            NumberOfResults = mBins->SearchObjectsInRadiusExclusive( nodes_array[Id],Radius,ResultsPointer,MaxNumberOfNodes);
            
            rResults.insert(rResults.begin(),localResults.begin(),localResults.begin()+NumberOfResults);     
                    
          
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
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "OpenMPDemSearch";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}
      
            
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
      NodeBinsType* mBins;

      bool mIsInitialized;
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
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
      OMP_NodeSearch& operator=(OMP_NodeSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      OMP_NodeSearch(OMP_NodeSearch const& rOther)
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


