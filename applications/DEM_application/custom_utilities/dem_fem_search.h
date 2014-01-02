//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_DEM_FEM_SEARCH_H_INCLUDED )
#define  KRATOS_DEM_FEM_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "utilities/openmp_utils.h"

// Configures
#include "dem_fem_configure.h"
#include "rigid_face_geometrical_object_configure.h"
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

class DEM_FEM_Search : public SpatialSearch
{
    public:
      ///@name Type Definitions
      ///@{
    
      /// Pointer definition of OMP_DEMSearch
      KRATOS_CLASS_POINTER_DEFINITION(DEM_FEM_Search);
      
      typedef PointType*                                PtrPointType;
      typedef std::vector<PtrPointType>*                PointVector;
      typedef std::vector<PtrPointType>::iterator       PointIterator;
      
      typedef double*                                   DistanceVector;
      typedef double*                                   DistanceIterator;
      
      //Configure Types
      typedef RigidFaceConfigure<3>                        ElementConfigureType;   //Element
      //Bin Types
      typedef BinsObjectDynamic<ElementConfigureType>   BinsType;
	  
	  //Configure Types
	  typedef RigidFaceGeometricalObjectConfigure<3>        RigidFaceGeometricalConfigureType;     	
      //Bin Types
	  typedef BinsObjectDynamic<RigidFaceGeometricalConfigureType>   GeometricalBinsType;
      typedef PointerVectorSet<GeometricalObject, IndexedObject>     GeometricalObjectType;
	  
	  
     
      
      ///@}
      ///@name Life Cycle 
      ///@{
      
      /// Default constructor.
      DEM_FEM_Search(){
      }

      /// Destructor.
      ~DEM_FEM_Search(){
      }
	  
	  
	   void SearchRigidFaceForDEMInRadiusExclusiveImplementation (
          ElementsContainerType   const& rStructureElements,
          ConditionsContainerType const& rElements,
          const RadiusArrayType & Radius, 
          VectorResultConditionsContainerType& rResults, 
          VectorDistanceType& rResultsDistance )
      {     
          KRATOS_TRY
          
          int MaxNumberOfElements = rStructureElements.size();

          ElementsContainerType::ContainerType& elements_sear   = const_cast<ElementsContainerType::ContainerType&>  (rStructureElements.GetContainer());
          ConditionsContainerType::ContainerType& elements_bins = const_cast<ConditionsContainerType::ContainerType&>(rElements.GetContainer());

          GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
          GeometricalObjectType::ContainerType BinsElementPointerToGeometricalObjecPointerTemporalVector;

          SearElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_sear.size());
          BinsElementPointerToGeometricalObjecPointerTemporalVector.reserve(elements_bins.size());
              
          for(ConditionsContainerType::ContainerType::iterator it = elements_bins.begin(); it != elements_bins.end(); it++)
              BinsElementPointerToGeometricalObjecPointerTemporalVector.push_back(*it);
          
          for(ElementsContainerType::ContainerType::iterator it = elements_sear.begin(); it != elements_sear.end(); it++)
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
                      Condition::Pointer elem = boost::dynamic_pointer_cast<Condition>(*it);
                      rResults[i].push_back(elem);
                      rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);
                  }
              }
          }
          
          KRATOS_CATCH("")
      }
	  
	  
	  
	  
	 void SearchConditionsInRadiusExclusive (
          ModelPart& rModelPart,
          ConditionsContainerType const& TConditions,
          const RadiusArrayType & Radius,
          VectorResultConditionsContainerType& rResults,
          VectorDistanceType& rResultsDistance )
      {
      
	      KRATOS_TRY
		  
		  ConditionsContainerType const& MConditions = rModelPart.GetCommunicator().LocalMesh().Conditions();
		
          
          int MaxNumberOfElements = TConditions.size();
          
          ConditionsContainerType::ContainerType& elements_array     = const_cast<ConditionsContainerType::ContainerType&>(MConditions.GetContainer());
          ConditionsContainerType::ContainerType& conditions_array   = const_cast<ConditionsContainerType::ContainerType&>(TConditions.GetContainer());
        
          BinsType bins(conditions_array.begin(), conditions_array.end());
		  
		  typedef ConditionsContainerType::ContainerType::value_type ConditionsPointerType;
		  
		  //KRATOS_WATCH(elements_array.size());
		  //KRATOS_WATCH(conditions_array.size());
                    
          #pragma omp parallel
          {
              ResultConditionsContainerType   localResults(MaxNumberOfElements);
              DistanceType                  localResultsDistances(MaxNumberOfElements);
              std::size_t                   NumberOfResults = 0;
              
              #pragma omp for
              for(int i = 0; i < static_cast<int>(elements_array.size()); i++)
              {
                  ResultConditionsContainerType::iterator ResultsPointer          = localResults.begin();
                  DistanceType::iterator                ResultsDistancesPointer = localResultsDistances.begin();
                
                  NumberOfResults = bins.SearchObjectsInRadiusExclusive(elements_array[i],Radius[i],ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);
                  
                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
                  rResultsDistance[i].insert(rResultsDistance[i].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);   

              }
          }
          
          KRATOS_CATCH("")	  
      }
      
	

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "DEM_FEM_Search" ;
          
          return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "DEM_FEM_Search";}

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
     DEM_FEM_Search& operator=(DEM_FEM_Search const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      DEM_FEM_Search(DEM_FEM_Search const& rOther)
      {
          *this = rOther;
      }

       
        
        
    }; // Class DEM_FEMSearch

  
}  // namespace Kratos.

#endif // KRATOS_DEM_FEM_SEARCH_H_INCLUDED  defined 


