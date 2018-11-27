//
//   Project Name:        Kratos
//   Last Modified by:    $Author: msantasusana, croig $
//   Date:                $Date: 2015-10-26 19:37:47 $
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
#include "rigid_face_geometrical_object_configure.h"
// Search
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"

// External includes
//#define CUSTOMTIMER

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

class KRATOS_API(DEM_APPLICATION) DEM_FEM_Search : public SpatialSearch
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

//       //Configure Types
//       typedef RigidFaceConfigure<3>                        ElementConfigureType;   //Element
//       //Bin Types
//       typedef BinsObjectDynamic<ElementConfigureType>   BinsType;
//
    //Configure Types
    typedef RigidFaceGeometricalObjectConfigure<3>        RigidFaceGeometricalConfigureType;
      //Bin Types
    typedef BinsObjectDynamic<RigidFaceGeometricalConfigureType>   GeometricalBinsType;
    typedef PointerVectorSet<GeometricalObject, IndexedObject>     GeometricalObjectType;
    //typedef PointerVector<GeometricalObject>     GeometricalObjectType;


      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      DEM_FEM_Search(){

        mBins = NULL;

      }

      /// Destructor.
      ~DEM_FEM_Search(){
      }

  void SearchRigidFaceForDEMInRadiusExclusiveImplementation (
      ElementsContainerType   const& rElements,
      ConditionsContainerType const& rConditions,
      VectorResultConditionsContainerType& rResults,
      VectorDistanceType& rResultsDistance)
    {
      KRATOS_TRY

      /*
        STEPS:
        ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
        1. INITIALIZE
        2. CALCULATE THE DEM BBX
        3. GATHER DEM BOUNDING BOX
        4. FIND THE FE INSIDE DEM_BB TO BUILD THE BINS AND CONSTRUCT THE GLOBAL BBX
        5. GATHER FEM ELEMENTS AND THE GLOBAL BBX
        6. BUILD THE BINS
        7. PERFORM THE SEARCH FOR THE SPHERES INSIDE THE BBX (amplified by the radius of each particle)
      */

      //1. INITIALIZE

      int MaxNumberOfElements = rConditions.size();

      ElementsContainerType::ContainerType& elements_sear   = const_cast<ElementsContainerType::ContainerType&>  (rElements.GetContainer());
      ConditionsContainerType::ContainerType& conditions_bins = const_cast<ConditionsContainerType::ContainerType&>(rConditions.GetContainer());

      //GeometricalObjectType::ContainerType SearElementPointerToGeometricalObjecPointerTemporalVector;
      GeometricalObjectType::ContainerType BinsConditionPointerToGeometricalObjecPointerTemporalVector;
      RadiusArrayType Radius_out;

      int num_of_threads = OpenMPUtils::GetNumThreads();
      DenseVector<unsigned int> total_dem_partition_index; vector<unsigned int> total_fem_partition_index;

      OpenMPUtils::CreatePartition(num_of_threads, elements_sear.size(), total_dem_partition_index);
      OpenMPUtils::CreatePartition(num_of_threads, conditions_bins.size(), total_fem_partition_index);

      //std::vector<GeometricalObjectType::ContainerType> Vector_SearElementPointerToGeometricalObjecPointerTemporalVector(num_of_threads);
      std::vector<GeometricalObjectType::ContainerType> Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector(num_of_threads);

      std::vector<array_1d<double, 3> > Vector_DEM_BB_LowPoint(num_of_threads); std::vector <array_1d<double, 3 > > Vector_DEM_BB_HighPoint(num_of_threads);
      std::vector<array_1d<double, 3> > Vector_GLOBAL_BB_LowPoint(num_of_threads); std::vector <array_1d<double, 3 > > Vector_GLOBAL_BB_HighPoint(num_of_threads);

      std::vector<double> Vector_Ref_Radius(num_of_threads);
      std::vector<RadiusArrayType> Vector_Radius_out(num_of_threads);

      double Global_Ref_Radius = 0.0;
      double inf = std::numeric_limits<double>::infinity();

      for (std::size_t i = 0; i < 3; i++) {
        DEM_BB_LowPoint[i]      = inf;
        DEM_BB_HighPoint[i]     = -inf;

        mGlobal_BB_LowPoint[i]  = inf;
        mGlobal_BB_HighPoint[i] = -inf;
      }

      typedef ElementsContainerType::ContainerType::iterator   Elem_iter;
      typedef ConditionsContainerType::ContainerType::iterator Cond_iter;

      //2. CALCULATE THE DEM BBX

      #pragma omp parallel
      {
        double radius = 0.0;
        int k = OpenMPUtils::ThisThread();

        for(std::size_t i = 0; i < 3; i++) {
            Vector_DEM_BB_LowPoint[k][i]  = inf;
            Vector_DEM_BB_HighPoint[k][i] = -inf;
        }

        #pragma omp for
        for (int p = 0; p <(int) elements_sear.size(); p++) {

            Elem_iter it = elements_sear.begin() + p;
            GeometryType& pGeometry = (*it)->GetGeometry();
            const array_1d<double, 3 >& aux_coor = pGeometry[0].Coordinates();

            SphericParticle* p_particle = dynamic_cast<SphericParticle*>((*it).get());
            radius = p_particle->GetSearchRadius();

            Vector_Ref_Radius[k]    = (Vector_Ref_Radius[k]  < radius) ? radius : Vector_Ref_Radius[k] ;

            for(std::size_t i = 0; i < 3; i++) {
              Vector_DEM_BB_LowPoint[k][i]   = (Vector_DEM_BB_LowPoint[k][i] > aux_coor[i]) ? aux_coor[i] : Vector_DEM_BB_LowPoint[k][i];
              Vector_DEM_BB_HighPoint[k][i]  = (Vector_DEM_BB_HighPoint[k][i] < aux_coor[i]) ? aux_coor[i] : Vector_DEM_BB_HighPoint[k][i];
            }
          } //pragma
      }//pragma parallel


      //3. GATHER DEM BOUNDING BOX
      for(int k = 0; k < num_of_threads; k++) {
        for(std::size_t i = 0; i < 3; i++) {
          DEM_BB_LowPoint[i]  = (DEM_BB_LowPoint[i] > Vector_DEM_BB_LowPoint[k][i]) ? Vector_DEM_BB_LowPoint[k][i] : DEM_BB_LowPoint[i];
          DEM_BB_HighPoint[i] = (DEM_BB_HighPoint[i] < Vector_DEM_BB_HighPoint[k][i]) ? Vector_DEM_BB_HighPoint[k][i] : DEM_BB_HighPoint[i];
        }

        Global_Ref_Radius = (Global_Ref_Radius < Vector_Ref_Radius[k]) ? Vector_Ref_Radius[k] : Global_Ref_Radius;
      }

      for(std::size_t i = 0; i < 3; i++) {
        DEM_BB_LowPoint[i]  -= 1.00f * Global_Ref_Radius;
        DEM_BB_HighPoint[i] += 1.00f * Global_Ref_Radius;
      }


       //4. FIND THE FE INSIDE DEM_BB TO BUILD THE BINS AND CONSTRUCT THE GLOBAL BBX

      #pragma omp parallel
      {
        int k = OpenMPUtils::ThisThread();
        Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector[k].reserve(total_fem_partition_index[k+1]);

        for(std::size_t i = 0; i < 3; i++) {
          Vector_GLOBAL_BB_LowPoint[k][i]  = inf;
          Vector_GLOBAL_BB_HighPoint[k][i] = -inf;
        }

        array_1d<double, 3> rHighPoint;
        array_1d<double, 3> rLowPoint;

        #pragma omp for private(rHighPoint,rLowPoint)
        for (int c = 0; c < (int)conditions_bins.size(); c++) {

          Cond_iter it = conditions_bins.begin() + c;

          const GeometryType& pGeometry = (*it)->GetGeometry();
          noalias(rLowPoint)  = pGeometry[0];
          noalias(rHighPoint) = pGeometry[0];

          for(unsigned int point = 1; point < pGeometry.size(); point++ ) {
            for(unsigned int i = 0; i < 3; i++ ) {
              rHighPoint[i] = ( rHighPoint[i] < pGeometry[point][i] ) ? pGeometry[point][i] : rHighPoint[i];
              rLowPoint[i]  = ( rLowPoint[i]  > pGeometry[point][i] ) ? pGeometry[point][i] : rLowPoint[i];
            }
          }

          bool add = true;

          for(unsigned int i = 0; i < 3; i++) {
            if(( rHighPoint[i]  < DEM_BB_LowPoint[i] ) || ( rLowPoint[i]  > DEM_BB_HighPoint[i] )) {
              add = false;
              break;
            }
          }

          if(add) {
            for(unsigned int i = 0; i < 3; i++ ) {
              Vector_GLOBAL_BB_LowPoint[k][i]   = (Vector_GLOBAL_BB_LowPoint[k][i] > rLowPoint[i]) ? rLowPoint[i] : Vector_GLOBAL_BB_LowPoint[k][i];
              Vector_GLOBAL_BB_HighPoint[k][i]  = (Vector_GLOBAL_BB_HighPoint[k][i] < rHighPoint[i]) ? rHighPoint[i] : Vector_GLOBAL_BB_HighPoint[k][i];
            }
            Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector[k].push_back(*it);
          }
        }//parallel for

      }//parallel omp

      //5. GATHER FEM ELEMENTS AND THE GLOBAL BBX
      int fem_total_size = 0;
      for(int k = 0; k < num_of_threads; k++) {
        fem_total_size += Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector[k].size();
      }

      BinsConditionPointerToGeometricalObjecPointerTemporalVector.reserve(fem_total_size);

      for(int k = 0; k < num_of_threads; k++) {
        BinsConditionPointerToGeometricalObjecPointerTemporalVector.insert(
          BinsConditionPointerToGeometricalObjecPointerTemporalVector.end(),
          Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector[k].begin(),
          Vector_BinsConditionPointerToGeometricalObjecPointerTemporalVector[k].end()
        );

        for(std::size_t i = 0; i < 3; i++) {
          mGlobal_BB_LowPoint[i]  = (mGlobal_BB_LowPoint[i] > Vector_GLOBAL_BB_LowPoint[k][i]) ? Vector_GLOBAL_BB_LowPoint[k][i] : mGlobal_BB_LowPoint[i];
          mGlobal_BB_HighPoint[i] = (mGlobal_BB_HighPoint[i] < Vector_GLOBAL_BB_HighPoint[k][i]) ? Vector_GLOBAL_BB_HighPoint[k][i] : mGlobal_BB_HighPoint[i];
        }
      }

      if(BinsConditionPointerToGeometricalObjecPointerTemporalVector.size() >0 ) {

      //6. CREATE THE BINS

      //if (mBins) free(mBins);
      delete mBins;
      mBins = new GeometricalBinsType(BinsConditionPointerToGeometricalObjecPointerTemporalVector.begin(), BinsConditionPointerToGeometricalObjecPointerTemporalVector.end());

      //7. PERFORM THE SEARCH ON THE SPHERES
      #pragma omp parallel
      {
        GeometricalObjectType::ContainerType  localResults(MaxNumberOfElements);
        DistanceType                          localResultsDistances(MaxNumberOfElements);
        std::size_t                           NumberOfResults = 0;

        #pragma omp for schedule(dynamic, 100)
        for (int p = 0; p < (int)elements_sear.size(); p++) {

          Elem_iter it = elements_sear.begin() + p;

            GeometricalObject::Pointer go_it(*it);
            bool search_particle = true;

            array_1d<double, 3 > & aux_coor = go_it->GetGeometry()[0].Coordinates();

            SphericParticle* p_particle = dynamic_cast<SphericParticle*>((*it).get());
            double Rad = p_particle->GetSearchRadius();

            for(unsigned int i = 0; i < 3; i++ ) {
              search_particle &= !(aux_coor[i]  < (mGlobal_BB_LowPoint[i] - Rad) ) || (aux_coor[i]  > (mGlobal_BB_HighPoint[i] + Rad) ); //amplify the BBX with the radius for every particle
            }

            if(search_particle) {

              GeometricalObjectType::ContainerType::iterator   ResultsPointer          = localResults.begin();
              DistanceType::iterator                           ResultsDistancesPointer = localResultsDistances.begin();

              NumberOfResults = (*mBins).SearchObjectsInRadiusExclusive(go_it,Rad,ResultsPointer,ResultsDistancesPointer,MaxNumberOfElements);

              rResults[p].reserve(NumberOfResults);

              for(GeometricalObjectType::ContainerType::iterator c_it = localResults.begin(); c_it != localResults.begin() + NumberOfResults; c_it++) {
                  Condition::Pointer condition = Kratos::dynamic_pointer_cast<Condition>(*c_it);
                  rResults[p].push_back(condition);
              }

              rResultsDistance[p].insert(rResultsDistance[p].begin(),localResultsDistances.begin(),localResultsDistances.begin()+NumberOfResults);

            }

          } //loop elements
        } //parallel

    }//if Bins is not empty--> search

    KRATOS_CATCH("")
  }


  array_1d<double, 3 > GetBBHighPoint() {
    return (mGlobal_BB_HighPoint);
  }
  array_1d<double, 3 > GetBBLowPoint() {
    return (mGlobal_BB_LowPoint);
  }

  /// Turn back information as a string.
  virtual std::string Info() const override
  {
      std::stringstream buffer;
      buffer << "DEM_FEM_Search" ;

      return buffer.str();
  }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "DEM_FEM_Search";}

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
      array_1d<double, 3 > mGlobal_BB_HighPoint;
      array_1d<double, 3 > mGlobal_BB_LowPoint;

      array_1d<double, 3 > DEM_BB_HighPoint;
      array_1d<double, 3 > DEM_BB_LowPoint;
      GeometricalBinsType* mBins;




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
