//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_object.h"
#include "interface_object_manager_serial.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "custom_configures/interface_object_configure.h"

namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  typedef matrix<int> GraphType; // GraphColoringProcess

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
  class InterfaceSearchStructure
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceSearchStructure
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceSearchStructure);

      ///@}
      ///@name Life Cycle
      ///@{

      InterfaceSearchStructure(InterfaceObjectManagerBase::Pointer pInterfaceObjectManager,
                               InterfaceObjectManagerBase::Pointer pInterfaceObjectManagerBins,
                               int EchoLevel) :
                               mpInterfaceObjectManager(pInterfaceObjectManager),
                               mpInterfaceObjectManagerBins(pInterfaceObjectManagerBins) {
          mEchoLevel = EchoLevel;
          Initialize();
      }


      /// Destructor.
      virtual ~InterfaceSearchStructure(){ }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Search(const double SearchRadius, const int MaxSearchIterations) {
          mSearchRadius = SearchRadius;
          mMaxSearchIterations = MaxSearchIterations;
          int increase_factor = 4;
          int num_iteration = 1;
          bool last_iteration = false;

          if (mMaxSearchIterations == 1) { // in case only one search iteration is conducted
              last_iteration = true;
          }

          // First Iteration is done outside the search loop bcs it has
          // to be done in any case
          // one search iteration should be enough in most cases (if the search
          // radius was either computed or specified properly)
          // only if some points did not find a neighbor or dont have a valid
          // projection, more search iterations are necessary
          ConductSearchIteration(last_iteration);

          while (num_iteration < mMaxSearchIterations && !mpInterfaceObjectManager->AllNeighborsFound()) {
              mSearchRadius *= increase_factor;
              ++num_iteration;

              if (num_iteration == mMaxSearchIterations) {
                  last_iteration = true;
              }

              if (mEchoLevel > 1 && mCommRank == 0) {
                  std::cout << "MAPPER WARNING, search radius was increased, "
                            << "another search iteration is conducted, "
                            << "search iteration " << num_iteration << " / "
                            << mMaxSearchIterations << ", search radius "
                            << mSearchRadius << std::endl;
              }

              ConductSearchIteration(last_iteration);
          }
          if (mEchoLevel > 1) {
              mpInterfaceObjectManager->CheckResults();
          }
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
        buffer << "InterfaceSearchStructure" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceSearchStructure";}

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

      InterfaceObjectManagerBase::Pointer mpInterfaceObjectManager;
      InterfaceObjectManagerBase::Pointer mpInterfaceObjectManagerBins;

      BinsObjectDynamic<InterfaceObjectConfigure>::Pointer mpLocalBinStructure;
      int mLocalBinStructureSize;

      double mSearchRadius;
      int mMaxSearchIterations;

      int mCommRank = 0; // default, for serial version
      int mCommSize = 1; // default, for serial version

      int mEchoLevel = 0;

      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

      void FindLocalNeighbors(InterfaceObjectConfigure::ContainerType& interface_objects,
                              const int interface_objects_size, std::vector<InterfaceObject::Pointer>& interface_object_results,
                              std::vector<double>& min_distances, std::vector<array_1d<double,2>>& local_coordinates,
                              std::vector<std::vector<double>>& shape_functions, const int comm_partner) {
          // This function finds neighbors of the InterfaceObjects in interface_objects in bin_structure
          // It must be executable by serial and parallel version!
          // interface_objects_size must be passed bcs interface_objects might contain old entries (it has the max receive buffer size as size)!

          if (mpLocalBinStructure) { // this partition has a bin structure
              InterfaceObjectConfigure::ResultContainerType neighbor_results(mLocalBinStructureSize);
              std::vector<double> neighbor_distances(mLocalBinStructureSize);

              InterfaceObjectConfigure::IteratorType interface_object_itr;
              InterfaceObjectConfigure::ResultIteratorType results_itr;
              std::vector<double>::iterator distance_itr;

            //   Searching the neighbors
              for (int i = 0; i < interface_objects_size; ++i){
                  interface_object_itr = interface_objects.begin() + i;
                  double search_radius = mSearchRadius; // reset search radius

                  results_itr = neighbor_results.begin();
                  distance_itr = neighbor_distances.begin();

                  std::size_t number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
                              *interface_object_itr, search_radius, results_itr,
                              distance_itr, mLocalBinStructureSize);

                  if (number_of_results > 0) { // neighbors were found
                      SelectBestResult(interface_object_itr, neighbor_results,
                                       neighbor_distances, number_of_results,
                                       interface_object_results[i], min_distances[i],
                                       local_coordinates[i], shape_functions[i]);
                  } else {
                      min_distances[i] = -1.0f; // indicates that the search was not succesful
                      interface_object_results[i].reset(); // Release an old pointer, that is probably existing from a previous search
                  }
              }
          } else { // this partition has no part of the point receiving interface, i.e. the origin of the mapped values
              for (int i = 0; i < interface_objects_size; ++i) { // no results in this partition
                  min_distances[i] = -1.0f; // indicates that the search was not succesful
                  interface_object_results[i].reset(); // Release an old pointer, that is probably existing from a previous search
              }
          }
      }

      void SelectBestResult(const InterfaceObjectConfigure::IteratorType& point,
                            const InterfaceObjectConfigure::ResultContainerType& result_list,
                            const std::vector<double>& distances, const std::size_t num_results,
                            InterfaceObject::Pointer& vec_closest_results, double& closest_distance,
                            array_1d<double,2>& local_coords, std::vector<double>& shape_functions_values) {

          double min_distance = std::numeric_limits<double>::max();
          closest_distance = -1.0f; // indicate a failed search in case no result is good
          array_1d<double,2> local_coords_temp;

          for (int i = 0; i < static_cast<int>(num_results); ++i) { // find index of best result
              if (result_list[i]->EvaluateResult((*point)->Coordinates(), min_distance, distances[i],
                                                 local_coords_temp, shape_functions_values)) {
                  closest_distance = min_distance;
                  local_coords = local_coords_temp;
                  vec_closest_results = result_list[i];
              }
          }
      }


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

      // int m_omp_threshold_num_nodes = 1000; // TODO constexpr???

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void Initialize() { // build the local bin-structure
          InterfaceObjectConfigure::ContainerType interface_objects_bins = mpInterfaceObjectManagerBins->GetInterfaceObjects();

          mLocalBinStructureSize = interface_objects_bins.size();

          if (mLocalBinStructureSize > 0) { // only construct the bins if the partition has a part of the interface
              mpLocalBinStructure = BinsObjectDynamic<InterfaceObjectConfigure>::Pointer(
                  new BinsObjectDynamic<InterfaceObjectConfigure>(interface_objects_bins.begin(), interface_objects_bins.end()));
          }
      }

      virtual void ConductSearchIteration(const bool last_iteration) {
          InterfaceObjectConfigure::ContainerType interface_objects;
          mpInterfaceObjectManager->GetInterfaceObjectsSerialSearch(interface_objects);

          int num_objects = interface_objects.size();

          std::vector<InterfaceObject::Pointer> interface_object_results(num_objects);
          std::vector<double> min_distances(num_objects);
          std::vector<array_1d<double,2>> local_coordinates(num_objects);
          std::vector<std::vector<double>> shape_functions(num_objects);

          FindLocalNeighbors(interface_objects, num_objects, interface_object_results,
                             min_distances, local_coordinates, shape_functions, mCommRank);

          mpInterfaceObjectManagerBins->StoreSearchResults(min_distances, interface_object_results, shape_functions, local_coordinates);
          mpInterfaceObjectManager->PostProcessReceivedResults(min_distances, interface_objects);
      }


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
      InterfaceSearchStructure& operator=(InterfaceSearchStructure const& rOther);

    //   /// Copy constructor.
    //   InterfaceSearchStructure(InterfaceSearchStructure const& rOther){}


      ///@}

    }; // Class InterfaceSearchStructure

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceSearchStructure& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceSearchStructure& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_STRUCTURE_H_INCLUDED  defined
