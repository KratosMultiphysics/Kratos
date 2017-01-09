//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher

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

  // TODO
  // typedef boost::shared_ptr<matrix<int> > GraphType; // GraphColoringProcess

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

      InterfaceSearchStructure(InterfaceObjectManagerBase::Pointer i_interface_object_manager,
                               InterfaceObjectManagerBase::Pointer i_interface_object_manager_bins,
                               int i_echo_level) :
                               m_interface_object_manager(i_interface_object_manager),
                               m_interface_object_manager_bins(i_interface_object_manager_bins) {
          m_echo_level = i_echo_level;
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

      void Search(const double i_search_radius, const int i_max_search_iterations) {
          m_search_radius = i_search_radius;
          m_max_search_iterations = i_max_search_iterations;
          // int increase_factor = 4;
          // int num_iteration = 1;

          // First Iteration is done outside the search loop bcs it has
          // to be done in any case
          // one search iteration should be enough in most cases (if the search
          // radius was either computed or specified properly)
          // only if some points dont have a valid projection, more search
          // iterations are necessary
          ConductSearchIteration();

          // while (num_iteration < m_max_search_iterations && !m_interface_object_manager->AllNeighborsFound()) {
          //     m_search_radius *= increase_factor;
          //     ++num_iteration;
          //
          //     if (m_comm_rank == 0) {
          //         std::cout << "MAPPER WARNING, search radius was increased, "
          //                   << "another search iteration is conducted, "
          //                   << "search iteration " << num_iteration
          //                   << ", search radius " << m_search_radius
          //                   << std::endl;
          //     }
          //
          //     // Clearing the managers is ugly, because all interface objects are
          //     // sent again, not just the ones that havent found a neighbor
          //     // m_interface_object_manager->Clear();
          //     // m_interface_object_manager_bins->Clear();
          //     ConductSearchIteration();
          // }

          m_interface_object_manager->CheckResults();

          // TODO
          // PrintPairs(); // makes only sense in serial, since there is no information locally existing abt the neighbor
          // Might be implemented at some point...
      }

      virtual void GetMPIData(int& max_send_buffer_size, int& max_receive_buffer_size,
                              GraphType& colored_graph, int& max_colors) {
          KRATOS_ERROR << "MappingApplication; InterfaceSearchStructure; "
                       << "\"GetMPIData\" of the base class (serial "
                       << "version) called!" << std::endl;
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

      InterfaceObjectManagerBase::Pointer m_interface_object_manager;
      InterfaceObjectManagerBase::Pointer m_interface_object_manager_bins;

      BinsObjectDynamic<InterfaceObjectConfigure>::Pointer m_local_bin_structure;
      int m_local_bin_structure_size;

      double m_search_radius;
      int m_max_search_iterations;

      int m_comm_rank = 0; // default, for serial version
      int m_comm_size = 1; // default, for serial version

      int m_echo_level = 0;

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

          if (m_local_bin_structure) { // this partition has a bin structure
              InterfaceObjectConfigure::ResultContainerType neighbor_results(m_local_bin_structure_size);
              std::vector<double> neighbor_distances(m_local_bin_structure_size);

              InterfaceObjectConfigure::IteratorType interface_object_itr;
              InterfaceObjectConfigure::ResultIteratorType results_itr;
              std::vector<double>::iterator distance_itr;

            //   Searching the neighbors
              for (int i = 0; i < interface_objects_size; ++i){
                  interface_object_itr = interface_objects.begin() + i;
                  double search_radius = m_search_radius; // reset search radius

                  results_itr = neighbor_results.begin();
                  distance_itr = neighbor_distances.begin();

                  std::size_t number_of_results = m_local_bin_structure->SearchObjectsInRadius(
                              *interface_object_itr, search_radius, results_itr,
                              distance_itr, m_local_bin_structure_size);

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
          InterfaceObjectConfigure::ContainerType interface_objects_bins = m_interface_object_manager_bins->GetInterfaceObjects();

          m_local_bin_structure_size = interface_objects_bins.size();

          if (m_local_bin_structure_size > 0) { // only construct the bins if the partition has a part of the interface
              m_local_bin_structure = BinsObjectDynamic<InterfaceObjectConfigure>::Pointer(
                  new BinsObjectDynamic<InterfaceObjectConfigure>(interface_objects_bins.begin(), interface_objects_bins.end()));
          }
      }

      virtual void ConductSearchIteration() {
          InterfaceObjectConfigure::ContainerType interface_objects;
          m_interface_object_manager->GetInterfaceObjectsSerialSearch(interface_objects);

          int num_objects = interface_objects.size();

          std::vector<InterfaceObject::Pointer> interface_object_results(num_objects);
          std::vector<double> min_distances(num_objects);
          std::vector<array_1d<double,2>> local_coordinates(num_objects);
          std::vector<std::vector<double>> shape_functions(num_objects);

          FindLocalNeighbors(interface_objects, num_objects, interface_object_results,
                             min_distances, local_coordinates, shape_functions, m_comm_rank);

          m_interface_object_manager_bins->StoreSearchResults(min_distances, interface_object_results, shape_functions, local_coordinates);
          m_interface_object_manager->PostProcessReceivedResults(min_distances, interface_objects);
      }

      virtual void PrintPairs() {
          std::vector<InterfaceObject::Pointer> interface_objects = m_interface_object_manager->GetDestinationInterfaceObjects();
          int num_objects = interface_objects.size();

          std::vector<InterfaceObject::Pointer> interface_objects_bins = m_interface_object_manager_bins->GetOriginInterfaceObjects();
          if (num_objects != static_cast<int>(interface_objects_bins.size()))
              KRATOS_ERROR << "MappingApplication; InterfaceSearchStructure; \"PrintPairs\" size mismatch!" << std::endl;
          for (int i = 0; i < num_objects; ++i) {
              if (interface_objects_bins[i] == nullptr) {
                  std::cout << "InterfaceObject 1 ";
                  interface_objects[i]->PrintMatchInfo();
                  std::cout << " has not found a match! " << std::endl;
              } else {
                  interface_objects[i]->PrintMatchInfo();
                  std::cout << " paired with ";
                  interface_objects_bins[i]->PrintMatchInfo();
                  std::cout << std::endl;
              }
          }
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
