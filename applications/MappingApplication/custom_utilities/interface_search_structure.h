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
#include "interface_object_manager.h"
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

      InterfaceSearchStructure(InterfaceObjectManager::Pointer i_interface_object_manager,
                               InterfaceObjectManager::Pointer i_interface_object_manager_bins,
                               ModelPart& i_model_part_bins) :
                               m_interface_object_manager(i_interface_object_manager),
                               m_interface_object_manager_bins(i_interface_object_manager_bins) {
          // Default Constructor, here the tolerances are computed in the search structure
          ComputeLocalSearchTolerance(i_model_part_bins);
          Initialize();
      }

      InterfaceSearchStructure(InterfaceObjectManager::Pointer i_interface_object_manager,
                               InterfaceObjectManager::Pointer i_interface_object_manager_bins,
                               double i_initial_search_radius, int i_max_search_iterations) :
                               m_interface_object_manager(i_interface_object_manager),
                               m_interface_object_manager_bins(i_interface_object_manager_bins),
                               m_initial_search_radius(i_initial_search_radius),
                               m_max_search_iterations(i_max_search_iterations) {
          // Custom Constructor, here the tolerances are passed from outside
          Initialize();
      }


      /// Destructor.
      virtual ~InterfaceSearchStructure(){}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      virtual void Search() {
          // TODO What happens in the case of remeshing? Reinitialize / redo the whole mapper?
          // TODO improve this fct towards a nicer structure / naming
          InterfaceObjectConfigure::ContainerType point_list = m_interface_object_manager->GetPointListSerialSearch();
          int num_points = point_list.size();

          std::vector<InterfaceObject::Pointer> points_results(num_points);
          std::vector<double> min_distances(num_points);

          FindLocalNeighbors(point_list, num_points, points_results, min_distances, m_comm_rank);

          m_interface_object_manager_bins->StoreTempSearchResults(points_results, m_comm_rank);
          m_interface_object_manager->PostProcessReceivedResults(min_distances, m_comm_rank);

          m_interface_object_manager->CheckResults();

          // PrintPairs();
      }

      virtual void Clear() {
          // TODO Implement


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

      InterfaceObjectManager::Pointer m_interface_object_manager;
      InterfaceObjectManager::Pointer m_interface_object_manager_bins;

      BinsObjectDynamic<InterfaceObjectConfigure>::Pointer m_local_bin_structure;
      int m_local_bin_structure_size;

      double m_initial_search_radius;
      int m_max_search_iterations;

    //   InterfaceObjectConfigure::ContainerType m_remote_p_point_list;
    //   std::vector<InterfaceObject::Pointer> m_candidate_receive_points;
    //   std::vector<double> m_min_distances;

      int m_comm_rank = 0; // default, for serial version
      int m_comm_size = 1; // default, for serial version

      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

      void FindLocalNeighbors(InterfaceObjectConfigure::ContainerType& point_list,
                              const int point_list_size, std::vector<InterfaceObject::Pointer>& points_results,
                              std::vector<double>& min_distances, const int comm_partner) {
          // This function finds neighbors of the points in point_list in bin_structure
          // It must be executable by serial and parallel version!
          // point_list_size must be passed bcs the point list might contain old entries (it has the max receive buffer size as size)!

          if (m_local_bin_structure == nullptr) { // this partition has no part of the sending interface, i.e. the origin of the mapped values
              for (int i = 0; i < point_list_size; ++i)
                  min_distances[i] = -1.0; // no results in this partition

          } else { // this partition has a bin structure

              InterfaceObjectConfigure::ResultContainerType neighbor_results(m_local_bin_structure_size);
              std::vector<double> neighbor_distances(m_local_bin_structure_size);

              InterfaceObjectConfigure::ResultIteratorType results_itr;
              std::vector<double>::iterator distance_itr;

              std::size_t number_of_results;

              InterfaceObjectConfigure::IteratorType point_itr;

              double search_radius;
              bool search_sucess;
            //   Searching the neighbors
            //   num_threads(variable)
            //   #pragma omp parallel for schedule(dynamic) if(point_list_size > m_omp_threshold_num_nodes) \
            //                            firstprivate(neighbor_results, neighbor_distances) \
            //                            private(point_itr, search_radius, search_sucess, \
            //                            results_itr, distance_itr, number_of_results)
              for (int i = 0; i < point_list_size; ++i){
                  point_itr = point_list.begin() + i;
                  search_radius = m_initial_search_radius; // reset search radius
                  search_sucess = false;

                  for (int j = 0; j < m_max_search_iterations; ++j) {
                      results_itr = neighbor_results.begin();
                      distance_itr = neighbor_distances.begin();

                      number_of_results = m_local_bin_structure->SearchObjectsInRadius(
                                  *point_itr, search_radius, results_itr,
                                  distance_itr, m_local_bin_structure_size);
                    //   std::cout << "number_of_results " << number_of_results << std::endl;
                      if (number_of_results > 0) { // neighbors were found
                          //std::cout << "Neighbor Search finished after " << j + 1 << " iterations" << std::endl;
                          search_sucess = true;
                          break;
                      }
                      search_radius *= 2; // double the search radius in case of an unsuccesful search
                  }

                  if (search_sucess) {
                      SelectBestResult(point_itr, neighbor_results, neighbor_distances, number_of_results,
                                             points_results[i], min_distances[i]);
                      // points_results[i].SelectBestResult(neighbor_results, neighbor_distances, number_of_results,
                      //                        points_results[i], min_distances[i]);
                  } else {
                      // TODO: Output some KRATOS Warning, but not here!
                      //std::cout << "WARNING: Node " << (*point_itr)->Id() << " did not find a neighbor!" << std::endl;
                      min_distances[i] = -1.0; // indicates that the search was not succesful
                      points_results[i] = nullptr; // TODO delete reference / unset pointer? "reset"
                  }
              }
          }
      }

      // TODO what abt elements? => how to destinguish nodes from elements here?
      // move this function to InterfaceNode and InterfaceElement?
      void SelectBestResult(InterfaceObjectConfigure::IteratorType point, const InterfaceObjectConfigure::ResultContainerType& result_list, const std::vector<double>& distances,
                                  const std::size_t& num_results, InterfaceObject::Pointer& vec_closest_results, double& closest_distance) {


          double min_distance = 1e15;
          int min_index = 0; // TODO this is sort of the safe option to avoid segfaults, discuss with Jordi

          for (int i = 0; i < static_cast<int>(num_results); ++i) { // find index of neighbor with smallest distance
              if (result_list[i]->CheckResult((*point)->Coordinates(),
                                              min_distance, distances[i])) {
                  min_index = i;
              }
          }
          vec_closest_results = result_list[min_index];
          closest_distance = min_distance;
      }

      void ComputeLocalSearchTolerance(ModelPart& model_part) {
          m_max_search_iterations = 10; // factor of 2^9 => 512 * m_initial_search_radius!
          double search_safety_factor = 1.2;
          double max_element_size = 0.0;

          // Loop through each edge of a geometrical entity ONCE
          for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
              for (std::size_t i = 0; i < (condition.GetGeometry().size() - 1); ++i) {
                  double node_1_x = condition.GetGeometry()[i].X();
                  double node_1_y = condition.GetGeometry()[i].Y();
                  double node_1_z = condition.GetGeometry()[i].Z();

                  for (std::size_t j = i + 1; j < condition.GetGeometry().size(); ++j) {
                      double node_2_x = condition.GetGeometry()[j].X();
                      double node_2_y = condition.GetGeometry()[j].Y();
                      double node_2_z = condition.GetGeometry()[j].Z();

                      double edge_length = sqrt(pow(node_1_x - node_2_x , 2) + pow(node_1_y - node_2_y , 2) + pow(node_1_z - node_2_z , 2));
                      max_element_size = std::max(max_element_size, edge_length);
                  }
              }
          }
          m_initial_search_radius = max_element_size * search_safety_factor;
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


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void Initialize() { // build the local bin-structure
          InterfaceObjectConfigure::ContainerType point_list_bins = m_interface_object_manager_bins->GetPointList();

          m_local_bin_structure_size = point_list_bins.size();

          if (m_local_bin_structure_size > 0) { // only construct the bins if the partition has a part of the interface
              m_local_bin_structure = BinsObjectDynamic<InterfaceObjectConfigure>::Pointer(
                  new BinsObjectDynamic<InterfaceObjectConfigure>(point_list_bins.begin(), point_list_bins.end()));
          }
      }

      void PrintPairs() {
          InterfaceObjectConfigure::ContainerType point_list = m_interface_object_manager->GetPointList();
          int num_points = point_list.size();

          InterfaceObjectConfigure::ContainerType point_list_bins = m_interface_object_manager_bins->GetPointListBins();
          // TODO warning size mismatch
          for (int i = 0; i < num_points; ++i) {
              std::cout << "Node1 " << point_list[i]->GetObjectId() <<
              " paired with Node2 " << point_list_bins[i]->GetObjectId() << std::endl;
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
