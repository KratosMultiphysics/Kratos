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

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_object_manager_base.h"
#include "mapper_utilities_mpi.h"


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
  class InterfaceObjectManagerParallel : public InterfaceObjectManagerBase
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceObjectManagerParallel
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerParallel);

      ///@}
      ///@name Life Cycle
      ///@{

      InterfaceObjectManagerParallel(ModelPart& i_model_part, int i_comm_rank, int i_comm_size,
                                     MapperUtilities::InterfaceObjectConstructionType i_interface_object_type,
                                     GeometryData::IntegrationMethod i_integration_method, int i_echo_level) :
                                     InterfaceObjectManagerBase(i_model_part, i_comm_rank, i_comm_size,
                                     i_interface_object_type, i_integration_method, i_echo_level) {}

      /// Destructor.
      virtual ~InterfaceObjectManagerParallel() {
          delete [] m_bounding_box;
      }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      // **********************************************************************
      // Side we want to find neighbors for aka destination *******************
      // **********************************************************************
      void CheckResults() override {
          for (auto& interface_obj : m_interface_objects) {
              if (!interface_obj->NeighborFound()) {
                  std::cout << "MAPPER WARNING, Rank " << m_comm_rank
                            << "\tPoint has not found a neighbor, [ "
                            << interface_obj->X() << " | "
                            << interface_obj->Y() << " | "
                            << interface_obj->Z() << " ]" << std::endl;
              }
          }
      }

      void ComputeCandidatePartitions(CandidateManager& candidate_manager, int* local_comm_list,
                                      int* local_memory_size_array,
                                      double* global_bounding_boxes,
                                      const bool last_iteration) override {
          double* bounding_box[6];
          std::vector<int> partition_list; // For debugging
          for (auto& interface_obj : m_interface_objects) {
              for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) { // loop over partitions
                  for (int j = 0; j < 6; ++j) { // retrieve bounding box of partition
                      bounding_box[j] = &global_bounding_boxes[(partition_index * 6) + j];
                  }
                  if (!interface_obj->NeighborFound()) { // check if the interface object already found a neighbor
                      if (interface_obj->IsInBoundingBox(bounding_box)) {
                          candidate_manager.candidate_send_objects[partition_index].push_back(interface_obj);
                          local_comm_list[partition_index] = 1;
                          ++local_memory_size_array[partition_index];
                          interface_obj->SetIsBeingSent();
                          partition_list.push_back(partition_index); // For debugging
                      }
                  }
              }

              if (m_echo_level > 2) {
                  PrintCandidatePartitions(interface_obj, partition_list); // For debugging
              }

              if (last_iteration) {
                  // Robustness check, if interface_obj is sent to at least one partition
                  if (!interface_obj->GetIsBeingSent()) {
                      // Send interface_obj to all Partitions
                      std::cout << "MAPPER WARNING, Rank " << m_comm_rank
                                << ", interface_obj [ " << interface_obj->X()
                                << " " << interface_obj->Y() << " "
                                << interface_obj->Z() << " ] was not in a "
                                << "bounding box and is sent to all partitions!"
                                << std::endl;

                      for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) {
                          candidate_manager.candidate_send_objects[partition_index].push_back(interface_obj);
                          local_comm_list[partition_index] = 1;
                          ++local_memory_size_array[partition_index];
                      }
                  }
              }
          }
      }

      // Function for debugging
      void PrintCandidatePartitions(const InterfaceObject::Pointer interface_obj,
                                    std::vector<int>& partition_list) override {
          std::cout << "Rank " << m_comm_rank << ", interface_obj [" << interface_obj->X() <<
          " " << interface_obj->Y() << " " << interface_obj->Z() << "] " <<
          " sent to ranks: [";
          for (size_t i = 0; i < partition_list.size(); ++i) {
              std::cout << partition_list[i] << " ";
          }
          std::cout << "]" << std::endl;
          partition_list.clear();
      }

      void PrepareMatching(CandidateManager& candidate_manager, int* local_comm_list,
                           int* local_memory_size_array) override {
          for (int partition_index = 0; partition_index < m_comm_size; ++partition_index) {
              if (candidate_manager.candidate_send_objects.count(partition_index) > 0) {
                  candidate_manager.matching_information.reserve(candidate_manager.candidate_send_objects.at(partition_index).size());
                  for (auto interface_obj : candidate_manager.candidate_send_objects.at(partition_index)) {
                      if (interface_obj->HasNeighborInPartition(partition_index)) {
                          m_send_objects[partition_index].push_back(interface_obj);
                          candidate_manager.matching_information[partition_index].push_back(1);
                          local_comm_list[partition_index] = 1;
                          ++local_memory_size_array[partition_index];
                      } else {
                          candidate_manager.matching_information[partition_index].push_back(0);
                      }
                  }
              }
          }
      }

      void FillSendBufferWithMatchInformation(CandidateManager& candidate_manager,
                                              int* send_buffer, int& send_buffer_size,
                                              const int comm_partner) override {
          int i = 0;
          if (candidate_manager.matching_information.count(comm_partner) > 0) {
              for (auto info : candidate_manager.matching_information.at(comm_partner)) {
                  send_buffer[i] = info;
                  ++i;
              }
          }
          send_buffer_size = i;
      }

      void FillBufferLocalSearch(CandidateManager& candidate_manager,
                                 InterfaceObjectConfigure::ContainerType& send_objects,
                                 int& num_objects) override {
          int i = 0;
          if (candidate_manager.candidate_send_objects.count(m_comm_rank) > 0) {
              for (auto interface_obj : candidate_manager.candidate_send_objects.at(m_comm_rank)) {
                  send_objects[i] = interface_obj;
                  ++i;
              }
          }
          num_objects = i;
      }

      void FillSendBufferRemoteSearch(CandidateManager& candidate_manager,
                                      double* send_buffer, int& send_buffer_size,
                                      const int comm_partner) override {
          int i = 0;
          if (candidate_manager.candidate_send_objects.count(comm_partner) > 0) {
              for (auto interface_obj : candidate_manager.candidate_send_objects.at(comm_partner)){
                  send_buffer[(i*3) + 0] = interface_obj->X();
                  send_buffer[(i*3) + 1] = interface_obj->Y();
                  send_buffer[(i*3) + 2] = interface_obj->Z();
                  ++i;
              }
          }
          send_buffer_size = 3 * i;
      }

      void PostProcessReceivedResults(CandidateManager& candidate_manager,
                                      const std::vector<double>& distances,
                                      const int comm_partner) override {
          int i = 0;
          if (candidate_manager.candidate_send_objects.count(comm_partner) > 0) {
              for (auto interface_obj : candidate_manager.candidate_send_objects.at(comm_partner)) {
                  if (distances[i] > -0.5f) // failed search has value "-1"
                      interface_obj->ProcessDistance(distances[i], comm_partner);
                  ++i;
              }
          }
      }

      void PostProcessReceivedResults(CandidateManager& candidate_manager,
                                      const double* distances,
                                      const int comm_partner) override {
          int i = 0;
          if (candidate_manager.candidate_send_objects.count(comm_partner) > 0) {
              for (auto interface_obj : candidate_manager.candidate_send_objects.at(comm_partner)) {
                  if (distances[i] > -0.5f) // failed search has value "-1"
                      interface_obj->ProcessDistance(distances[i], comm_partner);
                  ++i;
              }
          }
      }

      // **********************************************************************
      // Side where we search neighbors aka origin ****************************
      // **********************************************************************
      void ProcessReceiveBuffer(InterfaceObjectConfigure::ContainerType& remote_p_interface_object_list,
                                const double* coordinate_list, const int coordinate_list_size,
                                int& num_objects) override {
          num_objects = coordinate_list_size / 3;

          for (int i = 0; i < num_objects; ++i) { // create InterfaceObjects
              remote_p_interface_object_list[i] = InterfaceObject::Pointer(new InterfaceObject(
                  coordinate_list[(i*3)+0], coordinate_list[(i*3)+1], coordinate_list[(i*3)+2]));
          }
      }

      void FillSendBufferWithResults(double* send_buffer, const int send_buffer_size,
                                     const std::vector<double>& min_distances) override {
          for (int i = 0; i < send_buffer_size; ++i)
              send_buffer[i] = min_distances[i];
      }

      void StoreTempSearchResults(CandidateManager& candidate_manager,
                                  std::vector<InterfaceObject::Pointer> temp_closest_results,
                                  std::vector<std::vector<double>> temp_shape_functions,
                                  std::vector<array_1d<double,2>> temp_local_coordinates,
                                  const int comm_partner) override {
          MapInsertElement(candidate_manager.candidate_receive_objects, comm_partner, temp_closest_results);
          MapInsertElement(candidate_manager.candidate_local_coordinates, comm_partner, temp_local_coordinates);
          MapInsertElement(candidate_manager.candidate_shape_functions, comm_partner, temp_shape_functions);
      }

      void ProcessMatchInformation(CandidateManager& candidate_manager,
                                   int* buffer, const int buffer_size,
                                   const int comm_partner) override {
          for (int i = 0; i < buffer_size; ++i) {
              if (buffer[i] == 1) { // Match
                  if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
                      if (!candidate_manager.candidate_receive_objects.at(comm_partner)[i]) {
                          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                                       << "\"ProcessMatchInformation\": interface_obj pointer mismatch"
                                       << std::endl;
                      }
                  }
                  m_receive_objects[comm_partner].push_back(candidate_manager.candidate_receive_objects.at(comm_partner)[i]);
                  m_shape_functions[comm_partner].push_back(candidate_manager.candidate_shape_functions.at(comm_partner)[i]);
                  m_local_coordinates[comm_partner].push_back(candidate_manager.candidate_local_coordinates.at(comm_partner)[i]);
              }
          }
      }

      // **********************************************************************
      // Functions for Mapping ************************************************
      // **********************************************************************
      void ComputeBufferSizesAndCommunicationGraph(int& send_buffer_size,
                                                   int& receive_buffer_size,
                                                   GraphType& colored_graph,
                                                   int& max_colors) override {
          int* local_comm_list = new int[m_comm_size]();
          int* local_memory_size_array = new int[m_comm_size]();

          for (auto& interface_obj : m_interface_objects) {
              if (interface_obj->NeighborFound()) { // check if the interface object already found a neighbor
                  int neighbor_rank = interface_obj->GetNeighborRank();
                  local_comm_list[neighbor_rank] = 1;
                  ++local_memory_size_array[neighbor_rank];
              }
          }

          // sizes are switched bcs transfer direction is inverted from searching to mapping
          MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                                                    receive_buffer_size,
                                                    send_buffer_size,
                                                    m_comm_rank,
                                                    m_comm_size);

          MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, m_comm_size,
                                                   colored_graph, max_colors);

          delete [] local_comm_list;
          delete [] local_memory_size_array;
      }
      void FillBufferWithValues(double* buffer, int& buffer_size, const int comm_partner,
                                const Variable<double>& variable, Kratos::Flags& options) override {
          int i = 0;
          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_receive_objects.count(comm_partner) > 0) {
              interface_objects = m_receive_objects.at(comm_partner);
          }

          for (auto interface_obj : interface_objects) {
              if (options.Is(MapperFlags::INTERPOLATE_VALUES)) {
                  buffer[i] = interface_obj->GetObjectValueInterpolated(variable, m_shape_functions.at(m_comm_rank)[i]);
              } else {
                  buffer[i] = interface_obj->GetObjectValue(variable);
              }
              ++i;
          }

          buffer_size = static_cast<int>(interface_objects.size());

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (buffer_size != i) {
                  KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                               << "\"FillSendBufferWithValues, double\": size mismatch"
                               << std::endl;
              }
          }
      }

      void FillBufferWithValues(double* buffer, int& buffer_size, const int comm_partner,
                                const Variable< array_1d<double,3> >& variable,
                                Kratos::Flags& options) override {
          int i = 0;
          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_receive_objects.count(comm_partner) > 0) {
              interface_objects = m_receive_objects.at(comm_partner);
          }

          for (auto interface_obj : interface_objects) {
              if (options.Is(MapperFlags::INTERPOLATE_VALUES)) {
                  buffer[(i*3) + 0] = interface_obj->GetObjectValueInterpolated(variable, m_shape_functions.at(m_comm_rank)[i])[0];
                  buffer[(i*3) + 1] = interface_obj->GetObjectValueInterpolated(variable, m_shape_functions.at(m_comm_rank)[i])[1];
                  buffer[(i*3) + 2] = interface_obj->GetObjectValueInterpolated(variable, m_shape_functions.at(m_comm_rank)[i])[2];
              } else {
                  buffer[(i*3) + 0] = interface_obj->GetObjectValue(variable)[0];
                  buffer[(i*3) + 1] = interface_obj->GetObjectValue(variable)[1];
                  buffer[(i*3) + 2] = interface_obj->GetObjectValue(variable)[2];
              }
              ++i;
          }

          buffer_size = static_cast<int>(interface_objects.size()) * 3;

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (buffer_size != i*3) {
                  KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                               << "\"FillSendBufferWithValues, double<3>\": size mismatch"
                               << std::endl;
              }
          }
      }

      void ProcessValues(const double* buffer, const int buffer_size, const int comm_partner,
                         const Variable<double>& variable,
                         Kratos::Flags& options, const double factor) override {
          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_send_objects.count(comm_partner) > 0) {
             interface_objects = m_send_objects.at(comm_partner);
          }

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (static_cast<int>(interface_objects.size()) != buffer_size) {
                  KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                               << "\"ProcessValues, double\": Wrong number of results "
                               << "received!; interface_objects.size() = " << interface_objects.size()
                               << ", buffer_size = " << buffer_size << std::endl;
              }
          }

          for (int i = 0; i < buffer_size; ++i) {
              interface_objects[i]->SetObjectValue(variable, buffer[i],
                                                   options, factor);
          }
      }

      void ProcessValues(const double* buffer, const int buffer_size, const int comm_partner,
                         const Variable< array_1d<double,3> >& variable,
                         Kratos::Flags& options, const double factor) override {

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (buffer_size % 3 != 0) {
                  KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                               << "\"ProcessValues, double<3>\": Uneven number of results "
                               << "received!; buffer_size modulo 3 = "
                               << buffer_size % 3 << std::endl;
              }
          }

          const int num_values = buffer_size / 3;

          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_send_objects.count(comm_partner) > 0) {
             interface_objects = m_send_objects.at(comm_partner);
          }

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (static_cast<int>(interface_objects.size()) != num_values) {
                 KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerParallel; "
                              << "\"ProcessValues, double<3>\": Wrong number of results "
                              << "received!; interface_objects.size() = "
                              << interface_objects.size() << ", num_values = "
                              << num_values << std::endl;
              }
          }

          array_1d<double,3> value;

          for (int i = 0; i < num_values; ++i) {
              value[0] = buffer[(i*3) + 0];
              value[1] = buffer[(i*3) + 1];
              value[2] = buffer[(i*3) + 2];

              interface_objects[i]->SetObjectValue(variable, value,
                                                   options, factor);
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
        buffer << "InterfaceObjectManagerParallel" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceObjectManagerParallel";}

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
      double* m_bounding_box = new double[6] {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
                                                // xmax, xmin,  ymax, ymin,  zmax, zmin


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
      InterfaceObjectManagerParallel& operator=(InterfaceObjectManagerParallel const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerParallel(InterfaceObjectManagerParallel const& rOther){}


      ///@}

    }; // Class InterfaceObjectManagerParallel

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceObjectManagerParallel& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceObjectManagerParallel& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_PARALLEL_H_INCLUDED  defined
