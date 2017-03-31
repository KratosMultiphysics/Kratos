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

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map> // for CandidateManager


// External includes


// Project includes
#include "includes/define.h"

#include "interface_node.h"
#include "interface_condition.h"
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

  // This struct stores temporary variables, such that they don't have to be
  // stored as members in the InterfaceObjectManager
  struct CandidateManager {
      std::unordered_map<int, std::vector<InterfaceObject::Pointer> > candidate_send_objects;
      std::unordered_map<int, std::vector<InterfaceObject::Pointer> > candidate_receive_objects;
      std::unordered_map<int, std::vector<std::vector<double> > > candidate_shape_functions;
      std::unordered_map<int, std::vector<array_1d<double,2> > > candidate_local_coordinates;

      std::unordered_map<int, std::vector<int> > matching_information;
  };

  /// Short class definition.
  /** Detail class definition.
  */
  class InterfaceObjectManagerBase
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceObjectManagerBase
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerBase);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Destructor.
      virtual ~InterfaceObjectManagerBase() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Clear() {
          m_send_objects.clear();
          m_receive_objects.clear();
          m_shape_functions.clear();
          m_local_coordinates.clear();

          for (auto& interface_obj : m_interface_objects) {
              interface_obj->Reset(); // Set the object back to its initial state
              // i.e. reset its information whether it has been sent somewhere or
              // whether it already found a neighbor
          }
      }

      template <typename T>
      void MapInsertElement(std::unordered_map<int, T>& map, int key, T& value) {
          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (map.count(key) > 0) {
                  KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                               << "\"InsertElement\" Key already present in Map!"
                               << std::endl;
              }
          }
          map.emplace(key, value);
      }

      // **********************************************************************
      // Side we want to find neighbors for aka destination *******************
      // **********************************************************************
      // ** InterfaceObjectManagerSerial and InterfaceObjectManagerParallel **

      bool AllNeighborsFound() {
          int all_neighbors_found = 1; // set to "1" aka "true" by default in case
                                       // this partition doesn't have a part of the interface!

          for (auto& interface_obj : m_interface_objects) {
              if (!interface_obj->NeighborFound()) {
                  all_neighbors_found = 0;
              }
          }

          // This is necessary bcs not all partitions would start a new search iteration!
          m_model_part.GetCommunicator().MinAll(all_neighbors_found);

          return all_neighbors_found;
      }

      virtual void CheckResults() = 0;

      InterfaceObjectConfigure::ContainerType& GetInterfaceObjects() {
          return m_interface_objects;
      }

      // ***** InterfaceObjectManagerSerial *****
      virtual void GetInterfaceObjectsSerialSearch(InterfaceObjectConfigure::ContainerType& candidate_send_objects) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"GetInterfaceObjectsSerialSearch\" of the base class called!"
                       << std::endl;
      }

      virtual void PostProcessReceivedResults(const std::vector<double>& distances,
                                              const InterfaceObjectConfigure::ContainerType& candidate_send_objects) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"PostProcessReceivedResults\" of the base "
                       << "class called!" << std::endl;
      }

      // ***** InterfaceObjectManagerParallel *****
      virtual void ComputeCandidatePartitions(CandidateManager& candidate_manager, int* local_comm_list,
                                              int* local_memory_size_array,
                                              double* global_bounding_boxes,
                                              const bool last_iteration) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ComputeCandidatePartitions\" of the base "
                       << "class called!" << std::endl;
      }

      // Function for debugging
      virtual void PrintCandidatePartitions(const InterfaceObject::Pointer point,
                                            std::vector<int>& partition_list) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"PrintCandidatePartitions\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void PrepareMatching(CandidateManager& candidate_manager, int* local_comm_list,
                                   int* local_memory_size_array) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"PrepareMatching\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void FillSendBufferWithMatchInformation(CandidateManager& candidate_manager,
                                                      int* send_buffer, int& send_buffer_size,
                                                      const int comm_partner){
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"FillSendBufferWithMatchInformation\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void FillBufferLocalSearch(CandidateManager& candidate_manager,
                                         InterfaceObjectConfigure::ContainerType& send_points,
                                         int& num_objects) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"FillBufferLocalSearch\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void FillSendBufferRemoteSearch(CandidateManager& candidate_manager,
                                              double* send_buffer, int& send_buffer_size,
                                              const int comm_partner) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"FillSendBufferRemoteSearch\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void PostProcessReceivedResults(CandidateManager& candidate_manager,
                                              const std::vector<double>& distances,
                                              const int comm_partner) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"PostProcessReceivedResults, std::vector<double> \" "
                       << "of the base class called!" << std::endl;
      }

      virtual void PostProcessReceivedResults(CandidateManager& candidate_manager,
                                              const double* distances,
                                              const int comm_partner) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"PostProcessReceivedResults, C-array\" of the base "
                       << "class called!" << std::endl;
      }

      // **********************************************************************
      // Side where we search neighbors aka origin ****************************
      // **********************************************************************
      // ***** InterfaceObjectManagerSerial *****
      virtual void StoreSearchResults(const std::vector<double>& distances,
                                      const std::vector<InterfaceObject::Pointer> temp_closest_results,
                                      const std::vector<std::vector<double>> temp_shape_functions,
                                      const std::vector<array_1d<double,2>> temp_local_coordinates) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"StoreSearchResults\" of the base "
                       << "class called!" << std::endl;
      }

      virtual std::vector<InterfaceObject::Pointer>& GetDestinationInterfaceObjects() {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"GetDestinationInterfaceObjects\" of the base "
                       << "class called!" << std::endl;
      }

      virtual std::vector<InterfaceObject::Pointer>& GetOriginInterfaceObjects() {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"GetOriginInterfaceObjects\" of the base "
                       << "class called!" << std::endl;
      }

      // ***** InterfaceObjectManagerParallel *****
      virtual void ProcessReceiveBuffer(InterfaceObjectConfigure::ContainerType& remote_p_point_list,
                                const double* coordinate_list, const int coordinate_list_size,
                                int& num_objects){
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessReceiveBuffer\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void FillSendBufferWithResults(double* send_buffer, const int send_buffer_size,
                                             const std::vector<double>& min_distances){
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"FillSendBufferWithResults\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void StoreTempSearchResults(CandidateManager& candidate_manager,
                                          const std::vector<InterfaceObject::Pointer> temp_closest_results,
                                          const std::vector<std::vector<double>> temp_shape_functions,
                                          const std::vector<array_1d<double,2>> temp_local_coordinates,
                                          const int comm_partner) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"StoreTempSearchResults\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void ProcessMatchInformation(CandidateManager& candidate_manager,
                                           int* buffer, const int buffer_size,
                                           const int comm_partner){
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessMatchInformation\" of the base "
                       << "class called!" << std::endl;
      }

      // **********************************************************************
      // Functions for Mapping ************************************************
      // **********************************************************************
      // ** InterfaceObjectManagerSerial and InterfaceObjectManagerParallel **
      template <typename T>
      void FillBufferWithValues(std::vector< T >& buffer, const Variable< T >& variable,
                                Kratos::Flags& options) {
          int i = 0;

          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_receive_objects.count(m_comm_rank) > 0) {
              interface_objects = m_receive_objects.at(m_comm_rank);
          }

          buffer.resize(interface_objects.size());

          for (auto interface_obj : interface_objects) {
              if (options.Is(MapperFlags::INTERPOLATE_VALUES)) {
                  buffer[i] = interface_obj->GetObjectValueInterpolated(variable, m_shape_functions.at(m_comm_rank)[i]);
              } else {
                  buffer[i] = interface_obj->GetObjectValue(variable);
              }
              ++i;
          }
      }

      template <typename T>
      void ProcessValues(const std::vector< T >& buffer, const Variable< T >& variable,
                         Kratos::Flags& options, const double factor) {
          std::vector<InterfaceObject::Pointer> interface_objects;
          if (m_send_objects.count(m_comm_rank) > 0) {
             interface_objects = m_send_objects.at(m_comm_rank);
          }

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
             if (interface_objects.size() != buffer.size()) {
                 KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                              << "\"ProcessValues, T\": Wrong number of results "
                              << "received!;\" interface_objects.size() = "
                              << interface_objects.size() << ", buffer.size() = "
                              << buffer.size() << std::endl;
             }
          }

          for (std::size_t i = 0; i < interface_objects.size(); ++i) {
             interface_objects[i]->SetObjectValue(variable, buffer[i],
                                                  options, factor);
          }
      }

      // ***** InterfaceObjectManagerParallel *****
      virtual void ComputeBufferSizesAndCommunicationGraph(int& max_send_buffer_size,
                                                           int& max_receive_buffer_size,
                                                           GraphType& colored_graph,
                                                           int& max_colors) {
          // TODO error
      }
      virtual void FillBufferWithValues(double* buffer, int& buffer_size, const int comm_partner,
                                        const Variable<double>& variable, Kratos::Flags& options) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessMatchInformation\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void FillBufferWithValues(double* buffer, int& buffer_size, const int comm_partner,
                                        const Variable< array_1d<double,3> >& variable,
                                        Kratos::Flags& options) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessMatchInformation\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void ProcessValues(const double* buffer, const int buffer_size, const int comm_partner,
                         const Variable<double>& variable,
                         Kratos::Flags& options, const double factor) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessMatchInformation\" of the base "
                       << "class called!" << std::endl;
      }

      virtual void ProcessValues(const double* buffer, const int buffer_size, const int comm_partner,
                         const Variable< array_1d<double,3> >& variable,
                         Kratos::Flags& options, const double factor) {
          KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; "
                       << "\"ProcessMatchInformation\" of the base "
                       << "class called!" << std::endl;
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
        buffer << "InterfaceObjectManagerBase" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceObjectManagerBase";}

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

      InterfaceObjectManagerBase(ModelPart& i_model_part, int i_comm_rank, int i_comm_size,
                                 MapperUtilities::InterfaceObjectConstructionType i_interface_object_type,
                                 GeometryData::IntegrationMethod i_integration_method, int i_echo_level) :
                                 m_model_part(i_model_part) {

          m_comm_rank = i_comm_rank;
          m_comm_size = i_comm_size;

          m_echo_level = i_echo_level;

          if (i_interface_object_type == MapperUtilities::Node) {
              InitializeInterfaceNodeManager(i_model_part, i_comm_rank, i_comm_size);
          } else if (i_interface_object_type == MapperUtilities::Condition_Center ||
                     i_interface_object_type == MapperUtilities::Condition_Gauss_Point) {
              InitializeInterfaceConditionManager(i_model_part, i_comm_rank, i_comm_size, i_integration_method);
          } else {
              KRATOS_ERROR << "MappingApplication; InterfaceObjectManagerBase; \"CreateInterfaceObjectManager\" "
                           << "type of interface object construction not implemented" << std::endl;
          }
      }

      ModelPart& m_model_part;

      InterfaceObjectConfigure::ContainerType m_interface_objects;

      int m_comm_rank = 0;
      int m_comm_size = 0;
      int m_echo_level = 0;

      // point-sending interface (destination)
      std::unordered_map<int, std::vector<InterfaceObject::Pointer> > m_send_objects;

      // point-receiving interface
      std::unordered_map<int, std::vector<InterfaceObject::Pointer> > m_receive_objects;
      std::unordered_map<int, std::vector<std::vector<double> > > m_shape_functions;
      std::unordered_map<int, std::vector<array_1d<double,2> > > m_local_coordinates;


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

      void InitializeInterfaceNodeManager(ModelPart& model_part, int comm_rank, int comm_size) {
          m_interface_objects.resize(model_part.GetCommunicator().LocalMesh().NumberOfNodes());

          int i = 0;
          for (auto &local_node : model_part.GetCommunicator().LocalMesh().Nodes()) {
              m_interface_objects[i] = InterfaceObject::Pointer( new InterfaceNode(local_node) );
              ++i;
          }
      }

      void InitializeInterfaceConditionManager(ModelPart& model_part, int comm_rank, int comm_size,
                                               GeometryData::IntegrationMethod integration_method) {
          bool construct_with_center;
          int size_factor = 1;
          if (integration_method == GeometryData::NumberOfIntegrationMethods) {
              construct_with_center = true;
              size_factor = 1;
          } else {
              construct_with_center = false;
              if (model_part.GetCommunicator().LocalMesh().NumberOfConditions() > 0) {
                  size_factor = model_part.GetCommunicator().LocalMesh().ConditionsBegin()->GetGeometry().IntegrationPointsNumber(integration_method);
              }
          }

          m_interface_objects.reserve(size_factor * model_part.GetCommunicator().LocalMesh().NumberOfConditions());

          if (construct_with_center) { // construct with condition center point
              for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
                  m_interface_objects.push_back(InterfaceObject::Pointer( new InterfaceCondition(condition, condition.GetGeometry().Center()) ));
              }
          } else { // construct with condition gauss points
              for (auto& condition : model_part.GetCommunicator().LocalMesh().Conditions()) {
                  const Geometry< Node<3> >& condition_geometry = condition.GetGeometry();

                  Matrix shape_functions = condition_geometry.ShapeFunctionsValues(integration_method); // TODO "ShapeFunctionsValues" seems to not be implemented for every geometry!!!

                  const int num_gauss_points = shape_functions.size1();
                  const int num_nodes = shape_functions.size2();

                  array_1d<double, 3> gauss_point_global_coords;

                  for (int g = 0; g < num_gauss_points; ++g) {
                      gauss_point_global_coords[0] = 0.0f;
                      gauss_point_global_coords[1] = 0.0f;
                      gauss_point_global_coords[2] = 0.0f;

                      for (int n = 0; n < num_nodes; ++n) {
                          gauss_point_global_coords[0] += shape_functions(g, n) * condition_geometry[n].X();
                          gauss_point_global_coords[1] += shape_functions(g, n) * condition_geometry[n].Y();
                          gauss_point_global_coords[2] += shape_functions(g, n) * condition_geometry[n].Z();
                      }
                      // TODO check again if this is whole computation of the GPs is correct

                      m_interface_objects.push_back(InterfaceObject::Pointer( new InterfaceCondition(condition, gauss_point_global_coords) ));
                  }
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
      InterfaceObjectManagerBase& operator=(InterfaceObjectManagerBase const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerBase(InterfaceObjectManagerBase const& rOther){}


      ///@}

    }; // Class InterfaceObjectManagerBase

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceObjectManagerBase& rThis)
  {
        return rIStream;
  }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceObjectManagerBase& rThis)
  {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
  }
  ///@}

  ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_BASE_H_INCLUDED  defined
