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

#if !defined(KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "mapper_communicator.h"
#include "interface_search_structure_mpi.h"
#include "mpi_manager.h"


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
  class MapperMPICommunicator : public MapperCommunicator
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperMPICommunicator
      KRATOS_CLASS_POINTER_DEFINITION(MapperMPICommunicator);

      ///@}
      ///@name Life Cycle
      ///@{

      MapperMPICommunicator(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination) :
            MapperCommunicator(i_model_part_origin, i_model_part_destination) {
      }


      /// Destructor.
      virtual ~MapperMPICommunicator(){}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Initialize(InterfaceObjectManager::Pointer i_interface_object_manager_origin,
                      InterfaceObjectManager::Pointer i_interface_object_manager_destination,
                      bool bi_directional_search) override {

          m_interface_object_manager_origin = i_interface_object_manager_origin;
          m_interface_object_manager_destination = i_interface_object_manager_destination;
          m_bidirectional_search = bi_directional_search;

          m_search_structure_mpi_map = InterfaceSearchStructureMPI::Pointer ( new InterfaceSearchStructureMPI(
             i_interface_object_manager_destination, i_interface_object_manager_origin,
             m_model_part_origin, MyPID(), TotalProcesses()) );

          if (bi_directional_search) {
              m_search_structure_mpi_inverse_map = InterfaceSearchStructureMPI::Pointer ( new InterfaceSearchStructureMPI(
                 i_interface_object_manager_origin, i_interface_object_manager_destination,
                 m_model_part_destination, MyPID(), TotalProcesses()) );
          }
          KRATOS_WATCH("End of InitializeMPI")

          ComputeSearchStructure();
      }

      void Initialize(InterfaceObjectManager::Pointer i_interface_object_manager_origin,
                      InterfaceObjectManager::Pointer i_interface_object_manager_destination,
                      bool bi_directional_search, double i_initial_search_radius,
                      int i_max_search_iterations) override {

          m_interface_object_manager_origin = i_interface_object_manager_origin;
          m_interface_object_manager_destination = i_interface_object_manager_destination;
          m_bidirectional_search = bi_directional_search;

          m_search_structure_mpi_map = InterfaceSearchStructureMPI::Pointer ( new InterfaceSearchStructureMPI(
             i_interface_object_manager_destination, i_interface_object_manager_origin, m_model_part_origin,
             i_initial_search_radius, i_max_search_iterations, MyPID(), TotalProcesses()) );

          if (bi_directional_search) {
              m_search_structure_mpi_inverse_map = InterfaceSearchStructureMPI::Pointer ( new InterfaceSearchStructureMPI(
                 i_interface_object_manager_origin, i_interface_object_manager_destination, m_model_part_destination,
                 i_initial_search_radius, i_max_search_iterations, MyPID(), TotalProcesses()) );
          }
          KRATOS_WATCH("End of InitializeMPI with Input(SearchRadius and MaxIter)")

          ComputeSearchStructure();
      }

      void ComputeSearchStructure() override {
          m_search_structure_mpi_map->Clear();
          m_search_structure_mpi_map->Search();
          m_search_structure_mpi_map->GetMPIData(m_max_send_buffer_size_map,
                                             m_max_receive_buffer_size_map,
                                             m_colored_graph_map,
                                             m_max_colors_map);
          if (m_bidirectional_search)
              m_search_structure_mpi_inverse_map ->Clear();
              m_search_structure_mpi_inverse_map ->Search();
              m_search_structure_mpi_inverse_map->GetMPIData(m_max_send_buffer_size_inverse_map,
                                                 m_max_receive_buffer_size_inverse_map,
                                                 m_colored_graph_inverse_map,
                                                 m_max_colors_inverse_map);
      }

      void TransferData(bool direction, const Variable<double>& origin_variable,
                        const Variable<double>& destination_variable,
                        const bool add_value,
                        bool use_other_search_structure = false) override {

          TransferDataParallel(origin_variable, destination_variable,
                               direction, add_value, use_other_search_structure);
      }

      void TransferData(bool direction, const Variable< array_1d<double,3> >& origin_variable,
                        const Variable< array_1d<double,3> >& destination_variable,
                        const bool add_value,
                        bool use_other_search_structure = false) override {

          TransferDataParallel(origin_variable, destination_variable,
                               direction, add_value, use_other_search_structure);
      }

      int MyPID () override { // Copy from "kratos/includes/mpi_communicator.h"
          int rank;
          MPI_Comm_rank(MPI_COMM_WORLD, &rank);
          return rank;
      }

      int TotalProcesses() override { // Copy from "kratos/includes/mpi_communicator.h"
          int nproc;
          MPI_Comm_size(MPI_COMM_WORLD, &nproc);
          return nproc;
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
        buffer << "MapperMPICommunicator" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperMPICommunicator";}

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

      InterfaceSearchStructureMPI::Pointer m_search_structure_mpi_map;
      InterfaceSearchStructureMPI::Pointer m_search_structure_mpi_inverse_map;

      int m_max_send_buffer_size_map;
      int m_max_send_buffer_size_inverse_map;

      int m_max_receive_buffer_size_map;
      int m_max_receive_buffer_size_inverse_map;

      GraphType m_colored_graph_map;
      GraphType m_colored_graph_inverse_map;

      int m_max_colors_map;
      int m_max_colors_inverse_map;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      template <typename T>
      void TransferDataParallel(const Variable< T >& origin_variable,
                               const Variable< T >& destination_variable,
                               const bool direction, const bool add_value,
                               const bool use_other_search_structure) {
          if ((!m_bidirectional_search && !direction && !use_other_search_structure) ||
              (!m_bidirectional_search && direction && use_other_search_structure))
              KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; \"TransferDataParallel\" Bi-Directional Search not defined!" << std::endl;

          // TODO here some implicit conversion is done, is that ok?
          std::vector<InterfaceObject::Pointer> local_mapping_points;
          std::vector<InterfaceObject::Pointer> remote_mapping_points;

          SelectPointLists(remote_mapping_points, local_mapping_points,
                           direction, use_other_search_structure);

          ExchangeDataLocal(origin_variable, destination_variable, remote_mapping_points,
                            local_mapping_points, add_value);

          ExchangeDataRemote(origin_variable, destination_variable,
                             m_interface_object_manager_origin, m_interface_object_manager_destination,
                             direction, add_value, use_other_search_structure);

          MPI_Barrier(MPI_COMM_WORLD);
      }

      void SelectPointLists(std::vector<InterfaceObject::Pointer>& point_list_1,
                            std::vector<InterfaceObject::Pointer>& point_list_2,
                            const bool direction, const bool use_other_search_structure) override {
                                // TODO unify it in terms of naming and use fnct of base-class!
          // boolean "use_other_search_structure" is intended to be used for mortar mapping
          // to exchange data in both directions with the same match-partners
          if (direction) { // map; origin -> destination
              if (use_other_search_structure) { // use search structure of InverseMap
                  point_list_1 = m_interface_object_manager_origin->GetLocalMappingVector();
                  point_list_2 = m_interface_object_manager_destination->GetRemoteMappingVector();
              } else { // use search structure of Map
                  point_list_1 = m_interface_object_manager_origin->GetRemoteMappingVector();
                  point_list_2 = m_interface_object_manager_destination->GetLocalMappingVector();
              }
          } else { // inverse map; destination -> origin
              if (use_other_search_structure) { // use search structure of Map
                  point_list_1 = m_interface_object_manager_destination->GetLocalMappingVector();
                  point_list_2 = m_interface_object_manager_origin->GetRemoteMappingVector();
              } else { // use search structure of InverseMap
                  point_list_1 = m_interface_object_manager_destination->GetRemoteMappingVector();
                  point_list_2 = m_interface_object_manager_origin->GetLocalMappingVector();
              }
          }
          if (point_list_1.size() != point_list_2.size())
              KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; \"SelectPointLists\" Size Mismatch!" << std::endl;
      }

      template <typename T>
      void ExchangeDataRemote(const Variable< T >& origin_variable,
                              const Variable< T >& destination_variable,
                              InterfaceObjectManager::Pointer& object_manager_1,
                              InterfaceObjectManager::Pointer& object_manager_2,
                              const bool direction, const bool add_value,
                              const bool use_other_search_structure) {
          MPI_Status status; // TODO use this for some checks, e.g. if message has arrived,...
          int send_buffer_size = 0;
          int receive_buffer_size = 0;

          InterfaceObjectManager::Pointer interface_object_manager_1;
          InterfaceObjectManager::Pointer interface_object_manager_2;

          GraphType colored_graph;
          int max_colors;
          int max_send_buffer_size;
          int max_receive_buffer_size;

          if (direction) { // map; origin -> destination
              interface_object_manager_1 = object_manager_1;
              interface_object_manager_2 = object_manager_2;
              colored_graph = m_colored_graph_map;
              max_colors = m_max_colors_map;
              max_send_buffer_size = m_max_send_buffer_size_map;
              max_receive_buffer_size = m_max_receive_buffer_size_map;
          }
          else { // inverse map; destination -> origin
              interface_object_manager_1 = object_manager_2;
              interface_object_manager_2 = object_manager_1;
              colored_graph = m_colored_graph_inverse_map;
              max_colors = m_max_colors_inverse_map;
              max_send_buffer_size = m_max_send_buffer_size_inverse_map;
              max_receive_buffer_size = m_max_receive_buffer_size_inverse_map;
          }

          int buffer_size_factor = MPIManager::SizeOfVariable(T());

          max_send_buffer_size *= buffer_size_factor;
          max_receive_buffer_size *= buffer_size_factor;

          double* send_buffer = new double[max_send_buffer_size];
          double* receive_buffer = new double[max_receive_buffer_size];

          for (int i = 0; i < max_colors; ++i) { // loop over communication steps (aka. colour)
              int comm_partner = colored_graph(MyPID(), i); // get the partner rank
              if (comm_partner != -1) { // check if rank is communicating in this communication step (aka. colour)
                  interface_object_manager_1->FillSendBufferWithValues(send_buffer, send_buffer_size, comm_partner,
                                                                       origin_variable, use_other_search_structure);

                  MPIManager::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size, receive_buffer_size,
                              max_send_buffer_size, max_receive_buffer_size, comm_partner, status);

                  interface_object_manager_2->ProcessValues(receive_buffer, receive_buffer_size, comm_partner, destination_variable,
                                                            add_value, use_other_search_structure);
              } // if I am communicating in this loop (comm_partner != -1)
          } // loop colors

          delete [] send_buffer;
          delete [] receive_buffer;
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
      MapperMPICommunicator& operator=(MapperMPICommunicator const& rOther);

    //   /// Copy constructor.
    //   MapperMPICommunicator(MapperMPICommunicator const& rOther){}


      ///@}

    }; // Class MapperMPICommunicator

  ///@}


  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperMPICommunicator& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperMPICommunicator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED  defined
