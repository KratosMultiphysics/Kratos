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
#include "interface_object_manager_parallel.h"
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

      MapperMPICommunicator(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                            Parameters& i_json_parameters) :
            MapperCommunicator(i_model_part_origin, i_model_part_destination, i_json_parameters) { }

      /// Destructor.
      virtual ~MapperMPICommunicator() { }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void InitializeOrigin(MapperUtilities::InterfaceObjectConstructionType i_interface_object_type_origin,
                      GeometryData::IntegrationMethod i_integration_method_origin = GeometryData::NumberOfIntegrationMethods) override {

          m_p_interface_object_manager_origin = InterfaceObjectManagerParallel::Pointer(
              new InterfaceObjectManagerParallel(m_model_part_origin, MyPID(),
                                                 TotalProcesses(),
                                                 i_interface_object_type_origin,
                                                 i_integration_method_origin,
                                                 m_echo_level) );

          m_interface_object_type_origin = i_interface_object_type_origin;
          m_integration_method_origin = i_integration_method_origin;
      }

      void InitializeDestination(MapperUtilities::InterfaceObjectConstructionType i_interface_object_type_destination,
                      GeometryData::IntegrationMethod i_integration_method_destination = GeometryData::NumberOfIntegrationMethods) override {

          m_p_interface_object_manager_destination = InterfaceObjectManagerParallel::Pointer(
              new InterfaceObjectManagerParallel(m_model_part_destination, MyPID(),
                                                 TotalProcesses(),
                                                 i_interface_object_type_destination,
                                                 i_integration_method_destination,
                                                 m_echo_level) );

          m_interface_object_type_destination = i_interface_object_type_destination;
          m_integration_method_destination = i_integration_method_destination;
      }



      void TransferNodalData(const Variable<double>& origin_variable,
                             const Variable<double>& destination_variable,
                             Kratos::Flags& options,
                             double factor = 1.0f) override {
          TransferDataParallel(origin_variable, destination_variable, options, factor);
      }

      void TransferNodalData(const Variable< array_1d<double,3> >& origin_variable,
                             const Variable< array_1d<double,3> >& destination_variable,
                             Kratos::Flags& options,
                             double factor = 1.0f) override {
          TransferDataParallel(origin_variable, destination_variable, options, factor);
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

      int m_max_send_buffer_size;
      int m_max_receive_buffer_size;

      GraphType m_colored_graph;
      int m_max_colors;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void InitializeSearchStructure() override {
          m_p_search_structure = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructureMPI(
             m_p_interface_object_manager_destination, m_p_interface_object_manager_origin, m_model_part_origin,
             MyPID(), TotalProcesses(), m_echo_level) );
      }

      virtual void InvokeSearch(const double i_initial_search_radius,
                                const int i_max_search_iterations) {
          m_p_search_structure->Search(i_initial_search_radius,
                                       i_max_search_iterations);

          m_p_interface_object_manager_destination->ComputeBufferSizesAndCommunicationGraph(
                                           m_max_send_buffer_size,
                                           m_max_receive_buffer_size,
                                           m_colored_graph,
                                           m_max_colors);
      }

      template <typename T>
      void TransferDataParallel(const Variable< T >& origin_variable,
                                const Variable< T >& destination_variable,
                                Kratos::Flags& options,
                                double factor) {
          if (options.Is(MapperFlags::SWAP_SIGN)) {
              factor *= (-1);
          }

          ExchangeDataLocal(origin_variable, destination_variable,
                            options, factor);

          ExchangeDataRemote(origin_variable, destination_variable,
                             options, factor,
                             m_p_interface_object_manager_origin,
                             m_p_interface_object_manager_destination);

          MPI_Barrier(MPI_COMM_WORLD);
      }

      template <typename T>
      void ExchangeDataRemote(const Variable< T >& origin_variable,
                              const Variable< T >& destination_variable,
                              Kratos::Flags& options,
                              const double factor,
                              InterfaceObjectManagerBase::Pointer& object_manager_1,
                              InterfaceObjectManagerBase::Pointer& object_manager_2) {
          int send_buffer_size = 0;
          int receive_buffer_size = 0;

          int max_send_buffer_size= m_max_send_buffer_size;
          int max_receive_buffer_size = m_max_receive_buffer_size;

          int buffer_size_factor = MapperUtilitiesMPI::SizeOfVariable(T());

          max_send_buffer_size *= buffer_size_factor;
          max_receive_buffer_size *= buffer_size_factor;

          double* send_buffer = new double[max_send_buffer_size];
          double* receive_buffer = new double[max_receive_buffer_size];

          for (int i = 0; i < m_max_colors; ++i) { // loop over communication steps (aka. colour)
              int comm_partner = m_colored_graph(MyPID(), i); // get the partner rank
              if (comm_partner != -1) { // check if rank is communicating in this communication step (aka. colour)
                  m_p_interface_object_manager_origin->FillBufferWithValues(send_buffer, send_buffer_size, comm_partner,
                                                                            origin_variable, options);

                  MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size, receive_buffer_size,
                                                  max_send_buffer_size, max_receive_buffer_size, comm_partner);

                  m_p_interface_object_manager_destination->ProcessValues(receive_buffer, receive_buffer_size, comm_partner,
                                                                          destination_variable, options, factor);
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
