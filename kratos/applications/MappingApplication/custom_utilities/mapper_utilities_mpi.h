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

#if !defined(KRATOS_MAPPER_UTILITIES_MPI_H_INCLUDED )
#define  KRATOS_MAPPER_UTILITIES_MPI_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes
#include "mpi.h"

// Project includes
#include "includes/define.h"


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
  class MapperUtilitiesMPI
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MapperUtilitiesMPI
      KRATOS_CLASS_POINTER_DEFINITION(MapperUtilitiesMPI);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MapperUtilitiesMPI(){}

      /// Destructor.
      virtual ~MapperUtilitiesMPI(){}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      template<class T>
      static inline int SizeOfVariable(const T& variable);

      template<typename T>
      static void MpiSendRecv(T* send_buffer, T* receive_buffer, const int send_buffer_size,
                              int& receive_buffer_size, const int comm_partner) {
          // Determime data type
  		    MPI_Datatype DataType = MapperUtilitiesMPI::GetMPIDatatype(T());

          //   Exchange the information about the receiving buffer size
          MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                       1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          // Perform the actual exchange of the buffers
          MPI_Sendrecv(send_buffer, send_buffer_size, DataType, comm_partner, 0, receive_buffer,
                       receive_buffer_size, DataType, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }

      // This function checks if the buffer sizes are large enough
      template<typename T>
      static void MpiSendRecv(T* send_buffer, T* receive_buffer, const int send_buffer_size,
                              int& receive_buffer_size, const int max_send_buffer_size,
                              const int max_receive_buffer_size, const int comm_partner) {
          // Determime data type
          MPI_Datatype DataType = MapperUtilitiesMPI::GetMPIDatatype(T());

          //   Exchange the information about the receiving buffer size
          MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                       1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          if (MapperUtilities::MAPPER_DEBUG_LEVEL) {
              if (send_buffer_size > max_send_buffer_size) {
                  KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; "
                               << "\"MpiSendRecv\" Send Buffer too small!; "
                               << "send_buffer_size = " << send_buffer_size
                               << ", max_send_buffer_size = "
                               << max_send_buffer_size << std::endl;
              }

              if (receive_buffer_size > max_receive_buffer_size) {
                  KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; "
                               << "\"MpiSendRecv\" Receive Buffer too small!; "
                               << "receive_buffer_size = " << receive_buffer_size
                               << ", max_receive_buffer_size = "
                               << max_receive_buffer_size << std::endl;
              }
          }

          // Perform the actual exchange of the buffers
          MPI_Sendrecv(send_buffer, send_buffer_size, DataType, comm_partner, 0, receive_buffer,
                       receive_buffer_size, DataType, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
        buffer << "MapperUtilitiesMPI" ;
        return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MapperUtilitiesMPI";}

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
      /// An auxiliary function to determine the MPI_Datatype corresponding to a given C type
      // copied from "external_libraries/mpi_python/mpi_python.h"
      template<class T>
      static inline MPI_Datatype GetMPIDatatype(const T& value);


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
      MapperUtilitiesMPI& operator=(MapperUtilitiesMPI const& rOther);

    //   /// Copy constructor.
    //   MapperUtilitiesMPI(MapperUtilitiesMPI const& rOther){}


      ///@}

    }; // Class MapperUtilitiesMPI

  ///@}

  template<>
  inline MPI_Datatype MapperUtilitiesMPI::GetMPIDatatype<int>(const int& value) {
      return MPI_INT ;
  }

  template<>
  inline MPI_Datatype MapperUtilitiesMPI::GetMPIDatatype<double>(const double& value) {
      return MPI_DOUBLE ;
  }

  template<>
  inline int MapperUtilitiesMPI::SizeOfVariable< double >(const double& variable) {
      return 1;
  }

  template<>
  inline int MapperUtilitiesMPI::SizeOfVariable< array_1d<double,3> >(const array_1d<double,3>& variable) {
      return 3;
  }

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MapperUtilitiesMPI& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MapperUtilitiesMPI& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_UTILITIES_MPI_H_INCLUDED  defined
