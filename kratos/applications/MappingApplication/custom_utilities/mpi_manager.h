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

#if !defined(KRATOS_MPI_MANAGER_H_INCLUDED )
#define  KRATOS_MPI_MANAGER_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include "mpi.h"

// External includes

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
  class MPIManager
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of MPIManager
      KRATOS_CLASS_POINTER_DEFINITION(MPIManager);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      MPIManager();

      /// Destructor.
      virtual ~MPIManager();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{
      // template<class T>
      // static inline MPI_Datatype GetMPIDatatype(const T& value);

      template<class T>
      static inline int SizeOfVariable(const T& variable);

      template<typename T>
      static void MpiSendRecv(T* send_buffer, T* receive_buffer, const int send_buffer_size,
                              int& receive_buffer_size, const int comm_partner, MPI_Status& status) {
          // Determime data type
  		  MPI_Datatype DataType = MPIManager::GetMPIDatatype(T());

          //   Exchange the information about the receiving buffer size
          MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                       1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, &status);

          // Perform the actual exchange of the buffers
          MPI_Sendrecv(send_buffer, send_buffer_size, DataType, comm_partner, 0, receive_buffer,
                       receive_buffer_size, DataType, comm_partner, 0, MPI_COMM_WORLD, &status);
      }

      // This function checks if the buffer sizes are large enough
      template<typename T>
      static void MpiSendRecv(T* send_buffer, T* receive_buffer, const int send_buffer_size,
                              int& receive_buffer_size, const int max_send_buffer_size,
                              const int max_receive_buffer_size, const int comm_partner, MPI_Status& status) {
          // Determime data type
          MPI_Datatype DataType = MPIManager::GetMPIDatatype(T());

          //   Exchange the information about the receiving buffer size
          MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                       1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, &status);

          if (send_buffer_size > max_send_buffer_size)
              KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; \"MpiSendRecv\" Send Buffer too small!  \
              send_buffer_size = " << send_buffer_size << ", max_send_buffer_size = " << max_send_buffer_size << std::endl;

          if (receive_buffer_size > max_receive_buffer_size)
              KRATOS_ERROR << "MappingApplication; MapperMPICommunicator; \"MpiSendRecv\" Receive Buffer too small!  \
              receive_buffer_size = " << receive_buffer_size << ", max_receive_buffer_size = " << max_receive_buffer_size << std::endl;

          // std::cout << "send_buffer_size = " << send_buffer_size << ", max_send_buffer_size = " << max_send_buffer_size << std::endl;
          // std::cout << "receive_buffer_size = " << receive_buffer_size << ", max_receive_buffer_size = " << max_receive_buffer_size << std::endl;

          // Perform the actual exchange of the buffers
          MPI_Sendrecv(send_buffer, send_buffer_size, DataType, comm_partner, 0, receive_buffer,
                       receive_buffer_size, DataType, comm_partner, 0, MPI_COMM_WORLD, &status);
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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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
      MPIManager& operator=(MPIManager const& rOther);

      /// Copy constructor.
      MPIManager(MPIManager const& rOther);


      ///@}

    }; // Class MPIManager

  ///@}

  template<>
  inline MPI_Datatype MPIManager::GetMPIDatatype<int>(const int& value) {
      return MPI_INT ;
  }

  template<>
  inline MPI_Datatype MPIManager::GetMPIDatatype<double>(const double& value) {
      return MPI_DOUBLE ;
  }

  template<>
  inline int MPIManager::SizeOfVariable< double >(const double& variable) {
      return 1;
  }

  template<>
  inline int MPIManager::SizeOfVariable< array_1d<double,3> >(const array_1d<double,3>& variable) {
      return 3;
  }

  // template<>
  // inline int MPIManager::SizeOfVariable< Variable<double> >(const Variable<double>& variable) {
  //     return 1;
  // }
  //
  // template<>
  // inline int MPIManager::SizeOfVariable< Variable< array_1d<double,3> > >(const Variable< array_1d<double,3> >& variable) {
  //     return 3;
  // }
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    MPIManager& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const MPIManager& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPI_MANAGER_H_INCLUDED  defined
