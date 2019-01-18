//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_MPI_POINT_CONFIGURE_INCLUDED)
#define  KRATOS_MPI_POINT_CONFIGURE_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <limits>
#include <cmath>

#include "mpi.h"
#include "mpi/includes/mpi_communicator.h"
#include "spatial_containers/configures/point_configure.h"

namespace Kratos {

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

/** Configuration file for Points over MPI
 * This class provides extends the PointCofigure class to provide
 * the functions to work over MPI. Specifically this class adds the following functionality:
 * - ReduceIds
 * - TransferObjects
 * - UpdateLocalInterface
 * - UpdateGhostInterface
 */
class MpiPointConfigure : public PointConfigure {
public:

  /// Pointer definition of PointConfigure
  KRATOS_CLASS_POINTER_DEFINITION(PointConfigure);

  /// Base type
  typedef PointConfigure BaseType;

  ///@}
  ///@name Life Cycle
  ///@{

  ///@}
  ///@name Operators
  ///@{

  ///@}
  ///@name Operations
  ///@{

  /** Reassign the id's of all objects.
   * Reassign the id's of all objects so that they are consecutive among all processes.
   * @param total_elements Number of elements from all partitions
   * @param offset         desired offset for the id's which will be applied globaly
   * @param firstObject  Id of the first object of every partition after calculating the new id's
   */
  static inline void ReduceIds(int& total_elements,const int& offset, int& firstObject) {
    int mpi_rank;
    int mpi_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    std::vector<int> reduceArray(mpi_size+1);

    MPI_Allgather(&total_elements,1,MPI_INT,&reduceArray[1],1,MPI_INT,MPI_COMM_WORLD);

    reduceArray[0] = 1 + offset;
    for(int i = 1; i <= mpi_size; i++) {
      reduceArray[i] += reduceArray[i-1] + 1;
    }

    firstObject = reduceArray[mpi_rank];
  }

  /** Transfers the objects from one partition to another
   * Transfers the objects in the SendObjects buffer from the local modelpart to the remote modelpart
   * and recieves objects from the other modelparts into this one.
   * This function does not make use of the internal Communicator transfer function associated to the modelpart.
   * @param SendObjects Vector of containers of objects to be sent
   * @param RecvObjects Vector of containers of objects received
   */
  template<class TObjectType>
  static inline void TransferObjects(Communicator * Communicator, TObjectType& SendObjects, TObjectType& RecvObjects) {
    MpiPointConfigure::TransferObjects(SendObjects, RecvObjects);
  }

  /** Transfers the objects from one partition to another
   * Transfers the objects in the SendObjects buffer from the local modelpart to the remote modelpart
   * and recieves objects from the other modelparts into this one.
   * This function does not make use of the internal Communicator transfer function associated to the modelpart.
   * @param SendObjects Vector of containers of objects to be sent
   * @param RecvObjects Vector of containers of objects received
   */
  template<class TObjectType>
  static inline void TransferObjects(TObjectType& SendObjects, TObjectType& RecvObjects) {

    int mpi_rank;
    int mpi_size;

    Kratos::Serializer particleSerializer;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int * msgSendSize = new int[mpi_size];
    int * msgRecvSize = new int[mpi_size];

    char ** message = new char * [mpi_size];
    char ** mpi_send_buffer = new char * [mpi_size];

    for(int i = 0; i < mpi_size; i++) {
      msgSendSize[i] = 0;
      msgRecvSize[i] = 0;
    }

    for(int i = 0; i < mpi_size; i++) {
      if(mpi_rank != i) {
        Kratos::Serializer particleSerializer;

        particleSerializer.save("ObjectList",SendObjects[i]);

        std::stringstream * stream = (std::stringstream *)particleSerializer.pGetBuffer();
        const std::string & stream_str = stream->str();
        const char * cstr = stream_str.c_str();

        msgSendSize[i] = sizeof(char) * (stream_str.size()+1);
        mpi_send_buffer[i] = (char *)malloc(msgSendSize[i]);
        memcpy(mpi_send_buffer[i],cstr,msgSendSize[i]);
      }
    }

    MPI_Alltoall(msgSendSize,1,MPI_INT,msgRecvSize,1,MPI_INT,MPI_COMM_WORLD);

    int NumberOfCommunicationEvents      = 0;
    int NumberOfCommunicationEventsIndex = 0;

    for(int j = 0; j < mpi_size; j++) {
      if(j != mpi_rank && msgRecvSize[j]) NumberOfCommunicationEvents++;
      if(j != mpi_rank && msgSendSize[j]) NumberOfCommunicationEvents++;
    }

    MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
    MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

    //Set up all receive and send events
    for(int i = 0; i < mpi_size; i++) {
      if(i != mpi_rank && msgRecvSize[i]) {
        message[i] = (char *)malloc(sizeof(char) * msgRecvSize[i]);

        MPI_Irecv(message[i],msgRecvSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
      }

      if(i != mpi_rank && msgSendSize[i]) {
        MPI_Isend(mpi_send_buffer[i],msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
      }
    }

    //wait untill all communications finish
    int err = MPI_Waitall(NumberOfCommunicationEvents, reqs, stats);

    if(err != MPI_SUCCESS) {
      KRATOS_THROW_ERROR(std::runtime_error,"Error in mpi_communicator","")
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(int i = 0; i < mpi_size; i++) {
      if (i != mpi_rank && msgRecvSize[i]) {
        Kratos::Serializer particleSerializer;
        std::stringstream * serializer_buffer;

        serializer_buffer = (std::stringstream *)particleSerializer.pGetBuffer();
        serializer_buffer->write(message[i], msgRecvSize[i]);

        particleSerializer.load("ObjectList",RecvObjects[i]);
      }

      MPI_Barrier(MPI_COMM_WORLD);
    }

    // Free buffers
    for(int i = 0; i < mpi_size; i++) {
      if(msgRecvSize[i]) {
        free(message[i]);
      }

      if(msgSendSize[i]) {
        free(mpi_send_buffer[i]);
      }
    }

    delete [] reqs;
    delete [] stats;

    delete [] message;
    delete [] mpi_send_buffer;

    delete [] msgSendSize;
    delete [] msgRecvSize;
  }

  /** Updates the local interface.
   * Updates the local interface of the communicator.
   * This function is intentionally empty as sparse points are not part of modelParts
   * @param Communicator  The modelPart Communicator
   * @param Object        Object to be added to the interface
   * @param Destination   Destination of the object in the Remote (Becomes part of the Remote Modelpart Local interface)
   */
  template<class TObjectType>
  static inline void UpdateLocalInterface(Communicator * Communicator, TObjectType& Object, const int & Destination) {
    // THIS FUNCTION IS INTENTONALY EMPTY
  }

  /** Updates the ghost interface.
   * Updates the ghost interface of the communicator.
   * This function is intentionally empty as sparse points are not part of modelParts
   * @param Communicator  The modelPart Communicator
   * @param Object        Object to be added to the interface
   * @param Source        source of the object in the Local (Becomes part of the Local Modelpart Ghost interface)
   */
  template<class TObjectType>
  static inline void UpdateGhostInterface(Communicator * Communicator, TObjectType& Object, const int & Source) {
    // THIS FUNCTION IS INTENTONALY EMPTY
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

  /// Turns back information as a string.
  virtual std::string Info() const {
    return "Spatial Containers Configure for 'Points' over MPI";
  }

  /// Turns back data as a string.
  virtual std::string Data() const {
    return "Dimension: " + std::to_string(Dimension);
  }

  /// Prints object's information.
  virtual void PrintInfo(std::ostream& rOStream) const {
    rOStream << Info() << std::endl;
  }

  /// Prints object's data.
  virtual void PrintData(std::ostream& rOStream) const {
    rOStream << Data() << Dimension << std::endl;
  }

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
  MpiPointConfigure& operator=(MpiPointConfigure const& rOther);

  /// Copy constructor.
  MpiPointConfigure(MpiPointConfigure const& rOther);

  ///@}

}; // Class MpiPointConfigure

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MpiPointConfigure& rThis){
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const MpiPointConfigure& rThis){
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

} // namespace Kratos.
#endif /* KRATOS_MPI_POINT_CONFIGURE_INCLUDED */
