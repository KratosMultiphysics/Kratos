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
#include <iomanip>


// External includes
#include "mpi.h"

// Project includes
#include "includes/define.h"
#include "processes/graph_coloring_process.h"


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

      /// Destructor.
      virtual ~MapperUtilitiesMPI() { }


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
          // Determine data type
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
          // Determine data type
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

      static void ComputeMaxBufferSizes(int* local_memory_size_array,
                                        int& send_buffer_size,
                                        int& receive_buffer_size,
                                        const int comm_rank,
                                        const int comm_size) {
          int* global_memory_size_array = new int[comm_size]();

          // Go through list of how many objects will be sent and find maximum
          send_buffer_size = local_memory_size_array[0];
          for (int i = 1; i < comm_size; ++i) {
              send_buffer_size = std::max(local_memory_size_array[i], send_buffer_size);
          }

          // TODO comment that it could be done more efficient, but only with a lot of code => what Charlie explained
          // Basically splitting the matrix and doing partial reduces
          // Exchange information about how much is going to be sent around among the
          // partitions and take the maximum that my partition will receive
          MPI_Allreduce(local_memory_size_array, global_memory_size_array,
                        comm_size, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          receive_buffer_size = global_memory_size_array[comm_rank];

          delete [] global_memory_size_array;
      }

      static void ComputeColoringGraph(int* local_comm_list,
                                       int num_partitions,
                                       GraphType& domains_colored_graph,
                                       int& max_colors) {
          int comm_list_size = num_partitions * num_partitions;
          int* global_comm_list = new int[comm_list_size]();

          MPI_Allgather(local_comm_list, num_partitions, MPI_INT, global_comm_list,
                        num_partitions, MPI_INT, MPI_COMM_WORLD);

          GraphType domains_graph;
          domains_graph.resize(num_partitions, num_partitions,false);
          domains_graph = ScalarMatrix(num_partitions, num_partitions, -1.0f);

          // set up communication matrix as input for the coloring process
          for (int i = 0; i < num_partitions; ++i) {
              for (int j = 0; j < num_partitions; ++j) {
                  domains_graph(i,j) = global_comm_list[(i * num_partitions) + j];
              }
          }

        // make matrix symmetric, using the larger value (i.e. 1, meaning that communication is required)
          // needed if for some ranks only one-sided communication is required
          for(std::size_t i = 0; i < domains_graph.size1(); ++i) {
              for(std::size_t j = i + 1; j < domains_graph.size2(); ++j) {
                  int val = std::max(domains_graph(i,j), domains_graph(j,i));
                  domains_graph(i,j) = val;
                  domains_graph(j,i) = val;
              }
              domains_graph(i,i) = 0; // set main diagonal to zero, i.e. avoid communication with myself
          }

          GraphColoringProcess(num_partitions, domains_graph,
                               domains_colored_graph, max_colors).Execute();

          delete [] global_comm_list;
      }

      static void PrintGraph(GraphType& graph, int num_color) {
          std::cout << "\nCommunication Graph: " << std::endl;
          std::cout << std::setw(5);
          for(std::size_t i = 0; i < graph.size1(); ++i) {
              for(int j = 0; j < num_color; ++j) {
                  std::cout << graph(i,j) << std::setw(5);
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;
      }

      static void ComputeLocalBoundingBox(ModelPart& model_part,
                                          double* local_bounding_box) {
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          for (auto &node : model_part.GetCommunicator().LocalMesh().Nodes()) { // loop over local nodes
              local_bounding_box[0] = std::max(node.X(), local_bounding_box[0]);
              local_bounding_box[1] = std::min(node.X(), local_bounding_box[1]);
              local_bounding_box[2] = std::max(node.Y(), local_bounding_box[2]);
              local_bounding_box[3] = std::min(node.Y(), local_bounding_box[3]);
              local_bounding_box[4] = std::max(node.Z(), local_bounding_box[4]);
              local_bounding_box[5] = std::min(node.Z(), local_bounding_box[5]);
          }
      }

      static void ComputeLocalBoundingBoxWithTolerance(double* local_bounding_box,
                                                       const double tolerance,
                                                       const int comm_rank,
                                                       const int echo_level,
                                                       double* local_bounding_box_tol) {
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          local_bounding_box_tol[0] = local_bounding_box[0] + tolerance;
          local_bounding_box_tol[1] = local_bounding_box[1] - tolerance;
          local_bounding_box_tol[2] = local_bounding_box[2] + tolerance;
          local_bounding_box_tol[3] = local_bounding_box[3] - tolerance;
          local_bounding_box_tol[4] = local_bounding_box[4] + tolerance;
          local_bounding_box_tol[5] = local_bounding_box[5] - tolerance;

          if (echo_level > 2) {
              MapperUtilitiesMPI::PrintBoundingBox(local_bounding_box_tol, comm_rank);
          }
      }

      static void PrintBoundingBox(double* bounding_box, const int comm_rank) {
          std::cout << "\nBounding Box, Rank " << comm_rank << " [ "
                    << bounding_box[1] << " "     // xmin
                    << bounding_box[3] << " "     // ymin
                    << bounding_box[5] << " ] [ " // zmin
                    << bounding_box[0] << " "     // xmax
                    << bounding_box[2] << " "     // ymax
                    << bounding_box[4] << " ]"    // zmax
                    << std::endl;
          // TODO maybe write it such that gid can directly read it (=> batch file?)
      }

      static void ComputeGlobalBoundingBoxes(double* local_bounding_box,
                                             const double tolerance,
                                             const int comm_rank,
                                             const int echo_level,
                                             double* global_bounding_boxes) {
          double* local_bounding_box_tol = new double[6];
          MapperUtilitiesMPI::ComputeLocalBoundingBoxWithTolerance(local_bounding_box,
                                                                   tolerance,
                                                                   comm_rank,
                                                                   echo_level,
                                                                   local_bounding_box_tol);
          MPI_Allgather(local_bounding_box_tol, 6, MPI_DOUBLE,
                        global_bounding_boxes, 6, MPI_DOUBLE, MPI_COMM_WORLD);
          delete local_bounding_box_tol;
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

      /// Default constructor.
      MapperUtilitiesMPI() { }

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
