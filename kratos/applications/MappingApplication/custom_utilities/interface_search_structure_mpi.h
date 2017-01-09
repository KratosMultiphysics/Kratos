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

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <iomanip>

// External includes


// Project includes
#include "includes/define.h"
#include "interface_search_structure.h"
#include "processes/graph_coloring_process.h"
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
  class InterfaceSearchStructureMPI : public InterfaceSearchStructure
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterfaceSearchStructureMPI
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceSearchStructureMPI);

      ///@}
      ///@name Life Cycle
      ///@{

      InterfaceSearchStructureMPI(InterfaceObjectManagerBase::Pointer i_interface_object_manager,
                                  InterfaceObjectManagerBase::Pointer i_interface_object_manager_bins,
                                  ModelPart& i_model_part_bins, int i_comm_rank, int i_comm_size,
                                  int i_echo_level) :
                                  InterfaceSearchStructure( i_interface_object_manager,
                                  i_interface_object_manager_bins, i_echo_level) {
          m_comm_rank = i_comm_rank;
          m_comm_size = i_comm_size;

          ComputeLocalBoundingBox(i_model_part_bins);
      }

      /// Destructor.
      virtual ~InterfaceSearchStructureMPI() {
        delete [] m_local_bounding_box;
      }


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{


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
           buffer << "InterfaceSearchStructureMPI" ;
           return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfaceSearchStructureMPI";}

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

      double* m_local_bounding_box = new double[6] {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
                                                  // xmax, xmin,  ymax, ymin,  zmax, zmin

      int m_send_buffer_size;
      int m_receive_buffer_size;

      GraphType m_domains_colored_graph;  // TODO pointer?
      int m_max_colors;  // size aka the number of columns

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void ConductSearchIteration() override {
          CandidateManager candidate_manager;
          FindNeighborCandidates(candidate_manager);
          SelectNeighbors(candidate_manager);
      }

      void GetMPIData(int& max_send_buffer_size, int& max_receive_buffer_size,
                      GraphType& colored_graph, int& max_colors) override {
          // sizes are switched bcs transfer direction is inverted from searching to mapping
          max_send_buffer_size = m_receive_buffer_size;
          max_receive_buffer_size = m_send_buffer_size;

          colored_graph = m_domains_colored_graph;
          max_colors = m_max_colors;
      }

      void PrepareSearching(CandidateManager& candidate_manager, int& max_send_buffer_size,
                            int& max_receive_buffer_size, bool& local_search_required,
                            GraphType& domains_colored_graph, int& max_colors) {
          double* global_bounding_boxes = new double[6 * m_comm_size];
          ComputeGlobalBoundingBoxes(global_bounding_boxes);

          int* local_comm_list = new int[m_comm_size]();
          int* local_memory_size_array = new int[m_comm_size]();

          m_interface_object_manager->ComputeCandidatePartitions(candidate_manager,local_comm_list,
                                                                 local_memory_size_array,
                                                                 global_bounding_boxes);

          if (local_comm_list[m_comm_rank] == 1)
              local_search_required = true;

          ComputeMaxBufferSizes(local_memory_size_array, m_comm_size, m_comm_rank);
          // Buffers are switched for sending back the results, therefore we have to assign the maximum!
          // In one direction 3 coordinates are sent, whereras only one distance is sent back
          max_send_buffer_size = std::max(m_send_buffer_size * 3, m_receive_buffer_size);
          max_receive_buffer_size = std::max(m_receive_buffer_size * 3, m_send_buffer_size);

          ComputeColoringGraph(local_comm_list, m_comm_size, domains_colored_graph, max_colors);
          // Output the colored Graph
          if (m_comm_rank == 0 && m_echo_level > 1) {
              PrintGraph(domains_colored_graph, max_colors);
          }

          delete [] local_comm_list;
          delete [] local_memory_size_array;
          delete [] global_bounding_boxes;
      }

      void FindNeighborCandidates(CandidateManager& candidate_manager) {
          // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
          int max_send_buffer_size;
          int max_receive_buffer_size;
          bool local_search_required = false;
          GraphType domains_colored_graph;
          int max_colors;

          PrepareSearching(candidate_manager, max_send_buffer_size, max_receive_buffer_size,
                           local_search_required, domains_colored_graph, max_colors);

          double* send_buffer = new double[max_send_buffer_size];
          double* receive_buffer = new double[max_receive_buffer_size];

          InterfaceObjectConfigure::ContainerType remote_p_interface_object_list(m_receive_buffer_size);
          std::vector<InterfaceObject::Pointer> candidate_receive_objects(m_receive_buffer_size);
          std::vector<double> min_distances(m_receive_buffer_size);
          std::vector<array_1d<double,2>> local_coordinates(m_receive_buffer_size);
          std::vector<std::vector<double>> shape_functions(m_receive_buffer_size);
          // auto local_coordinates = boost::shared_ptr<std::vector<array_1d<double,2>>>(new std::vector<array_1d<double,2>>(max_receive_buffer_size));

          int send_buffer_size;
          int receive_buffer_size;
          int num_objects;

          // search in local partition
          if (local_search_required) { // search in own partition is required
              m_interface_object_manager->FillBufferLocalSearch(candidate_manager, remote_p_interface_object_list, num_objects);

              FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                                 min_distances, local_coordinates, shape_functions, m_comm_rank);

              m_interface_object_manager_bins->StoreTempSearchResults(candidate_manager, candidate_receive_objects,
                                                                      shape_functions, local_coordinates, m_comm_rank);
              m_interface_object_manager->PostProcessReceivedResults(candidate_manager, min_distances, m_comm_rank);
          }

          // search in remote partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(m_comm_rank, i); // get the partner rank
              if (comm_partner != -1) {
                  m_interface_object_manager->FillSendBufferRemoteSearch(candidate_manager, send_buffer,
                                                                         send_buffer_size, comm_partner);

                  MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                  receive_buffer_size, max_send_buffer_size,
                                                  max_receive_buffer_size, comm_partner);

                  m_interface_object_manager_bins->ProcessReceiveBuffer(remote_p_interface_object_list, receive_buffer,
                                                                        receive_buffer_size, num_objects);

                  // Perform local Search
                  FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                                     min_distances, local_coordinates, shape_functions, m_comm_rank);

                  m_interface_object_manager_bins->StoreTempSearchResults(candidate_manager, candidate_receive_objects,
                                                                          shape_functions, local_coordinates, comm_partner);

                  // Send results back / receive results (distances)
                  int tmp_var = receive_buffer_size / 3;
                  receive_buffer_size = send_buffer_size / 3;
                  send_buffer_size = tmp_var;

                  m_interface_object_manager_bins->FillSendBufferWithResults(send_buffer, send_buffer_size, min_distances);

                  MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, comm_partner, 0, receive_buffer,
                               receive_buffer_size, MPI_DOUBLE, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                  m_interface_object_manager->PostProcessReceivedResults(candidate_manager, receive_buffer, comm_partner);

              } // if I am communicating in this loop (comm_partner != -1)
          } // loop communications

          delete [] send_buffer;
          delete [] receive_buffer;

          MPI_Barrier(MPI_COMM_WORLD);
      }

      void PrepareSelection(CandidateManager& candidate_manager,
                            bool& local_search_required,
                            GraphType& domains_colored_graph,
                            int& max_colors) {
          // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
          int* local_comm_list = new int[m_comm_size]();
          int* local_memory_size_array = new int[m_comm_size]();

          m_interface_object_manager->PrepareMatching(candidate_manager,
                                                      local_comm_list,
                                                      local_memory_size_array);

          if (local_comm_list[m_comm_rank] == 1)
              local_search_required = true;

          // Compute Buffer Sizes for Mapping
          ComputeMaxBufferSizes(local_memory_size_array, m_comm_size, m_comm_rank);

          ComputeColoringGraph(local_comm_list, m_comm_size, domains_colored_graph, max_colors);

          m_domains_colored_graph = domains_colored_graph; // save it for the mapping
          m_max_colors = max_colors;
          // Output the colored Graph
          if (m_comm_rank == 0 && m_echo_level > 1) {
              PrintGraph(domains_colored_graph, max_colors);
          }
          delete [] local_comm_list;
          delete [] local_memory_size_array;
      }

      void SelectNeighbors(CandidateManager& candidate_manager) {
          // // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
          int* send_buffer = new int[m_send_buffer_size];
          int* receive_buffer = new int[m_receive_buffer_size];
          // for debugging, i.e. to check the max buffer sizes
          int max_send_buffer_init_size = m_send_buffer_size;
          int max_receive_buffer_init_size = m_receive_buffer_size;

          bool local_search_required = false;
          GraphType domains_colored_graph;
          int max_colors;
          PrepareSelection(candidate_manager, local_search_required,
                           domains_colored_graph, max_colors);

          int send_buffer_size;
          int receive_buffer_size;

          // search in local partition
          if (local_search_required) { // search in own partition is required
              m_interface_object_manager->FillSendBufferWithMatchInformation(candidate_manager, send_buffer,
                                                                             send_buffer_size, m_comm_rank);
              m_interface_object_manager_bins->ProcessMatchInformation(candidate_manager, send_buffer,
                                                                       send_buffer_size, m_comm_rank);
          }

          // search in remote partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(m_comm_rank, i); // get the partner rank
              if (comm_partner != -1) {
                  m_interface_object_manager->FillSendBufferWithMatchInformation(candidate_manager, send_buffer,
                                                                                 send_buffer_size, comm_partner);

                  MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                  receive_buffer_size, max_send_buffer_init_size,
                                                  max_receive_buffer_init_size, comm_partner);

                  m_interface_object_manager_bins->ProcessMatchInformation(candidate_manager, receive_buffer,
                                                                           receive_buffer_size, comm_partner);
              } // if I am communicating in this loop (comm_partner != -1)
          } // loop communications

          delete [] send_buffer;
          delete [] receive_buffer;

          MPI_Barrier(MPI_COMM_WORLD);
      }

      void ComputeMaxBufferSizes(int* local_memory_size_array, const int comm_size,
                                 const int comm_rank) {
          int* global_memory_size_array = new int[comm_size]();

          // Go through list of how many objects will be sent and find maximum
          m_send_buffer_size = local_memory_size_array[0];
          for (int i = 1; i < comm_size; ++i) {
              m_send_buffer_size = std::max(local_memory_size_array[i], m_send_buffer_size);
          }

          // TODO comment that it could be done more efficient, but only with a lot of code => what Charlie explained
          // Basically splitting the matrix and doing partial reduces
          // Exchange information about how much is going to be sent around among the
          // partitions and take the maximum that my partition will receive
          MPI_Allreduce(local_memory_size_array, global_memory_size_array,
                        comm_size, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          m_receive_buffer_size = global_memory_size_array[comm_rank];

          delete [] global_memory_size_array;
      }

      void ComputeColoringGraph(int* local_comm_list, int num_partitions,
                                GraphType& domains_colored_graph, int& max_colors) {
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

      void PrintGraph(GraphType& graph, int num_color) {
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

      void ComputeLocalBoundingBox(ModelPart& model_part) {
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          for (auto &node : model_part.GetCommunicator().LocalMesh().Nodes()) { // loop over local nodes
              m_local_bounding_box[0] = std::max(node.X(), m_local_bounding_box[0]);
              m_local_bounding_box[1] = std::min(node.X(), m_local_bounding_box[1]);
              m_local_bounding_box[2] = std::max(node.Y(), m_local_bounding_box[2]);
              m_local_bounding_box[3] = std::min(node.Y(), m_local_bounding_box[3]);
              m_local_bounding_box[4] = std::max(node.Z(), m_local_bounding_box[4]);
              m_local_bounding_box[5] = std::min(node.Z(), m_local_bounding_box[5]);
          }
      }

      void ComputeLocalBoundingBoxWithTolerance(double* local_bounding_box_tol) {
          // xmax, xmin,  ymax, ymin,  zmax, zmin
          double tolerance = m_search_radius;

          local_bounding_box_tol[0] = m_local_bounding_box[0] + tolerance;
          local_bounding_box_tol[1] = m_local_bounding_box[1] - tolerance;
          local_bounding_box_tol[2] = m_local_bounding_box[2] + tolerance;
          local_bounding_box_tol[3] = m_local_bounding_box[3] - tolerance;
          local_bounding_box_tol[4] = m_local_bounding_box[4] + tolerance;
          local_bounding_box_tol[5] = m_local_bounding_box[5] - tolerance;

          // TODO Echo
          // PrintLocalBoundingBox(); // uncomment for printing the bounding boxes
      }

      void PrintLocalBoundingBox() {
          double* local_bounding_box_tol = new double[6];
          ComputeLocalBoundingBoxWithTolerance(local_bounding_box_tol);

          std::cout << "\nBounding Box, Rank " << m_comm_rank << " [ "
                    << local_bounding_box_tol[1] << " "     // xmin
                    << local_bounding_box_tol[3] << " "     // ymin
                    << local_bounding_box_tol[5] << " ] [ " // zmin
                    << local_bounding_box_tol[0] << " "     // xmax
                    << local_bounding_box_tol[2] << " "     // ymax
                    << local_bounding_box_tol[4] << " ]"    // zmax
                    << std::endl;
          // TODO maybe write it such that gid can directly read it (=> batch file?)
      }

      void ComputeGlobalBoundingBoxes(double* global_bounding_boxes) {
          double* local_bounding_box_tol = new double[6];
          ComputeLocalBoundingBoxWithTolerance(local_bounding_box_tol);
          MPI_Allgather(local_bounding_box_tol, 6, MPI_DOUBLE,
                        global_bounding_boxes, 6, MPI_DOUBLE, MPI_COMM_WORLD);
          delete local_bounding_box_tol;
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
      InterfaceSearchStructureMPI& operator=(InterfaceSearchStructureMPI const& rOther);

    //   /// Copy constructor.
    //   InterfaceSearchStructureMPI(InterfaceSearchStructureMPI const& rOther){}


      ///@}

    }; // Class InterfaceSearchStructureMPI

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceSearchStructureMPI& rThis)
    {
        return rIStream;
    }

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceSearchStructureMPI& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED  defined
