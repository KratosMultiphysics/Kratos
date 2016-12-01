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


// External includes


// Project includes
#include "includes/define.h"
// #include "mpi.h" // TODO needed?
#include "interface_search_structure.h"
#include "processes/graph_coloring_process.h"
#include <iomanip>


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

      InterfaceSearchStructureMPI(InterfaceObjectManager::Pointer i_interface_object_manager,
                               InterfaceObjectManager::Pointer i_interface_object_manager_bins,
                               ModelPart& i_model_part_bins, double i_initial_search_radius,
                               int i_max_search_iterations, int i_comm_rank, int i_comm_size) :
                               InterfaceSearchStructure( i_interface_object_manager,
                               i_interface_object_manager_bins, i_initial_search_radius,
                               i_max_search_iterations) {
          // Custom Constructor, here the tolerances are passed from outside
          m_comm_rank = i_comm_rank;
          m_comm_size = i_comm_size;

          ComputeLocalBoundingBox(i_model_part_bins);
      }

      /// Destructor.
      virtual ~InterfaceSearchStructureMPI(){}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Search() override {
          // TODO What happens in the case of remeshing? Reinitialize / redo the whole mapper?
          FindNeighborCandidates();
          SelectNeighbors();
      }

      void Clear() override {
          // TODO Implement


      }

      void GetMPIData(int& max_send_buffer_size, int& max_receive_buffer_size,
                      GraphType& colored_graph, int& max_colors) {
          // sizes are switched bcs transfer direction is inverted from searching to mapping
          max_send_buffer_size = m_max_receive_buffer_size;
          max_receive_buffer_size = m_max_send_buffer_size;

          colored_graph = m_domains_colored_graph;
          max_colors = m_max_colors;
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

      int m_max_send_buffer_size;
      int m_max_receive_buffer_size;

      GraphType m_domains_colored_graph;  // TODO pointer?
      int m_max_colors;  // size aka the number of columns

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void FindNeighborCandidates() {
          // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins

          InterfaceObjectConfigure::ContainerType point_list = m_interface_object_manager->GetPointList();

          int point_list_size = point_list.size();

          MPI_Status status;

          double* global_bounding_boxes = new double[6 * m_comm_size];
        //   MPI_Allgather(m_local_bounding_box, 6, MPI_DOUBLE, global_bounding_boxes, 6, MPI_DOUBLE, MPI_COMM_WORLD);
          ComputeGlobalBoundingBoxes(global_bounding_boxes);

          std::vector<std::vector<double>> candidate_send_points_distances(
                m_comm_size, std::vector<double>(point_list_size, -1.0));
          int* local_comm_list = new int[m_comm_size]();
          int* local_memory_size_array = new int[m_comm_size]();

          m_interface_object_manager->ComputeCandidatePartitions(local_comm_list, local_memory_size_array, global_bounding_boxes);

          int max_send_buffer_size = 0;
          int max_receive_buffer_size = 0;

          ComputeMaxBufferSizes(local_memory_size_array, max_send_buffer_size, max_receive_buffer_size, m_comm_size, m_comm_rank);
          m_max_send_buffer_size = max_send_buffer_size;
          m_max_receive_buffer_size = max_receive_buffer_size;
          // Buffers are switched for sending back the results, therefore we have to assign the maximum!
          // In one direction 3 coordinates are sent, whereras only one distance is sent back
          double* send_buffer = new double[std::max(max_send_buffer_size * 3, max_receive_buffer_size)];
          double* receive_buffer = new double[std::max(max_receive_buffer_size * 3, max_send_buffer_size)];

          InterfaceObjectConfigure::ContainerType remote_p_point_list(max_receive_buffer_size);
          std::vector<InterfaceObject::Pointer> candidate_receive_points(max_receive_buffer_size);
          std::vector<double> min_distances(max_receive_buffer_size);

          int max_colors;
          GraphType domains_colored_graph;
          ComputeColoringGraph(local_comm_list, m_comm_size, domains_colored_graph, max_colors);

          // Output the colored Graph
          if (m_comm_rank == 0){
              PrintGraph(domains_colored_graph, max_colors);
          }
          MPI_Barrier(MPI_COMM_WORLD);

          int send_buffer_size;
          int receive_buffer_size;
          int num_points;

          // search in same partition
          if (local_comm_list[m_comm_rank] == 1) { // search in own partition is required
            //   std::cout << "Rank " << m_comm_rank << " performs search in own partition" << std::endl;

              m_interface_object_manager->FillBufferLocalSearch(remote_p_point_list, num_points);

              FindLocalNeighbors(remote_p_point_list, num_points, candidate_receive_points, min_distances, m_comm_rank);

              m_interface_object_manager_bins->StoreTempSearchResults(candidate_receive_points, m_comm_rank);
              m_interface_object_manager->PostProcessReceivedResults(min_distances, m_comm_rank);
          }

          // search in other partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(m_comm_rank, i); // get the partner rank
              if (comm_partner != -1) {

                  m_interface_object_manager->FillSendBufferRemoteSearch(send_buffer, send_buffer_size, comm_partner);

                //   Exchange the information about the receiving buffer size
                  MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                               1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, &status);

                  // Perform the actual exchange of the buffers
                  MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, comm_partner, 0, receive_buffer,
                               receive_buffer_size, MPI_DOUBLE, comm_partner, 0, MPI_COMM_WORLD, &status);

                  m_interface_object_manager_bins->ProcessReceiveBuffer(remote_p_point_list, receive_buffer, receive_buffer_size, num_points);

                  // Perform local Search
                  FindLocalNeighbors(remote_p_point_list, num_points, candidate_receive_points, min_distances, m_comm_rank);

                  m_interface_object_manager_bins->StoreTempSearchResults(candidate_receive_points, comm_partner);

                  // Send results back / receive results (distances)
                  int tmp_var = receive_buffer_size / 3;
                  receive_buffer_size = send_buffer_size / 3;
                  send_buffer_size = tmp_var;

                  m_interface_object_manager_bins->FillSendBufferWithResults(send_buffer, send_buffer_size, min_distances);

                  MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, comm_partner, 0, receive_buffer,
                               receive_buffer_size, MPI_DOUBLE, comm_partner, 0, MPI_COMM_WORLD, &status);

                  m_interface_object_manager->PostProcessReceivedResults(receive_buffer, comm_partner);

              } // if I am communicating in this loop (comm_partner != -1)
          } // loop communications

          m_interface_object_manager->CheckResults();

          delete [] send_buffer;
          delete [] receive_buffer;
          delete [] local_comm_list;
          delete [] local_memory_size_array;
          delete [] global_bounding_boxes;

          MPI_Barrier(MPI_COMM_WORLD);
      }

      void SelectNeighbors() {
          // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins

          MPI_Status status;

          int* local_comm_list = new int[m_comm_size]();
          int* local_memory_size_array = new int[m_comm_size]();

          m_interface_object_manager->PrepareMatching(local_comm_list, local_memory_size_array);

          int max_send_buffer_size = 0;
          int max_receive_buffer_size = 0;

          ComputeMaxBufferSizes(local_memory_size_array, max_send_buffer_size, max_receive_buffer_size, m_comm_size, m_comm_rank);

          int* send_buffer = new int[m_max_send_buffer_size];
          int* receive_buffer = new int[m_max_receive_buffer_size];

          m_max_send_buffer_size = max_send_buffer_size;
          m_max_receive_buffer_size = max_receive_buffer_size;

          int max_colors;
          GraphType domains_colored_graph;
          ComputeColoringGraph(local_comm_list, m_comm_size, domains_colored_graph, max_colors);

          m_domains_colored_graph = domains_colored_graph; // save it for the mapping
          m_max_colors = max_colors;
          // Output the colored Graph
          if (m_comm_rank == 0){
              PrintGraph(domains_colored_graph, max_colors);
          }
          MPI_Barrier(MPI_COMM_WORLD);

          int send_buffer_size;
          int receive_buffer_size;

          // search in same partition
          if (local_comm_list[m_comm_rank] == 1) { // search in own partition is required
            //   std::cout << "Rank " << m_comm_rank << " performs search in own partition" << std::endl;

              m_interface_object_manager->FillSendBufferWithMatchInformation(send_buffer, send_buffer_size, m_comm_rank);
              m_interface_object_manager_bins->ProcessMatchInformation(send_buffer, send_buffer_size, m_comm_rank);
          }

          // search in other partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(m_comm_rank, i); // get the partner rank
              if (comm_partner != -1) {
                  m_interface_object_manager->FillSendBufferWithMatchInformation(send_buffer, send_buffer_size, comm_partner);

                  //   Exchange the information about the receiving buffer size
                  MPI_Sendrecv(&send_buffer_size, 1, MPI_INT, comm_partner, 0, &receive_buffer_size,
                               1, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, &status);

                  // Perform the actual exchange of the buffers
                  MPI_Sendrecv(send_buffer, send_buffer_size, MPI_INT, comm_partner, 0, receive_buffer,
                               receive_buffer_size, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, &status);

                  m_interface_object_manager_bins->ProcessMatchInformation(receive_buffer, receive_buffer_size, comm_partner);
              } // if I am communicating in this loop (comm_partner != -1)
          } // loop communications

          delete [] send_buffer;
          delete [] receive_buffer;
          delete [] local_comm_list;
          delete [] local_memory_size_array;

          MPI_Barrier(MPI_COMM_WORLD);
      }

      void ComputeMaxBufferSizes(int* local_memory_size_array, int& max_send_buffer_size, int& max_receive_buffer_size, const int comm_size, const int comm_rank) {
          int* global_memory_size_array = new int[comm_size]();

          // Go through list of how many objects will be sent and find maximum
          max_send_buffer_size = local_memory_size_array[0];
          for (int i = 1; i < comm_size; ++i) {
              max_send_buffer_size = std::max(local_memory_size_array[i], max_send_buffer_size);
          }

          // TODO comment that it could be done more efficient, but only with a lot of code => what Charlie explained
          // Basically splitting the matrix and doing partial reduces
          // Exchange information about how much is going to be sent around among the
          // partitions and take the maximum that my partition will receive
          MPI_Allreduce(local_memory_size_array, global_memory_size_array, comm_size, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

          max_receive_buffer_size = global_memory_size_array[comm_rank];

          delete [] global_memory_size_array;
      }

      void ComputeColoringGraph(int* local_comm_list, int num_partitions, GraphType& domains_colored_graph, int& max_colors) {
          int comm_list_size = num_partitions * num_partitions;
          int* global_comm_list = new int[comm_list_size]();

          MPI_Allgather(local_comm_list, num_partitions, MPI_INT, global_comm_list, num_partitions, MPI_INT, MPI_COMM_WORLD);

          GraphType domains_graph;
          domains_graph.resize(num_partitions, num_partitions,false);
          domains_graph = ScalarMatrix(num_partitions, num_partitions, -1.00);

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

          GraphColoringProcess(num_partitions, domains_graph, domains_colored_graph, max_colors).Execute();

          delete [] global_comm_list;
      }

      void PrintGraph(GraphType& graph, int num_color) {
          std::cout << std::endl;
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
          for (auto &point : model_part.GetCommunicator().LocalMesh().Nodes()) { // loop over local nodes
              m_local_bounding_box[0] = std::max(point.X(), m_local_bounding_box[0]);
              m_local_bounding_box[1] = std::min(point.X(), m_local_bounding_box[1]);
              m_local_bounding_box[2] = std::max(point.Y(), m_local_bounding_box[2]);
              m_local_bounding_box[3] = std::min(point.Y(), m_local_bounding_box[3]);
              m_local_bounding_box[4] = std::max(point.Z(), m_local_bounding_box[4]);
              m_local_bounding_box[5] = std::min(point.Z(), m_local_bounding_box[5]);
          }

          double max_element_size = m_initial_search_radius; // already has a tolerance, see base class

          m_local_bounding_box[0] += max_element_size;
          m_local_bounding_box[1] -= max_element_size;
          m_local_bounding_box[2] += max_element_size;
          m_local_bounding_box[3] -= max_element_size;
          m_local_bounding_box[4] += max_element_size;
          m_local_bounding_box[5] -= max_element_size;

          // PrintLocalBoundingBox(); // uncomment for printing the bounding boxes
      }

      void PrintLocalBoundingBox() {
          std::cout << "Bounding Box rank " << m_comm_rank << " [ " << m_local_bounding_box[1] <<
          " " << m_local_bounding_box[3] << " " << m_local_bounding_box[5] << " ] [ " << m_local_bounding_box[0] <<
          " " << m_local_bounding_box[2] << " " << m_local_bounding_box[4] << " ]" << std::endl;
          // TODO maybe write it such that gid can directly read it (=> batch file?)
          // see "applications/wind_turbine_application/custom_utilities/wind_turbine_rotation_utilities.cpp"
          // and search for "gid" to find an example of how to print to a modelpart
      }

      void ComputeGlobalBoundingBoxes(double* global_bounding_boxes) {
          MPI_Allgather(m_local_bounding_box, 6, MPI_DOUBLE, global_bounding_boxes, 6, MPI_DOUBLE, MPI_COMM_WORLD);
      }

      void ComputeGlobalSearchTolerance(ModelPart& model_part) {
          model_part.GetCommunicator().MaxAll(m_initial_search_radius);
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
