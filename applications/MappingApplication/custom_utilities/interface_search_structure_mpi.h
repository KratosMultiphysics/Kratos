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
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_MPI_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "interface_search_structure.h"
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

      InterfaceSearchStructureMPI(InterfaceObjectManagerBase::Pointer pInterfaceObjectManager,
                                  InterfaceObjectManagerBase::Pointer pInterfaceObjectManagerBins,
                                  ModelPart& rModelPartBins, int CommRank, int CommSize,
                                  int EchoLevel) :
                                  InterfaceSearchStructure( pInterfaceObjectManager,
                                  pInterfaceObjectManagerBins, EchoLevel) {
          mCommRank = CommRank;
          mCommSize = CommSize;

          MapperUtilitiesMPI::ComputeLocalBoundingBox(rModelPartBins,
                                                      mLocalBoundingBox);
      }

      /// Destructor.
      virtual ~InterfaceSearchStructureMPI() {
          delete [] mLocalBoundingBox;
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

      double* mLocalBoundingBox = new double[6] {-1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
                                                  // xmax, xmin,  ymax, ymin,  zmax, zmin

      int m_send_buffer_size;
      int m_receive_buffer_size;

      GraphType m_domains_colored_graph;
      int m_max_colors;  // size aka the number of columns

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void ConductSearchIteration(const bool LastIteration) override {
          CandidateManager candidate_manager;
          FindNeighborCandidates(candidate_manager, LastIteration);
          SelectNeighbors(candidate_manager);
      }

      void PrepareSearching(CandidateManager& candidate_manager, int& max_send_buffer_size,
                            int& max_receive_buffer_size, bool& local_search_required,
                            GraphType& domains_colored_graph, int& max_colors,
                            const bool LastIteration) {
          double* global_bounding_boxes = new double[6 * mCommSize];
          MapperUtilitiesMPI::ComputeGlobalBoundingBoxes(mLocalBoundingBox,
                                                         mSearchRadius,
                                                         mCommRank,
                                                         mEchoLevel,
                                                         global_bounding_boxes);

          int* local_comm_list = new int[mCommSize]();
          int* local_memory_size_array = new int[mCommSize]();

          mpInterfaceObjectManager->ComputeCandidatePartitions(candidate_manager,local_comm_list,
                                                                   local_memory_size_array,
                                                                   global_bounding_boxes,
                                                                   LastIteration);

          if (local_comm_list[mCommRank] == 1)
              local_search_required = true;

          MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                                                    m_send_buffer_size,
                                                    m_receive_buffer_size,
                                                    mCommRank, mCommSize);
          // Buffers are switched for sending back the results, therefore we have to assign the maximum!
          // In one direction 3 coordinates are sent, whereras only one distance is sent back
          max_send_buffer_size = std::max(m_send_buffer_size * 3, m_receive_buffer_size);
          max_receive_buffer_size = std::max(m_receive_buffer_size * 3, m_send_buffer_size);

          MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, mCommSize,
                                                   domains_colored_graph, max_colors);
          // Output the colored Graph
          if (mCommRank == 0 && mEchoLevel > 2) {
              MapperUtilitiesMPI::PrintGraph(domains_colored_graph, max_colors);
          }

          delete [] local_comm_list;
          delete [] local_memory_size_array;
          delete [] global_bounding_boxes;
      }

      void FindNeighborCandidates(CandidateManager& candidate_manager,
                                  const bool LastIteration) {
          // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
          int max_send_buffer_size;
          int max_receive_buffer_size;
          bool local_search_required = false;
          GraphType domains_colored_graph;
          int max_colors;

          PrepareSearching(candidate_manager, max_send_buffer_size,
                           max_receive_buffer_size, local_search_required,
                           domains_colored_graph, max_colors, LastIteration);

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
              mpInterfaceObjectManager->FillBufferLocalSearch(candidate_manager, remote_p_interface_object_list, num_objects);

              FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                                 min_distances, local_coordinates, shape_functions);

              mpInterfaceObjectManagerBins->StoreTempSearchResults(candidate_manager, candidate_receive_objects,
                                                                      shape_functions, local_coordinates, mCommRank);
              mpInterfaceObjectManager->PostProcessReceivedResults(candidate_manager, min_distances, mCommRank);
          }

          // search in remote partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(mCommRank, i); // get the partner rank
              if (comm_partner != -1) {
                  mpInterfaceObjectManager->FillSendBufferRemoteSearch(candidate_manager, send_buffer,
                                                                         send_buffer_size, comm_partner);

                  MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                  receive_buffer_size, max_send_buffer_size,
                                                  max_receive_buffer_size, comm_partner);

                  mpInterfaceObjectManagerBins->ProcessReceiveBuffer(remote_p_interface_object_list, receive_buffer,
                                                                        receive_buffer_size, num_objects);

                  // Perform local Search
                  FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                                     min_distances, local_coordinates, shape_functions);

                  mpInterfaceObjectManagerBins->StoreTempSearchResults(candidate_manager, candidate_receive_objects,
                                                                          shape_functions, local_coordinates, comm_partner);

                  // Send results back / receive results (distances)
                  int tmp_var = receive_buffer_size / 3;
                  receive_buffer_size = send_buffer_size / 3;
                  send_buffer_size = tmp_var;

                  mpInterfaceObjectManagerBins->FillSendBufferWithResults(send_buffer, send_buffer_size, min_distances);

                  MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, comm_partner, 0, receive_buffer,
                               receive_buffer_size, MPI_DOUBLE, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                  mpInterfaceObjectManager->PostProcessReceivedResults(candidate_manager, receive_buffer, comm_partner);

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
          int* local_comm_list = new int[mCommSize]();
          int* local_memory_size_array = new int[mCommSize]();

          mpInterfaceObjectManager->PrepareMatching(candidate_manager,
                                                      local_comm_list,
                                                      local_memory_size_array);

          if (local_comm_list[mCommRank] == 1) {
              local_search_required = true;
          }

          // Compute Buffer Sizes for Mapping
          MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                                                    m_send_buffer_size,
                                                    m_receive_buffer_size,
                                                    mCommRank, mCommSize);

          MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, mCommSize, domains_colored_graph, max_colors);

          m_domains_colored_graph = domains_colored_graph; // save it for the mapping
          m_max_colors = max_colors;
          // Output the colored Graph
          if (mCommRank == 0 && mEchoLevel > 2) {
              MapperUtilitiesMPI::PrintGraph(domains_colored_graph, max_colors);
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
              mpInterfaceObjectManager->FillSendBufferWithMatchInformation(candidate_manager, send_buffer,
                                                                             send_buffer_size, mCommRank);
              mpInterfaceObjectManagerBins->ProcessMatchInformation(candidate_manager, send_buffer,
                                                                       send_buffer_size, mCommRank);
          }

          // search in remote partitions
          for (int i = 0; i < max_colors; ++i) { // loop over communication steps
              int comm_partner = domains_colored_graph(mCommRank, i); // get the partner rank
              if (comm_partner != -1) {
                  mpInterfaceObjectManager->FillSendBufferWithMatchInformation(candidate_manager, send_buffer,
                                                                                 send_buffer_size, comm_partner);

                  MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                  receive_buffer_size, max_send_buffer_init_size,
                                                  max_receive_buffer_init_size, comm_partner);

                  mpInterfaceObjectManagerBins->ProcessMatchInformation(candidate_manager, receive_buffer,
                                                                           receive_buffer_size, comm_partner);
              } // if I am communicating in this loop (comm_partner != -1)
          } // loop communications

          delete [] send_buffer;
          delete [] receive_buffer;

          MPI_Barrier(MPI_COMM_WORLD);
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
