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

// External includes

// Project includes
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
                                  pInterfaceObjectManagerBins, EchoLevel)
    {
        mCommRank = CommRank;
        mCommSize = CommSize;

        MapperUtilitiesMPI::ComputeLocalBoundingBox(rModelPartBins,
                mLocalBoundingBox);
    }

    /// Destructor.
    virtual ~InterfaceSearchStructureMPI()
    {
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "InterfaceSearchStructureMPI" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "InterfaceSearchStructureMPI";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    double* mLocalBoundingBox = new double[6] { -1e10, 1e10, -1e10, 1e10, -1e10, 1e10}; // initialize "inverted"
    // xmax, xmin,  ymax, ymin,  zmax, zmin

    int mSendBufferSize;
    int mReceiveBufferSize;

    GraphType mDomainsColoredGraph;
    int mMaxColors;  // size aka the number of columns

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void ConductSearchIteration(const bool LastIteration) override
    {
        CandidateManager candidate_manager;
        FindNeighborCandidates(candidate_manager, LastIteration);
        SelectNeighbors(candidate_manager);
    }

    void PrepareSearching(CandidateManager& rCandidateManager, int& rMaxSendBufferSize,
                          int& rMaxReceiveBufferSize, bool& rLocalSearchRequired,
                          GraphType& rDomainsColoredGraph, int& rMaxColors,
                          const bool LastIteration)
    {
        double* global_bounding_boxes = new double[6 * mCommSize];
        MapperUtilitiesMPI::ComputeGlobalBoundingBoxes(mLocalBoundingBox,
                mSearchRadius,
                mCommRank,
                mEchoLevel,
                global_bounding_boxes);

        int* local_comm_list = new int[mCommSize]();
        int* local_memory_size_array = new int[mCommSize]();

        mpInterfaceObjectManager->ComputeCandidatePartitions(rCandidateManager, local_comm_list,
                local_memory_size_array,
                global_bounding_boxes,
                LastIteration);

        if (local_comm_list[mCommRank] == 1)
            rLocalSearchRequired = true;

        MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                mSendBufferSize,
                mReceiveBufferSize,
                mCommRank, mCommSize);
        // Buffers are switched for sending back the results, therefore we have to assign the maximum!
        // In one direction 3 coordinates are sent, whereras only one distance is sent back
        rMaxSendBufferSize = std::max(mSendBufferSize * 3, mReceiveBufferSize);
        rMaxReceiveBufferSize = std::max(mReceiveBufferSize * 3, mSendBufferSize);

        MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, mCommSize,
                rDomainsColoredGraph, rMaxColors);
        // Output the colored Graph
        if (mCommRank == 0 && mEchoLevel > 2)
        {
            MapperUtilitiesMPI::PrintGraph(rDomainsColoredGraph, rMaxColors);
        }

        delete [] local_comm_list;
        delete [] local_memory_size_array;
        delete [] global_bounding_boxes;
    }

    void FindNeighborCandidates(CandidateManager& rCandidateManager,
                                const bool LastIteration)
    {
        // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
        int max_send_buffer_size;
        int max_receive_buffer_size;
        bool local_search_required = false;
        GraphType domains_colored_graph;
        int max_colors;

        PrepareSearching(rCandidateManager, max_send_buffer_size,
                         max_receive_buffer_size, local_search_required,
                         domains_colored_graph, max_colors, LastIteration);

        double* send_buffer = new double[max_send_buffer_size];
        int* send_buffer_int = new int[max_send_buffer_size];
        double* receive_buffer = new double[max_receive_buffer_size];
        int* receive_buffer_int = new int[max_receive_buffer_size];

        InterfaceObjectConfigure::ContainerType remote_p_interface_object_list(mReceiveBufferSize);
        std::vector<InterfaceObject::Pointer> candidate_receive_objects(mReceiveBufferSize);
        std::vector<double> min_distances(mReceiveBufferSize);
        std::vector<std::vector<double>> shape_function_values(mReceiveBufferSize);
        std::vector<int> pairing_indices(mReceiveBufferSize);
        // auto local_coordinates = boost::shared_ptr<std::vector<array_1d<double,2>>>(new std::vector<array_1d<double,2>>(max_receive_buffer_size));

        int send_buffer_size;
        int receive_buffer_size;
        int num_objects;

        // search in local partition
        if (local_search_required)   // search in own partition is required
        {
            mpInterfaceObjectManager->FillBufferLocalSearch(rCandidateManager, remote_p_interface_object_list, num_objects);

            FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                               min_distances, shape_function_values, pairing_indices);

            mpInterfaceObjectManagerBins->StoreTempSearchResults(rCandidateManager, candidate_receive_objects,
                    shape_function_values, mCommRank);
            mpInterfaceObjectManager->PostProcessReceivedResults(rCandidateManager, min_distances,
                    pairing_indices, mCommRank);
        }

        // search in remote partitions
        for (int i = 0; i < max_colors; ++i)   // loop over communication steps
        {
            int comm_partner = domains_colored_graph(mCommRank, i); // get the partner rank
            if (comm_partner != -1)
            {
                mpInterfaceObjectManager->FillSendBufferRemoteSearch(rCandidateManager, send_buffer,
                        send_buffer_size, comm_partner);

                MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                receive_buffer_size, max_send_buffer_size,
                                                max_receive_buffer_size, comm_partner);

                mpInterfaceObjectManagerBins->ProcessReceiveBuffer(remote_p_interface_object_list, receive_buffer,
                        receive_buffer_size, num_objects);

                // Perform local Search
                FindLocalNeighbors(remote_p_interface_object_list, num_objects, candidate_receive_objects,
                                   min_distances, shape_function_values, pairing_indices);

                mpInterfaceObjectManagerBins->StoreTempSearchResults(rCandidateManager, candidate_receive_objects,
                        shape_function_values, comm_partner);

                // Send results back / receive results (distances)
                int tmp_var = receive_buffer_size / 3;
                receive_buffer_size = send_buffer_size / 3;
                send_buffer_size = tmp_var;

                mpInterfaceObjectManagerBins->FillSendBufferWithResults(send_buffer, send_buffer_size, min_distances);
                MPI_Sendrecv(send_buffer, send_buffer_size, MPI_DOUBLE, comm_partner, 0, receive_buffer,
                             receive_buffer_size, MPI_DOUBLE, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                mpInterfaceObjectManagerBins->FillSendBufferWithResults(send_buffer_int, send_buffer_size, pairing_indices);
                MPI_Sendrecv(send_buffer_int, send_buffer_size, MPI_INT, comm_partner, 0, receive_buffer_int,
                             receive_buffer_size, MPI_INT, comm_partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                mpInterfaceObjectManager->PostProcessReceivedResults(rCandidateManager, receive_buffer,
                        receive_buffer_int, comm_partner);

            } // if I am communicating in this loop (comm_partner != -1)
        } // loop communications

        delete [] send_buffer;
        delete [] send_buffer_int;
        delete [] receive_buffer;
        delete [] receive_buffer_int;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void PrepareSelection(CandidateManager& rCandidateManager,
                          bool& rLocalSearchRequired,
                          GraphType& rDomainsColoredGraph,
                          int& rMaxColors)
    {
        // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
        int* local_comm_list = new int[mCommSize]();
        int* local_memory_size_array = new int[mCommSize]();

        mpInterfaceObjectManager->PrepareMatching(rCandidateManager,
                local_comm_list,
                local_memory_size_array);

        if (local_comm_list[mCommRank] == 1)
        {
            rLocalSearchRequired = true;
        }

        // Compute Buffer Sizes for Mapping
        MapperUtilitiesMPI::ComputeMaxBufferSizes(local_memory_size_array,
                mSendBufferSize,
                mReceiveBufferSize,
                mCommRank, mCommSize);

        MapperUtilitiesMPI::ComputeColoringGraph(local_comm_list, mCommSize, rDomainsColoredGraph, rMaxColors);

        mDomainsColoredGraph = rDomainsColoredGraph; // save it for the mapping
        mMaxColors = rMaxColors;
        // Output the colored Graph
        if (mCommRank == 0 && mEchoLevel > 2)
        {
            MapperUtilitiesMPI::PrintGraph(rDomainsColoredGraph, rMaxColors);
        }
        delete [] local_comm_list;
        delete [] local_memory_size_array;
    }

    void SelectNeighbors(CandidateManager& rCandidateManager)
    {
        // // This function finds neighbors of the nodes in point_comm_manager in point_comm_manager_bins
        int* send_buffer = new int[mSendBufferSize];
        int* receive_buffer = new int[mReceiveBufferSize];
        // for debugging, i.e. to check the max buffer sizes
        int max_send_buffer_init_size = mSendBufferSize;
        int max_receive_buffer_init_size = mReceiveBufferSize;

        bool local_search_required = false;
        GraphType domains_colored_graph;
        int max_colors;
        PrepareSelection(rCandidateManager, local_search_required,
                         domains_colored_graph, max_colors);

        int send_buffer_size;
        int receive_buffer_size;

        // search in local partition
        if (local_search_required)   // search in own partition is required
        {
            mpInterfaceObjectManager->FillSendBufferWithMatchInformation(rCandidateManager, send_buffer,
                    send_buffer_size, mCommRank);
            mpInterfaceObjectManagerBins->ProcessMatchInformation(rCandidateManager, send_buffer,
                    send_buffer_size, mCommRank);
        }

        // search in remote partitions
        for (int i = 0; i < max_colors; ++i)   // loop over communication steps
        {
            int comm_partner = domains_colored_graph(mCommRank, i); // get the partner rank
            if (comm_partner != -1)
            {
                mpInterfaceObjectManager->FillSendBufferWithMatchInformation(rCandidateManager, send_buffer,
                        send_buffer_size, comm_partner);

                MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size,
                                                receive_buffer_size, max_send_buffer_init_size,
                                                max_receive_buffer_init_size, comm_partner);

                mpInterfaceObjectManagerBins->ProcessMatchInformation(rCandidateManager, receive_buffer,
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
