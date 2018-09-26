//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "processes/graph_coloring_process.h"
#include "mpi.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos
{

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

/// This function recomputes the communication plan for MPI

/** The objective of this class is to read the mesh owned by each node in a distributed context
 * and to fill the communication plan (coloring) so to allow the communication to be performed correctly
 * It fills the Ghost and Local lists and performs the coloring, then it updates the MPI communicator
 */
class ParallelFillCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParallelFillCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(ParallelFillCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    ParallelFillCommunicator(ModelPart& r_model_part)
        : mrBaseModelPart(r_model_part)
    {
    }

    /// Destructor.

    virtual ~ParallelFillCommunicator()
    {
    }

    void Execute()
    {
        KRATOS_TRY

        //use epetra to compute the communication plan
        ComputeCommunicationPlan(mrBaseModelPart);


        KRATOS_CATCH("");
    }

    ///************************************************************************************************
    ///************************************************************************************************
    ///function to print DETAILED mesh information. WARNING: to be used for debugging only as many informations
    ///are plotted

    void PrintDebugInfo()
    {
        PrintModelPartDebugInfo(mrBaseModelPart);
    }

    void PrintModelPartDebugInfo(ModelPart& rModelPart)
    {
        KRATOS_TRY

        std::cout.flush();
        MPI_Barrier(MPI_COMM_WORLD);



        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        //get number of processors
        int num_processors = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

        for (int i = 0; i < num_processors; i++)
        {
            if (rank == i)
            {

                std::cout << " *************************************** " << std::endl;

                std::cout << " proc = " << rank << "communication colors " << rModelPart.GetCommunicator().NeighbourIndices() << std::endl;

                //print ghost mesh
                std::cout << " proc = " << rank << " ghost mesh" << std::endl;
                for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().GhostMesh().NodesBegin();
                        it != rModelPart.GetCommunicator().GhostMesh().NodesEnd();
                        it++)
                {
                    if(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank)
                        KRATOS_THROW_ERROR(std::logic_error,"error partition index can not be = to rank for ghost nodes","")
                        std::cout << it->Id() << " " ;
                }
                std::cout << std::endl;

                //print local mesh
                std::cout << " proc = " << rank << " local mesh" << std::endl;
                for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh().NodesBegin();
                        it != rModelPart.GetCommunicator().LocalMesh().NodesEnd();
                        it++)
                {
                    if(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                        KRATOS_THROW_ERROR(std::logic_error,"error partition index can not be != from rank for local nodes","")
                        std::cout << it->Id() << " " ;
                }
                std::cout << std::endl;

                //print interface mesh
                std::cout << " proc = " << rank << " interface mesh" << std::endl;
                for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().InterfaceMesh().NodesBegin();
                        it != rModelPart.GetCommunicator().InterfaceMesh().NodesEnd();
                        it++)
                {
                    std::cout << it->Id() << " " ;
                }
                std::cout << std::endl;

                //now print everything color by color

                int destination = 0;
                std::cout << "NeighbourIndices " ;
                const vector<int>& neighbours_indices = rModelPart.GetCommunicator().NeighbourIndices();
                for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                    std::cout << neighbours_indices[i_color] << " " ;
                std::cout << std::endl;
                for (unsigned int i_color = 0; i_color < neighbours_indices.size(); i_color++)
                {
                    std::cout << "color = " << i_color << std::endl;
                    if ((destination = neighbours_indices[i_color]) >= 0)
                    {
                        std::cout << "ghost mesh for color --> " << i_color << std::endl;
                        for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().GhostMesh(i_color).NodesBegin();
                                it != rModelPart.GetCommunicator().GhostMesh(i_color).NodesEnd();
                                it++)
                        {
                            if(it->FastGetSolutionStepValue(PARTITION_INDEX)==rank)
                                KRATOS_THROW_ERROR(std::logic_error,"error partition index can not be = to rank for ghost nodes","")
                                std::cout << it->Id() << " " ;
                        }

                        std::cout << "finished printing ghost mesh for color --> " << i_color<< std::endl;

                        std::cout << "local mesh for color --> " << i_color << std::endl;
                        for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh(i_color).NodesBegin();
                                it != rModelPart.GetCommunicator().LocalMesh(i_color).NodesEnd();
                                it++)
                        {
                            if(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                                KRATOS_THROW_ERROR(std::logic_error,"error partition index can not be != from rank for local nodes","")
                                std::cout << it->Id() << " " ;
                        }
                        std::cout << "finished printing local mesh for color --> " << i_color<< std::endl;

                        std::cout << "interface mesh for color --> " << i_color << std::endl;
                        for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().InterfaceMesh(i_color).NodesBegin();
                                it != rModelPart.GetCommunicator().InterfaceMesh(i_color).NodesEnd();
                                it++)
                        {
                            std::cout << it->Id() << " " ;
                        }
                        std::cout << "finished printing interface mesh for color --> " << i_color<< std::endl;
                    }
                    else
                    {
                        if(rModelPart.GetCommunicator().GhostMesh(i_color).Nodes().size()!=0)
                        {
                            std::cout << "rank = " << rank << " color = " << i_color << std::endl;
                            KRATOS_THROW_ERROR(std::logic_error,"nodes found in ghost mesh when communication is not expected","")
                        }
                        if(rModelPart.GetCommunicator().LocalMesh(i_color).Nodes().size()!=0)
                        {
                            std::cout << "local mesh for color --> " << i_color << "*********************************" <<  std::endl;
                            for (ModelPart::NodesContainerType::iterator it = rModelPart.GetCommunicator().LocalMesh(i_color).NodesBegin();
                                    it != rModelPart.GetCommunicator().LocalMesh(i_color).NodesEnd();
                                    it++)
                            {
                                if(it->FastGetSolutionStepValue(PARTITION_INDEX)!=rank)
                                    KRATOS_THROW_ERROR(std::logic_error,"error partition index can not be != from rank for local nodes","")
                                    std::cout << it->Id() << " " << it->FastGetSolutionStepValue(PARTITION_INDEX) << std::endl ;
                            }
                            std::cout << "finished printing local mesh for color --> " << i_color<< std::endl;
                            std::cout << "nodes found in local mesh when communication is not expected" << std::endl;
                            KRATOS_THROW_ERROR(std::logic_error,"nodes found in local mesh when communication is not expected","")
                        }
                        if(rModelPart.GetCommunicator().InterfaceMesh(i_color).Nodes().size()!=0)
                            KRATOS_THROW_ERROR(std::logic_error,"nodes found in interface mesh when communication is not expected","")
                        }
                }

                std::cout << "finished printing proc -> " << rank << "*********************" << std::endl;
                std::cout << std::endl;
                std::cout.flush();

            }

            MPI_Barrier(MPI_COMM_WORLD);
        }
        KRATOS_CATCH("");

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
        buffer << "ParallelFillCommunicator";
        return buffer.str();
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ParallelFillCommunicator" << std::endl;
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
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

    void ComputeCommunicationPlan(ModelPart& rModelPart)
    {
        constexpr unsigned root_id = 0;

        Communicator::Pointer pnew_comm = Kratos::make_shared< MPICommunicator >(&rModelPart.GetNodalSolutionStepVariablesList());
        rModelPart.SetCommunicator(pnew_comm);

        // Get rank of current processor.
        int mpi_rank;
        KRATOS_ERROR_IF_NOT(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank) == MPI_SUCCESS) << "MPI_Comm_rank() failed." << std::endl;
        const unsigned my_rank = mpi_rank;
        MPI_Status status;

        // Get number of processors.
        int mpi_num_processors;
        KRATOS_ERROR_IF_NOT(MPI_Comm_size(MPI_COMM_WORLD, &mpi_num_processors) == MPI_SUCCESS) << "MPI_Comm_size() failed." << std::endl;
        const unsigned num_processors = mpi_num_processors;

        // Find all ghost nodes on this process and mark the corresponding neighbour process for communication.
        vector<bool> receive_from_neighbour(num_processors, false);
        for (const auto& rNode : rModelPart.Nodes())
        {
            const unsigned partition_index = rNode.FastGetSolutionStepValue(PARTITION_INDEX);
            KRATOS_ERROR_IF(partition_index >= num_processors) << "The partition index is out of range. Invalid model part." << std::endl;
            if(partition_index != my_rank)
                receive_from_neighbour[partition_index] = true;
        }

        // Make a list of my receive process ids.
        std::vector<unsigned> my_receive_neighbours;
        my_receive_neighbours.reserve(30);
        for (unsigned p_id = 0; p_id < num_processors; ++p_id)
        {
            if (receive_from_neighbour[p_id])
                my_receive_neighbours.push_back(p_id);
        }

        // Initialize arrays for all neighbour id lists on root process.
        std::vector<unsigned> number_of_receive_neighbours;
        std::vector<std::vector<unsigned>> receive_neighbours;
        if (my_rank == root_id)
        {
            number_of_receive_neighbours.resize(num_processors);
            receive_neighbours.resize(num_processors);
        }
        {
            unsigned send_buf = my_receive_neighbours.size();
            int ierr = MPI_Gather(&send_buf, 1, MPI_UNSIGNED, number_of_receive_neighbours.data(), 1, MPI_UNSIGNED, root_id, MPI_COMM_WORLD);
            KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Gather() failed." << std::endl;
        }
        if (my_rank == root_id)
            for (unsigned p_id = 0; p_id < num_processors; ++p_id)
                receive_neighbours[p_id].resize(number_of_receive_neighbours[p_id]);

        // Fill the neighbour id lists of the partitions on root.
        if (my_rank == root_id) // On root we directly copy the data without calling MPI.
            std::copy(my_receive_neighbours.begin(), my_receive_neighbours.end(), receive_neighbours[root_id].begin());
        // Gather the remaining id lists to root.
        for (unsigned p_id = 1; p_id < num_processors; ++p_id)
            if (my_rank == root_id)
            {
                int ierr = MPI_Recv(receive_neighbours[p_id].data(), number_of_receive_neighbours[p_id], MPI_UNSIGNED,
                                    p_id, p_id, MPI_COMM_WORLD, &status);
                KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Recv() failed." << std::endl;
            }
            else if (my_rank == p_id)
            {
                int ierr = MPI_Send(my_receive_neighbours.data(), my_receive_neighbours.size(), MPI_UNSIGNED,
                                    root_id, p_id, MPI_COMM_WORLD);
                KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Send() failed." << std::endl;
            }

        // Create the colored graph for communication.
        DenseMatrix<int> domains_colored_graph;
        int max_color_found = -1;
        if (my_rank == root_id)
        {
            ///@TODO for large problems, this should use a compressed matrix.
            DenseMatrix<int> domains_graph = ZeroMatrix(num_processors, num_processors);
            for (unsigned index1 = 0; index1 < num_processors; ++index1)
                for (unsigned index2 : receive_neighbours[index1])
                {
                    KRATOS_ERROR_IF(index1 == index2) << "Trying to communicate with the node itself." << std::endl;
                    domains_graph(index1, index2) = 1;
                    domains_graph(index2, index1) = 1;
                }

            // max_color is overwritten by the GraphColoringProcess.
            int max_color = 2 * num_processors; // Max. number of one-directional communications (this has no effect).
            GraphColoringProcess coloring_process(num_processors, domains_graph, domains_colored_graph, max_color);
            coloring_process.Execute();
            // Count max colors.
            for (unsigned p_id = 0; p_id < num_processors; ++p_id)
                for (int j = 0; j < max_color; ++j)
                    if (domains_colored_graph(p_id, j) != -1 && max_color_found < j) max_color_found = j;

            max_color_found += 1;
        }

        // Scatter max_color_found.
        if (my_rank == root_id)
        {
            std::vector<int> send_buf(num_processors, max_color_found);
            int ierr = MPI_Scatter(send_buf.data(), 1, MPI_INT, &max_color_found, 1, MPI_INT, root_id, MPI_COMM_WORLD);
            KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Scatter() failed." << std::endl;
        }
        else
        {
            int ierr = MPI_Scatter(nullptr, 1, MPI_INT, &max_color_found, 1, MPI_INT, root_id, MPI_COMM_WORLD);
            KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Scatter() failed." << std::endl;
        }

        // Now send the colors of the communication to the processors.
        std::vector<int> colors(max_color_found);
        if (my_rank == root_id) // On root we directly copy the data.
            for (int j = 0; j < max_color_found; ++j)
                colors[j] = domains_colored_graph(root_id, j);
        // Send the remaining color patterns to processes.
        std::vector<int> send_colors(max_color_found);
        for (unsigned p_id = 1; p_id < num_processors; ++p_id)
        {
            if (my_rank == root_id)
            {
                for (int j = 0; j < max_color_found; ++j)
                    send_colors[j] = domains_colored_graph(p_id, j);
                int ierr = MPI_Send(send_colors.data(), max_color_found, MPI_INT, p_id, p_id, MPI_COMM_WORLD);
                KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Send() failed." << std::endl;
            }
            else if (my_rank == p_id)
            {
                int ierr = MPI_Recv(colors.data(), max_color_found, MPI_INT, root_id, p_id, MPI_COMM_WORLD, &status);
                KRATOS_ERROR_IF_NOT(ierr == MPI_SUCCESS) << "MPI_Recv() failed." << std::endl;
            }
        }

        InitializeParallelCommunicationMeshes(rModelPart, colors, my_rank);
    }

    /// Initialize the communicator's ghost, local and interface meshes for all communication pairs (colors).
    void InitializeParallelCommunicationMeshes(ModelPart& rModelPart,
                                               const std::vector<int>& rColors,
                                               unsigned MyRank)
    {
        KRATOS_TRY;
        // Allocate space needed in the communicator.
        rModelPart.GetCommunicator().SetNumberOfColors(rColors.size());
        rModelPart.GetCommunicator().NeighbourIndices().resize(rColors.size());
        for (unsigned color = 0; color < rColors.size(); ++color)
        {
            rModelPart.GetCommunicator().LocalMesh(color).Nodes().clear();
            rModelPart.GetCommunicator().GhostMesh(color).Nodes().clear();
            rModelPart.GetCommunicator().InterfaceMesh(color).Nodes().clear();
        }

        // For each color fill the list of ghost and local nodes and the
        // interface
        // mesh.
        for (unsigned color = 0; color < rColors.size(); ++color)
        {
            rModelPart.GetCommunicator().NeighbourIndices()[color] = rColors[color];
            GenerateMeshes(rColors[color], MyRank, color, rModelPart);
        }

        // Fill the list of all of the nodes to be communicated.
        ModelPart::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().LocalMesh().Nodes();
        ModelPart::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().Nodes();
        ModelPart::NodesContainerType& r_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().Nodes();
        r_local_nodes.clear();
        r_ghost_nodes.clear();
        r_interface_nodes.clear();

        KRATOS_ERROR_IF(r_local_nodes.size() != 0)
            << "Local size can't be zero." << std::endl;

        // Fill nodes for LocalMesh and GhostMesh.
        for (auto it_node = rModelPart.NodesBegin(); it_node != rModelPart.NodesEnd(); ++it_node)
        {
            const unsigned index = it_node->FastGetSolutionStepValue(PARTITION_INDEX);
            if (index == MyRank)
                r_local_nodes.push_back(*(it_node.base()));
            else
                r_ghost_nodes.push_back(*(it_node.base()));
        }

        // Fill nodes for the InterfaceMesh.
        for (ModelPart::MeshType& r_color_interface_mesh : rModelPart.GetCommunicator().InterfaceMeshes())
        {
            ModelPart::NodesContainerType& r_color_interface_nodes =
                r_color_interface_mesh.Nodes();
            for (auto it = r_color_interface_nodes.begin(); it != r_color_interface_nodes.end(); ++it)
                r_interface_nodes.push_back(*(it.base()));
        }
        r_interface_nodes.Unique();
        r_local_nodes.Unique();
        r_ghost_nodes.Unique();

        // Assign elements and conditions for LocalMesh.
        rModelPart.GetCommunicator().LocalMesh().Elements().clear();
        rModelPart.GetCommunicator().LocalMesh().Conditions().clear();
        rModelPart.GetCommunicator().LocalMesh().Elements() = rModelPart.Elements();
        rModelPart.GetCommunicator().LocalMesh().Conditions() = rModelPart.Conditions();

        // Call the sub model part.
        for (ModelPart& r_sub_model_part : rModelPart.SubModelParts())
            ComputeCommunicationPlan(r_sub_model_part);

        KRATOS_CATCH("");
    }

    /// Generate the ghost, local and interface meshes for processes of a communication pair (color).
    void GenerateMeshes(int NeighbourPID, int MyPID, unsigned Color, ModelPart& rModelPart)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(NeighbourPID == MyPID)
            << "Neighbour process coincides with rank! this should not happen."
            << std::endl;

        if (NeighbourPID == -1) // Don't communicate with this process.
            return;

        ModelPart::NodesContainerType& r_ghost_nodes =
            rModelPart.GetCommunicator().GhostMesh(Color).Nodes();
        r_ghost_nodes.clear();

        // Fill nodes for GhostMesh(Color).
        for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            const int index = it->FastGetSolutionStepValue(PARTITION_INDEX);
            if (index == NeighbourPID)
                r_ghost_nodes.push_back(*(it.base()));
        }
        unsigned num_ghost_nodes = r_ghost_nodes.size();
        r_ghost_nodes.Unique();
        KRATOS_ERROR_IF(num_ghost_nodes != r_ghost_nodes.size())
            << "The list of nodes to receive has repeated nodes." << std::endl;

        std::vector<int> ids_to_receive(r_ghost_nodes.size());

        { // Fill receive ids (ids of ghost nodes).
            int i = 0;
            for (const ModelPart::NodeType& rNode : r_ghost_nodes)
                ids_to_receive[i++] = rNode.Id();
        }

        std::vector<int> ids_to_send;
        { // Syncronize how many nodes need to be sent/received.
            MPI_Status status;
            int send_tag = Color;
            int receive_tag = Color;
            unsigned send_buf = ids_to_receive.size();
            unsigned recv_buf;
            MPI_Sendrecv(&send_buf, 1, MPI_UNSIGNED, NeighbourPID, send_tag,
                         &recv_buf, 1, MPI_UNSIGNED, NeighbourPID, receive_tag,
                         MPI_COMM_WORLD, &status);
            ids_to_send.resize(recv_buf);
        }

        { // Send/receive node ids.
            MPI_Status status;
            int send_tag = Color;
            int receive_tag = Color;
            MPI_Sendrecv(ids_to_receive.data(), ids_to_receive.size(), MPI_INT, NeighbourPID,
                         send_tag, ids_to_send.data(), ids_to_send.size(), MPI_INT,
                         NeighbourPID, receive_tag, MPI_COMM_WORLD, &status);
        }

        // Fill nodes for LocalMesh(Color).
        ModelPart::NodesContainerType& r_local_nodes =
            rModelPart.GetCommunicator().LocalMesh(Color).Nodes();
        r_local_nodes.clear();
        for (int id : ids_to_send)
            r_local_nodes.push_back(rModelPart.Nodes()(id));

        for (const ModelPart::NodeType& r_node : r_local_nodes)
            KRATOS_ERROR_IF(r_node.FastGetSolutionStepValue(PARTITION_INDEX) != MyPID) << "A node in the local mesh is trying to communicate to the wrong partition."
                                                                                       << std::endl;

        r_local_nodes.Unique();
        KRATOS_ERROR_IF(r_local_nodes.size() != ids_to_send.size())
            << "Impossible situation. Something went wrong." << std::endl;

        // Fill InterfaceMesh(Color) with local and ghost nodes.
        ModelPart::NodesContainerType& r_interface_nodes =
            rModelPart.GetCommunicator().InterfaceMesh(Color).Nodes();
        r_interface_nodes.clear();

        for (auto it = r_ghost_nodes.begin(); it != r_ghost_nodes.end(); ++it)
            r_interface_nodes.push_back(*(it.base()));

        for (auto it = r_local_nodes.begin(); it != r_local_nodes.end(); it++)
            r_interface_nodes.push_back(*(it.base()));

        unsigned num_interface_nodes = r_interface_nodes.size();
        r_interface_nodes.Unique();
        KRATOS_ERROR_IF(num_interface_nodes != r_interface_nodes.size())
            << "Something went wrong in the interface nodes." << std::endl;

        KRATOS_CATCH("");
    }

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
    ModelPart& mrBaseModelPart;


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

    ParallelFillCommunicator & operator=(ParallelFillCommunicator const& rOther)
    {
        return *this;
    }

    /// Copy constructor.

    ParallelFillCommunicator(ParallelFillCommunicator const& rOther) = delete;
    //: rModelPart(rOther.mrBaseModelPart)
    //{
    //}


    ///@}

}; // Class ParallelFillCommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

inline std::istream & operator >>(std::istream& rIStream,
                                  ParallelFillCommunicator& rThis)
{
    return rIStream;
}

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const ParallelFillCommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_PARALLEL_FILL_COMMUNICATOR_H_INCLUDED  defined


