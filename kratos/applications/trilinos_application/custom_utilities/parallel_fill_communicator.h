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


// External includes


// Project includes
#include "includes/define.h"
#include "processes/graph_coloring_process.h"
#include "mpi.h"
#include "includes/mpi_communicator.h"

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

    void ComputeCommunicationPlan(ModelPart& rModelPart)
    {
        int root_id = 0;
        
        Communicator::Pointer pnew_comm = boost::make_shared< MPICommunicator >(&rModelPart.GetNodalSolutionStepVariablesList());
        rModelPart.SetCommunicator(pnew_comm);


        //get rank of current processor
        int my_rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Status status;

        //get number of processors
        int num_processors = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

        //generate (on each processor) a list of what the processor needs to receive (initialize it to -1 if no receive is needed)
        vector<int> my_receive_list_full(num_processors, -1);
        for (ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); it++)
        {
            int index = it->FastGetSolutionStepValue(PARTITION_INDEX);

            if(index < 0)
                KRATOS_THROW_ERROR(std::logic_error,"the partition index can not be lesser than one. something failed","");

            if(index != my_rank)
                my_receive_list_full[index] = 1;
        }

        std::vector<int> receive_list_compact;
        receive_list_compact.reserve(30);
        for (int i = 0; i < num_processors; i++)
        {
            if (my_receive_list_full[i] == 1)
                receive_list_compact.push_back(i);
        }

        //HERE WE COMMUNICATE HOW MANY "SEND" ARE NEEDED
        //communicate to node 0 the nneighbours
        int *number_or_receive_needed = NULL;
        int **neighbours = NULL;
        if (my_rank == root_id)
        {
            number_or_receive_needed = new int[num_processors];
            neighbours = new int*[num_processors];
        }

        // 	MPI_Barrier(MPI_COMM_WORLD);

        //communicate to node 0 the indices;
        int nrecv = receive_list_compact.size();
        MPI_Gather(&nrecv, 1, MPI_INT, number_or_receive_needed, 1, MPI_INT, root_id, MPI_COMM_WORLD);
        // 	MPI_Barrier(MPI_COMM_WORLD);



        //allocate memory on root
        if (my_rank == root_id)
        {
            for (int i = 0; i < num_processors; i++)
            {
                neighbours[i] = new int[number_or_receive_needed[i] ];
            }

//                for (int i = 0; i < num_processors; i++)
//                    std::cout << number_or_receive_needed[i] << std::endl;
        }

        //now gather on node 0 all of the neighbours to the partitions
        if (my_rank == root_id) //on node 0 we directly copy the data without mpi call
            for (int i = 0; i<int(receive_list_compact.size()); i++)
                neighbours[0][i] = receive_list_compact[i];

        for (int i = 1; i < num_processors; i++)
        {
            if (my_rank == root_id)
            {
                MPI_Recv(neighbours[i], number_or_receive_needed[i], MPI_INT, i, i, MPI_COMM_WORLD, &status);
                //KRATOS_WATCH(status);
            }
            else if (my_rank == i)
            {
                int nreceives = receive_list_compact.size();
                int* temp = new int[nreceives];
                for (unsigned int k = 0; k < receive_list_compact.size(); k++)
                    temp[k] = receive_list_compact[k];
                // 	      MPI_Send(receive_list_compact.data(), receive_list_compact.size(), MPI_INT, 0, i, MPI_COMM_WORLD);
                MPI_Send(temp, receive_list_compact.size(), MPI_INT, 0, i, MPI_COMM_WORLD);
                delete [] temp;
            }
            // 	  MPI_Barrier(MPI_COMM_WORLD);
        }


        // 	if(my_rank == root_id)
        // 	{
        // 	  for(unsigned int i=0; i<num_processors; i++)
        // 	  {
        // 	    for(unsigned int j=0; j<number_or_receive_needed[i]; j++)
        // 	      std::cout << neighbours[i][j] << " ";
        // 	    std::cout << std::endl;
        // 	  }
        // 	}




        //*************************************************************************************
        //here do the coloring - the scalar part should be improved quite a lot!!
        matrix<int> dense_colored_graph;
        matrix<int> dense_graph;
        int max_color_found = -1;
        if (my_rank == root_id)
        {
            dense_graph = scalar_matrix<int>(num_processors, num_processors, 0);

            // 	    KRATOS_WATCH(dense_graph);

            for (int i = 0; i < num_processors; i++)
            {
                for (int j = 0; j < number_or_receive_needed[i]; j++)
                {
                    int index1 = i;
                    int index2 = neighbours[i][j];

                    if(index1 == index2)
                        KRATOS_THROW_ERROR(std::logic_error,"trying to communicate with the node itself","");

                    if (index1 != index2)
                    {
                        dense_graph(index1, index2) = 1;
                        dense_graph(index2, index1) = 1;
                    }
                }
            }

            // 	    KRATOS_WATCH(dense_graph);

            int max_color = 2 * num_processors;
            GraphColoringProcess coloring_process(num_processors, dense_graph, dense_colored_graph, max_color);
            coloring_process.Execute();
            //KRATOS_WATCH(dense_colored_graph);
            //count max colors

            for (int i = 0; i< static_cast<int> (num_processors); i++)
                for (int j = 0; j<static_cast<int>( max_color); j++)
                    if (dense_colored_graph(i, j) != -1 && max_color_found < j) max_color_found = j;

            max_color_found += 1;


// 		//verify that the communication graph is correct
// 		for(unsigned int i=0; i<num_processors; i++)
// 			for(unsigned int j=0; j<max_color_found; j++)
// 			{
//
// 				int ij_entry = dense_colored_graph(i, j);
//
// 				//communication needed
// 				if( ij_entry != -1)
// 				{
// 					if(dense_colored_graph(ij_entry,j) != i)
// 						KRATOS_THROW_ERROR(std::logic_error,"communication is not symmetric - case A. Error!!","");
// 				}
// 				else
// 				{
// 					for(unsigned int k=0; k<num_processors;k++)
// 						if(dense_colored_graph(k, j) == i)
// 							KRATOS_THROW_ERROR(std::logic_error,"communication is not symmetric - case B. Error!!","");
// 				}
//
//
// 			}

        }



        //scatter the max_number_of_colors found
        int* aux = NULL;
        if (my_rank == root_id)
        {
            aux = new int[num_processors];
            for (int i = 0; i < num_processors; i++)
                aux[i] = max_color_found;
        }

        //here send an array of size max_color_found
        MPI_Scatter(aux, 1, MPI_INT, &max_color_found, 1, MPI_INT, root_id, MPI_COMM_WORLD);
        //KRATOS_WATCH(max_color_found);
        if (my_rank == root_id)
            delete [] aux;

        //now spread the colors of the communication to the processors.
        int* colors = new int[max_color_found ];
        int* send_colors = new int[max_color_found ];
        if (my_rank == root_id)
            for (int j = 0; j < max_color_found ; j++)
                colors[j] = dense_colored_graph(0, j);

        for (int i = 1; i < num_processors; i++)
        {
            if (my_rank == root_id)
            {
                for (int j = 0; j < max_color_found ; j++)
                {
                    send_colors[j] = dense_colored_graph(i, j);
//                        std::cout << send_colors[j] << " ";
                }
//                    std::cout << std::endl;

                MPI_Send(send_colors, max_color_found , MPI_INT, i, i, MPI_COMM_WORLD);

                //KRATOS_WATCH(status);
            }
            else if (my_rank == i)
            {
                MPI_Recv(colors, max_color_found , MPI_INT, 0, i, MPI_COMM_WORLD, &status);
            }
        }

        // 	for(unsigned int i=0; i<num_processors; i++)
        // 	  {
        // 	    if(my_rank == i)
        // 	    {
        // 	      std::cout << "processor i ="<< i << " colors -->" ;
        // 	      for(unsigned int j=0; j<max_color_found; j++)
        // 		std::cout << colors[j] << " ";
        // 	      std::cout << std::endl;
        // 	    }
        // 	    MPI_Barrier(MPI_COMM_WORLD);
        // 	  }

        //allocate space needed in the communicator
        rModelPart.GetCommunicator().SetNumberOfColors(max_color_found );
        rModelPart.GetCommunicator().NeighbourIndices().resize(max_color_found);

        for (int i = 0; i<max_color_found; i++)
        {
            rModelPart.GetCommunicator().LocalMesh(i).Nodes().clear();
            rModelPart.GetCommunicator().GhostMesh(i).Nodes().clear();
            rModelPart.GetCommunicator().InterfaceMesh(i).Nodes().clear();
        }

        //for each color fill the list of ghost and local nodes and the interface mesh
        for (int i = 0; i < max_color_found ; i++)
        {
            rModelPart.GetCommunicator().NeighbourIndices()[i] = colors[i];
            GenerateMeshes(colors[i], my_rank, i, rModelPart);
        }

        //now fill the list of all of the nodes to be communicated
        ModelPart::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().LocalMesh().Nodes();
        ModelPart::NodesContainerType& r_ghost_nodes = rModelPart.GetCommunicator().GhostMesh().Nodes();
        ModelPart::NodesContainerType& r_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh().Nodes();

        r_local_nodes.clear();
        r_ghost_nodes.clear();
        r_interface_nodes.clear();

        if(r_local_nodes.size() != 0)
            KRATOS_THROW_ERROR(std::logic_error, "local size can not be zero","")

            // filling all locals and all ghosts (independently of the color)
//            for (int i = 0; i < max_color_found + 1; i++)
//            {
//                 for (ModelPart::NodesContainerType::iterator it =  rModelPart.GetCommunicator().LocalMesh(i).NodesBegin();
//                                                              it != rModelPart.GetCommunicator().LocalMesh(i).NodesEnd(); it++)
//                    r_local_nodes.push_back(*(it.base()));
//
//                 for (ModelPart::NodesContainerType::iterator it =  rModelPart.GetCommunicator().GhostMesh(i).NodesBegin();
//                                                              it != rModelPart.GetCommunicator().GhostMesh(i).NodesEnd(); it++)
//                    r_ghost_nodes.push_back(*(it.base()));
//            }

            for (ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); it++)
            {
                int index = it->FastGetSolutionStepValue(PARTITION_INDEX);
                if (index == my_rank)
                    r_local_nodes.push_back(*(it.base()));
                else
                    r_ghost_nodes.push_back(*(it.base()));
            }

        //finally fill the interface - first add all of the local
        for (int i = 0; i < max_color_found ; i++)
        {
            ModelPart::NodesContainerType& r_interface_nodes_color = rModelPart.GetCommunicator().InterfaceMesh(i).Nodes();

            for (ModelPart::NodesContainerType::iterator it = r_interface_nodes_color.begin(); it != r_interface_nodes_color.end(); it++)
                r_interface_nodes.push_back(*(it.base()));
            // 	    for(ModelPart::NodesContainerType::ptr_iterator it = r_interface_nodes_color.ptr_begin(); it!=r_interface_nodes_color.ptr_end(); it++)
            // 	      r_interface_nodes.push_back(*it);
        }
        r_interface_nodes.Unique();
        r_local_nodes.Unique();
        r_ghost_nodes.Unique();

        //finally assign elements and conditions
        rModelPart.GetCommunicator().LocalMesh().Elements().clear();
        rModelPart.GetCommunicator().LocalMesh().Conditions().clear();
        rModelPart.GetCommunicator().LocalMesh().Elements() = rModelPart.Elements();
        rModelPart.GetCommunicator().LocalMesh().Conditions() = rModelPart.Conditions();





        //*******************************************
        //on each processor initialize the MPI communicator

        //deallocate memory
        delete [] colors;
        delete [] send_colors;
        if (my_rank == root_id)
        {
            delete [] number_or_receive_needed;
            for (int i = 0; i < num_processors; i++)
            {
                delete [] neighbours[i];
            }
            delete [] neighbours;
        }
        
        
        //now call the submodelpart
        for (auto i_sub_model_part = rModelPart.SubModelPartsBegin(); i_sub_model_part != rModelPart.SubModelPartsEnd(); i_sub_model_part++)
        {
            ComputeCommunicationPlan( *i_sub_model_part );
        }


    }




















    ///this function is designed to obtain on each of the processors, and on each of the colors,
    ///the list of ids this processor will need to receive from the "communicate_processor"
    ///and the list of ids it will need to send to that processor
    void GenerateMeshes(int communicate_processor, int my_rank, int color, ModelPart& rModelPart)
    {
        KRATOS_TRY

        if(communicate_processor == my_rank)
            KRATOS_THROW_ERROR(std::logic_error,"communicate_processor coincides with rank! this should not happen","");

        if (communicate_processor != -1)
        {
            int nnodes_to_send = -1;
            int nnodes_to_receive = -1;

            int* ids_to_receive = NULL;
            int* ids_to_send = NULL;

            ModelPart::NodesContainerType& kratos_nodes_to_receive = rModelPart.GetCommunicator().GhostMesh(color).Nodes();
            kratos_nodes_to_receive.clear();

            //loop on all the nodes and find how many nodes we need to receive from the other processors
            // 	    ModelPart::NodesContainerType kratos_nodes_to_receive;
            for (ModelPart::NodesContainerType::iterator it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); it++)
            {
                int index = it->FastGetSolutionStepValue(PARTITION_INDEX);
                if (index == communicate_processor)
                    kratos_nodes_to_receive.push_back(*(it.base()));
            }
            int temp  = kratos_nodes_to_receive.size();
            kratos_nodes_to_receive.Unique();
            if(temp != static_cast<int>(kratos_nodes_to_receive.size()))
                KRATOS_THROW_ERROR(std::logic_error,"the list of nodes to receive has repeated nodes","");

            nnodes_to_receive = kratos_nodes_to_receive.size();
            ids_to_receive = new int[ nnodes_to_receive ];

            int i = 0;
            for (ModelPart::NodesContainerType::iterator it = kratos_nodes_to_receive.begin(); it != kratos_nodes_to_receive.end(); it++)
                ids_to_receive[i++] = it->Id();

            //syncronize how many nodes need to be sent/received
            MPI_Status status;
            int send_tag = color;
            int receive_tag = color;
            MPI_Sendrecv(&nnodes_to_receive, 1, MPI_INT, communicate_processor, send_tag, &nnodes_to_send, 1, MPI_INT, communicate_processor, receive_tag, MPI_COMM_WORLD, &status);

            //get the list of nodes the other processor needs to be sent, and send the processor that we need to receive_list_compact
            ids_to_send = new int[ nnodes_to_send ];
            MPI_Sendrecv(ids_to_receive, nnodes_to_receive, MPI_INT, communicate_processor, send_tag, ids_to_send, nnodes_to_send, MPI_INT, communicate_processor, receive_tag, MPI_COMM_WORLD, &status);



            //allocate memory in the communicator
//                std::cout << "my_rank = " << my_rank << std::endl;
//                std::cout << "color = " << color << " ids to receive from " << communicate_processor << " -->";
//                for (int i = 0; i < nnodes_to_receive; i++)
//                    std::cout << ids_to_receive[i] << " ";
//                std::cout << std::endl;
//                std::cout << "color = " << color << " ids to send to " << communicate_processor << " -->";
//                for (int i = 0; i < nnodes_to_send; i++)
//                    std::cout << ids_to_send[i] << " ";
//                std::cout << std::endl;

            ModelPart::NodesContainerType& r_local_nodes = rModelPart.GetCommunicator().LocalMesh(color).Nodes();
            r_local_nodes.clear();

//                KRATOS_WATCH("ln111111");
            //fill the list of nodes to be sent
            for (int i = 0; i < nnodes_to_send; i++)
            {
                r_local_nodes.push_back(rModelPart.Nodes()(ids_to_send[i]));
            }

            for (ModelPart::NodesContainerType::iterator it = r_local_nodes.begin(); it != r_local_nodes.end(); it++)
                if(it->FastGetSolutionStepValue(PARTITION_INDEX)!=my_rank)
                    KRATOS_THROW_ERROR(std::logic_error,"a node in the local mesh is trying to communicate to the wrong partition","");


//                KRATOS_WATCH("ln222222");
            r_local_nodes.Unique();
            if(static_cast<int>(r_local_nodes.size()) != nnodes_to_send)
                KRATOS_THROW_ERROR(std::logic_error,"impossible situation. Something wrong happend","");

            //add local and ghost to the interface mesh
            ModelPart::NodesContainerType& r_interface_nodes = rModelPart.GetCommunicator().InterfaceMesh(color).Nodes();
            r_interface_nodes.clear();

            // 	    for(ModelPart::NodesContainerType::ptr_iterator it = kratos_nodes_to_receive.ptr_begin(); it!=kratos_nodes_to_receive.ptr_end(); it++)
            // 	      r_interface_nodes.push_back(*it);
            //
            // 	    for(ModelPart::NodesContainerType::ptr_iterator it = r_local_nodes.ptr_begin(); it!=r_local_nodes.ptr_end(); it++)
            // 	      r_interface_nodes.push_back(*it);

            for (ModelPart::NodesContainerType::iterator it = kratos_nodes_to_receive.begin(); it != kratos_nodes_to_receive.end(); it++)
                r_interface_nodes.push_back(*(it.base()));

            for (ModelPart::NodesContainerType::iterator it = r_local_nodes.begin(); it != r_local_nodes.end(); it++)
                r_interface_nodes.push_back(*(it.base()));

            int size = r_interface_nodes.size();
            r_interface_nodes.Unique();
            if(size != static_cast<int>(r_interface_nodes.size()))
                KRATOS_THROW_ERROR(std::logic_error,"something went wrong in the interface nodes","");


            delete [] ids_to_receive;
            delete [] ids_to_send;
        }

        KRATOS_CATCH("");
    }

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


