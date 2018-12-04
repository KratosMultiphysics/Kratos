//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//

#ifndef KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_IN_MEMORY_PROCESS_H
#define KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_IN_MEMORY_PROCESS_H

#ifdef KRATOS_USE_METIS_5
  #include "metis.h"
#else
  #include <parmetis.h>
#endif

// Project includes
#include "includes/define.h"
#include "includes/io.h"
#include "includes/model_part_io.h"

#include "processes/process.h"
#include "processes/graph_coloring_process.h"

#include "custom_processes/metis_divide_input_to_partitions_process.h"

// This one needs mpi
#include "mpi.h"

#ifndef KRATOS_USE_METIS_5
  extern "C" {
  extern void METIS_PartGraphKway(int*,  //int* n
                                  int*,  //idxtype* xadj
                                  int*,  //idxtype* adjcncy
                                  int*,  //idxtype* vwgt
                                  int*,  //idxtype* adjwgt
                                  int*,  //int* wgtflag
                                  int*,  //int* numflag
                                  int*,  //int* nparts
                                  int*,  //int* options
                                  int*,  //int* edgecut
                                  int*); //indxtype* part
  }
#endif

namespace Kratos
{
///@addtogroup MetisApplication
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

/// Call Metis to divide an heterogeneous mesh, by partitioning its nodal graph.
class MetisDivideHeterogeneousInputInMemoryProcess : public MetisDivideHeterogeneousInputProcess
{
public:
    ///@name Type Definitions
    ///@{

    #ifdef KRATOS_USE_METIS_5
      typedef idx_t idxtype;
    #else
      typedef int idxtype;
    #endif

    /// Pointer definition of MetisDivideHeterogeneousInputInMemoryProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisDivideHeterogeneousInputInMemoryProcess);

    typedef MetisDivideHeterogeneousInputProcess BaseType;

    typedef std::size_t SizeType;
    typedef std::size_t IndexType;
    typedef matrix<int> GraphType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetisDivideHeterogeneousInputInMemoryProcess(IO& rIO, ModelPartIO& rSerialIO, SizeType NumberOfPartitions, int Dimension = 3, int Verbosity = 0, bool SynchronizeConditions = false):
        BaseType(rIO,NumberOfPartitions,Dimension,Verbosity,SynchronizeConditions), mrSerialIO(rSerialIO)
    {
    }

    /// Copy constructor.
    MetisDivideHeterogeneousInputInMemoryProcess(MetisDivideHeterogeneousInputInMemoryProcess const& rOther):
        BaseType(rOther.mrIO,rOther.mNumberOfPartitions,rOther.mDimension,rOther.mVerbosity,rOther.mSynchronizeConditions), mrSerialIO(rOther.mrSerialIO)
    {
    }

    /// Destructor.
    virtual ~MetisDivideHeterogeneousInputInMemoryProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        this->Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    /// Generate a partition using Metis.
    /** Partitioned input is written as <problem name>_<mpi rank>.mdpa
     */
    void Execute() override
    {
        int mpi_rank;
        int mpi_size;

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        int * msgSendSize = new int[mpi_size];
        int * msgRecvSize = new int[mpi_size];

        const char ** mpi_send_buffer = new const char * [mpi_size];
        char ** mpi_recv_buffer = new char * [mpi_size];
        std::string * str = new std::string[mpi_size];

        // Set size
        for(int i = 0; i < mpi_size; i++) {
            msgSendSize[i] = 0;
            msgRecvSize[i] = 0;
        }

        // Transfer Streams
        Kratos::shared_ptr<std::iostream> * streams = new Kratos::shared_ptr<std::iostream>[mpi_size];
        std::stringbuf * stringbufs = new std::stringbuf[mpi_size];

        for(auto i = 0; i < mpi_size; i++) {
          streams[i] = Kratos::shared_ptr<std::iostream>(new std::iostream(&stringbufs[i]));
        }

        // Main process calculates the partitions and writes the result into temporal streams
        if(mpi_rank == 0) {
            // Read nodal graph from input

            IO::ConnectivitiesContainerType KratosFormatNodeConnectivities;

            SizeType NumNodes = BaseType::mrIO.ReadNodalGraph(KratosFormatNodeConnectivities);

            // Write connectivity data in CSR format
            idxtype* NodeIndices = 0;
            idxtype* NodeConnectivities = 0;

            ConvertKratosToCSRFormat(KratosFormatNodeConnectivities, &NodeIndices, &NodeConnectivities);

            std::vector<idxtype> NodePartition;
            PartitionNodes(NumNodes,NodeIndices,NodeConnectivities,NodePartition);

            // Free some memory we no longer need
            delete [] NodeIndices;
            delete [] NodeConnectivities;

            // Partition elements
            IO::ConnectivitiesContainerType ElementConnectivities;
            SizeType NumElements =  BaseType::mrIO.ReadElementsConnectivities(ElementConnectivities);
            if (NumElements != ElementConnectivities.size())
            {
                std::stringstream Msg;
                Msg << std::endl;
                Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
                Msg << "Read " << NumElements << " elements, but element list has " << ElementConnectivities.size() << " entries." << std::endl;
                Msg << "Elements are most likely not correlatively numbered." << std::endl;

                KRATOS_ERROR << Msg.str();
            }

            std::vector<idxtype> ElementPartition;

            if (mSynchronizeConditions)
                PartitionElementsSynchronous(NodePartition,ElementConnectivities,ElementPartition);
            else
                PartitionMesh(NodePartition,ElementConnectivities,ElementPartition);

            // Partition conditions
            IO::ConnectivitiesContainerType ConditionConnectivities;
            SizeType NumConditions = BaseType::mrIO.ReadConditionsConnectivities(ConditionConnectivities);
            if (NumConditions != ConditionConnectivities.size())
            {
                std::stringstream Msg;
                Msg << std::endl;
                Msg << "ERROR in MetisDivideHeterogenousInputProcess:" << std::endl;
                Msg << "Read " << NumConditions << " conditions, but condition list has " << ConditionConnectivities.size() << " entries." << std::endl;
                Msg << "Conditions are most likely not correlatively numbered." << std::endl;

                KRATOS_ERROR << Msg.str();
            }

            std::vector<idxtype> ConditionPartition;

            if (mSynchronizeConditions)
                PartitionConditionsSynchronous(NodePartition,ElementPartition,ConditionConnectivities,ElementConnectivities,ConditionPartition);
            else
                PartitionMesh(NodePartition,ConditionConnectivities,ConditionPartition);

            // Detect hanging nodes (nodes that belong to a partition where no local elements have them) and send them to another partition.
            // Hanging nodes should be avoided, as they can cause problems when setting the Dofs
            RedistributeHangingNodes(NodePartition,ElementPartition,ElementConnectivities,ConditionPartition,ConditionConnectivities);

            // Coloring
            GraphType DomainGraph = zero_matrix<int>(mNumberOfPartitions);
            CalculateDomainsGraph(DomainGraph,NumElements,ElementConnectivities,NodePartition,ElementPartition);
            CalculateDomainsGraph(DomainGraph,NumConditions,ConditionConnectivities,NodePartition,ConditionPartition);

            int NumColors;
            GraphType ColoredDomainGraph;
            GraphColoringProcess(mNumberOfPartitions,DomainGraph,ColoredDomainGraph,NumColors).Execute();

            if (mVerbosity > 0) {
                KRATOS_INFO("NumColors") << NumColors << std::endl;
            }

            if (mVerbosity > 2) {
                KRATOS_INFO("ColoredDomainGraph") << ColoredDomainGraph << std::endl;
            }

            // Write partition info into separate input files
            IO::PartitionIndicesContainerType nodes_all_partitions;
            IO::PartitionIndicesContainerType elements_all_partitions;
            IO::PartitionIndicesContainerType conditions_all_partitions;

            // Create lists containing all nodes/elements/conditions known to each partition
            DividingNodes(nodes_all_partitions, ElementConnectivities, ConditionConnectivities, NodePartition, ElementPartition, ConditionPartition);
            DividingElements(elements_all_partitions, ElementPartition);
            DividingConditions(conditions_all_partitions, ConditionPartition);

            if (mVerbosity > 1) {
                std::cout << "Final list of nodes known by each partition" << std::endl;
                for(SizeType i = 0 ; i < NumNodes ; i++) {
                    std::cout << "Node #" << i+1 << "->";
                    for(std::vector<std::size_t>::iterator j = nodes_all_partitions[i].begin() ; j != nodes_all_partitions[i].end() ; j++) {
                        std::cout << *j << ",";
                    }
                    std::cout << std::endl;
                }
            }

            IO::PartitionIndicesType io_nodes_partitions(NodePartition.begin(), NodePartition.end());
            IO::PartitionIndicesType io_elements_partitions(ElementPartition.begin(), ElementPartition.end());
            IO::PartitionIndicesType io_conditions_partitions(ConditionPartition.begin(), ConditionPartition.end());

            // Write files
            mrIO.DivideInputToPartitions(
                streams, mNumberOfPartitions, ColoredDomainGraph,
                io_nodes_partitions, io_elements_partitions, io_conditions_partitions,
                nodes_all_partitions, elements_all_partitions, conditions_all_partitions
            );
        }

        // Calculate the message and prepare the buffers
        if(mpi_rank == 0) {
            for(auto i = 0; i < mpi_size; i++) {
                str[i] = stringbufs[i].str();
                msgSendSize[i] = str[i].size();
                mpi_send_buffer[i] = str[i].c_str();
            }
        }

        // Send the message size to all processes
        MPI_Scatter(msgSendSize,1,MPI_INT,&msgRecvSize[mpi_rank],1,MPI_INT,0,MPI_COMM_WORLD);

        // Calculate the number of events:
        auto NumberOfCommunicationEvents = 1 + mpi_size * !mpi_rank;
        auto NumberOfCommunicationEventsIndex = 0;

        // Prepare the communication events
        MPI_Request * reqs = new MPI_Request[NumberOfCommunicationEvents];
        MPI_Status * stats = new MPI_Status[NumberOfCommunicationEvents];

        // Set up all receive and send events
        if( mpi_rank == 0) {
            for(auto i = 0; i < mpi_size; i++) {
                char* aux_char = const_cast<char*>(mpi_send_buffer[i]);
                MPI_Isend(aux_char,msgSendSize[i],MPI_CHAR,i,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);
            }
        }

        // Recieve the buffers
        mpi_recv_buffer[mpi_rank] = (char *)malloc(sizeof(char) * msgRecvSize[mpi_rank]);
        MPI_Irecv(mpi_recv_buffer[mpi_rank],msgRecvSize[mpi_rank],MPI_CHAR,0,0,MPI_COMM_WORLD,&reqs[NumberOfCommunicationEventsIndex++]);

        // Wait untill all communications finish
        if( MPI_Waitall(NumberOfCommunicationEvents, reqs, stats) != MPI_SUCCESS ) {
            KRATOS_ERROR << "Error in metis_partition_mem" << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_rank != 0) {
            streams[mpi_rank]->write(mpi_recv_buffer[mpi_rank], msgRecvSize[mpi_rank]);
        }

#ifdef KRATOS_DEBUG
        // Print the partitions in debug mode
        std::ofstream debug_ofstream("debug_modelpart_"+std::to_string(mpi_rank)+".mdpa");
        debug_ofstream << stringbufs[mpi_rank].str() << std::endl;
#endif

        // TODO: Try to come up with a better way to change the buffer.
        mrSerialIO.SwapStreamSource(streams[mpi_rank]);

        // Free buffers
        free(mpi_recv_buffer[mpi_rank]);

        delete [] reqs;
        delete [] stats;

        delete [] mpi_recv_buffer;
        delete [] mpi_send_buffer;
        delete [] str;

        delete [] msgSendSize;
        delete [] msgRecvSize;
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
    std::string Info() const override
    {
        return "MetisDivideHeterogeneousInputInMemoryProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetisDivideHeterogeneousInputInMemoryProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPartIO& mrSerialIO;

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

    // Copy constructor.
    //MetisDivideHeterogeneousInputInMemoryProcess(MetisDivideHeterogeneousInputInMemoryProcess const& rOther);


    ///@}

}; // Class MetisDivideHeterogeneousInputInMemoryProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} // addtogroup block

}

#endif // KRATOS_METIS_DIVIDE_HETEROGENEOUS_INPUT_IN_MEMORY_PROCESS_H
