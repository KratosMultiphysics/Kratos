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

#if !defined(KRATOS_MAPPER_UTILITIES_MPI_H_INCLUDED )
#define  KRATOS_MAPPER_UTILITIES_MPI_H_INCLUDED

// System includes
#include <iomanip> // for "std::setw"

// External includes
#include "mpi.h"

// Project includes
#include "includes/define.h"
#include "mapper_utilities.h"
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

/// Auxiliary functions for the MappingApplication in MPI
/** This class provides a set of auxiliary functions for MPI that are used by several other functions / classes
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
    static inline int SizeOfVariable(const T& rVariable);

    template<typename T>
    static void MpiSendRecv(T* pSendBuffer, T* pReceiveBuffer, const int SendBufferSize,
                            int& ReceiveBufferSize, const int CommPartner)
    {
        // Determine data type
        MPI_Datatype DataType = MapperUtilitiesMPI::GetMPIDatatype(T());

        //   Exchange the information about the receiving buffer size
        MPI_Sendrecv(&SendBufferSize, 1, MPI_INT, CommPartner, 0, &ReceiveBufferSize,
                     1, MPI_INT, CommPartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Perform the actual exchange of the buffers
        MPI_Sendrecv(pSendBuffer, SendBufferSize, DataType, CommPartner, 0, pReceiveBuffer,
                     ReceiveBufferSize, DataType, CommPartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // This function checks if the buffer sizes are large enough
    template<typename T>
    static void MpiSendRecv(T* pSendBuffer, T* pReceiveBuffer, const int SendBufferSize,
                            int& rReceiveBufferSize, const int MaxSendBufferSize,
                            const int MaxReceiveBufferSize, const int CommPartner)
    {
        // Determine data type
        MPI_Datatype DataType = MapperUtilitiesMPI::GetMPIDatatype(T());

        //   Exchange the information about the receiving buffer size
        MPI_Sendrecv(&SendBufferSize, 1, MPI_INT, CommPartner, 0, &rReceiveBufferSize,
                     1, MPI_INT, CommPartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        KRATOS_DEBUG_ERROR_IF(SendBufferSize > MaxSendBufferSize)
            << "Send Buffer too small!; "
            << "SendBufferSize = " << SendBufferSize
            << ", MaxSendBufferSize = "
            << MaxSendBufferSize << std::endl;

        KRATOS_DEBUG_ERROR_IF(rReceiveBufferSize > MaxReceiveBufferSize)
            << "Receive Buffer too small!; "
            << "rReceiveBufferSize = " << rReceiveBufferSize
            << ", MaxReceiveBufferSize = "
            << MaxReceiveBufferSize << std::endl;

        // Perform the actual exchange of the buffers
        MPI_Sendrecv(pSendBuffer, SendBufferSize, DataType, CommPartner, 0, pReceiveBuffer,
                     rReceiveBufferSize, DataType, CommPartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    static void ComputeMaxBufferSizes(int* pLocalMemorySizeArray,
                                      int& rSendBufferSize,
                                      int& rReceiveBufferSize,
                                      const int CommRank,
                                      const int CommSize)
    {
        int* global_memory_size_array = new int[CommSize]();

        // Go through list of how many objects will be sent and find maximum
        rSendBufferSize = pLocalMemorySizeArray[0];
        for (int i = 1; i < CommSize; ++i)
        {
            rSendBufferSize = std::max(pLocalMemorySizeArray[i], rSendBufferSize);
        }

        // TODO comment that it could be done more efficient, but only with a lot of code => what Charlie explained
        // Basically splitting the matrix and doing partial reduces
        // Exchange information about how much is going to be sent around among the
        // partitions and take the maximum that my partition will receive
        MPI_Allreduce(pLocalMemorySizeArray, global_memory_size_array,
                      CommSize, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        rReceiveBufferSize = global_memory_size_array[CommRank];

        delete [] global_memory_size_array;
    }

    static void ComputeColoringGraph(int* pLocalCommList,
                                     int NumPartitions,
                                     GraphType& rDomainsColoredGraph,
                                     int& rMaxColors)
    {
        int comm_list_size = NumPartitions * NumPartitions;
        int* global_comm_list = new int[comm_list_size]();

        MPI_Allgather(pLocalCommList, NumPartitions, MPI_INT, global_comm_list,
                      NumPartitions, MPI_INT, MPI_COMM_WORLD);

        GraphType domains_graph;
        domains_graph.resize(NumPartitions, NumPartitions, false);
        domains_graph = ScalarMatrix(NumPartitions, NumPartitions, -1.0f);

        // set up communication matrix as input for the coloring process
        for (int i = 0; i < NumPartitions; ++i)
        {
            for (int j = 0; j < NumPartitions; ++j)
            {
                domains_graph(i, j) = global_comm_list[(i * NumPartitions) + j];
            }
        }

        // make matrix symmetric, using the larger value (i.e. 1, meaning that communication is required)
        // needed if for some ranks only one-sided communication is required
        for(std::size_t i = 0; i < domains_graph.size1(); ++i)
        {
            for(std::size_t j = i + 1; j < domains_graph.size2(); ++j)
            {
                int val = std::max(domains_graph(i, j), domains_graph(j, i));
                domains_graph(i, j) = val;
                domains_graph(j, i) = val;
            }
            domains_graph(i, i) = 0; // set main diagonal to zero, i.e. avoid communication with myself
        }

        GraphColoringProcess(NumPartitions, domains_graph,
                             rDomainsColoredGraph, rMaxColors).Execute();

        delete [] global_comm_list;
    }

    static void PrintGraph(GraphType& rGraph, int NumColors)
    {
        std::cout << "Mapper Communication Graph: " << std::endl;
        std::cout << std::setw(5);
        for(std::size_t i = 0; i < rGraph.size1(); ++i)
        {
            for(int j = 0; j < NumColors; ++j)
            {
                std::cout << rGraph(i, j) << std::setw(5);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    static void PrintRankLocalBoundingBox(double* pBoundingBox, const int CommRank)
    {
        std::cout << "Rank-local Bounding Box, Rank " << CommRank << ": "
                  << MapperUtilities::BoundingBoxStingStream(pBoundingBox)
                  << std::endl;
        // TODO maybe write it such that gid can directly read it
    }

    static void ComputeGlobalBoundingBoxes(double* pLocalBoundingBox,
                                           const double Tolerance,
                                           const int CommRank,
                                           const int EchoLevel,
                                           double* pGlobalBoundingBoxes)
    {
        double* local_bounding_box_tol = new double[6];
        MapperUtilities::ComputeBoundingBoxWithTolerance(pLocalBoundingBox,
                Tolerance,
                local_bounding_box_tol);
        
        if (EchoLevel >= 3)
        {
            MapperUtilitiesMPI::PrintRankLocalBoundingBox(local_bounding_box_tol, CommRank);
            // TODO maybe write it such that gid can directly read it
        }

        MPI_Allgather(local_bounding_box_tol, 6, MPI_DOUBLE,
                      pGlobalBoundingBoxes, 6, MPI_DOUBLE, MPI_COMM_WORLD);
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
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperUtilitiesMPI";
    }

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
    static inline MPI_Datatype GetMPIDatatype(const T& rValue);


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
inline MPI_Datatype MapperUtilitiesMPI::GetMPIDatatype<int>(const int& rValue)
{
    return MPI_INT ;
}

template<>
inline MPI_Datatype MapperUtilitiesMPI::GetMPIDatatype<double>(const double& rValue)
{
    return MPI_DOUBLE ;
}

template<>
inline int MapperUtilitiesMPI::SizeOfVariable< double >(const double& rVariable)
{
    return 1;
}

template<>
inline int MapperUtilitiesMPI::SizeOfVariable< array_1d<double, 3> >(const array_1d<double, 3>& rVariable)
{
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