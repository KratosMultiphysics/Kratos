//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Vicente Mataix Ferrandiz

// System includes
#include <cmath>

// External includes
#include "mpi.h"

// Project includes
#include "includes/exception.h"
#include "mpi/utilities/mpi_search_utilities.h"

namespace Kratos
{

void MPISearchData::Initialize()
{
    // Set up the buffers
    MPI_Comm_rank(MPI_COMM_WORLD, &CommRank);
    MPI_Comm_size(MPI_COMM_WORLD, &CommSize);

    SendSizes.resize(CommSize);
    RecvSizes.resize(CommSize);

    SendBufferDouble.resize(CommSize);
    RecvBufferDouble.resize(CommSize);

    SendBufferChar.resize(CommSize);
    RecvBufferChar.resize(CommSize);
}

/***********************************************************************************/
/***********************************************************************************/

void MPISearchUtilities::ComputeBoundingBoxesWithTolerance(
    const std::vector<double>& rBoundingBoxes,
    const double Tolerance,
    std::vector<double>& rBoundingBoxesWithTolerance
    )
{
    const SizeType size_vec = rBoundingBoxes.size();

    KRATOS_DEBUG_ERROR_IF_NOT(std::fmod(size_vec, 6) == 0) << "Bounding Boxes size has to be a multiple of 6!" << std::endl;

    if (rBoundingBoxesWithTolerance.size() != size_vec) {
        rBoundingBoxesWithTolerance.resize(size_vec);
    }

    // Apply Tolerances
    for (IndexType i=0; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] + Tolerance;
    }

    for (IndexType i=1; i<size_vec; i+=2) {
        rBoundingBoxesWithTolerance[i] = rBoundingBoxes[i] - Tolerance;
    }
}

/***********************************************************************************/
/***********************************************************************************/

inline MPI_Datatype GetMPIDatatype(const double& rValue) { return MPI_DOUBLE; }
inline MPI_Datatype GetMPIDatatype(const char& rValue)   { return MPI_CHAR; }

template<typename TDataType>
int MPISearchUtilities::ExchangeDataAsync(
    const std::vector<std::vector<TDataType>>& rSendBuffer,
    std::vector<std::vector<TDataType>>& rRecvBuffer,
    MPISearchData& rSearchData
    )
{
    // Exchange the buffer sizes
    MPI_Alltoall(rSearchData.SendSizes.data(), 1, MPI_INT, rSearchData.RecvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Send Information to Candidate Partitions
    int num_comm_events     = 0;
    int num_comm_events_idx = 0;

    for(int i=0; i<rSearchData.CommSize; ++i) {
        if(i != rSearchData.CommRank && rSearchData.RecvSizes[i]) num_comm_events++;
        if(i != rSearchData.CommRank && rSearchData.SendSizes[i]) num_comm_events++;
    }

    // TODO make members?
    std::vector<MPI_Request> reqs(num_comm_events);
    std::vector<MPI_Status> stats(num_comm_events);

    const MPI_Datatype mpi_datatype(GetMPIDatatype(TDataType()));

    // Exchange the data
    for (int i=0; i<rSearchData.CommSize; ++i) {
        if (i != rSearchData.CommRank && rSearchData.RecvSizes[i]) { // TODO check what "rSearchData.RecvSizes[i]" returns
            if (rRecvBuffer[i].size() != static_cast<SizeType>(rSearchData.RecvSizes[i])) {
                rRecvBuffer[i].resize(rSearchData.RecvSizes[i]);
            }

            MPI_Irecv(rRecvBuffer[i].data(), rSearchData.RecvSizes[i],
                      mpi_datatype, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }

        if (i != rSearchData.CommRank && rSearchData.SendSizes[i]) {
            MPI_Isend(rSendBuffer[i].data(), rSearchData.SendSizes[i],
                      mpi_datatype, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }
    }

    //wait until all communications finish
    const int err = MPI_Waitall(num_comm_events, reqs.data(), stats.data());

    return err;
}

// Explicit instantiation for function
template int MPISearchUtilities::ExchangeDataAsync<double>(
    const std::vector<std::vector<double>>& rSendBuffer,
    std::vector<std::vector<double>>& rRecvBuffer,
    MPISearchData& rSearchData
    );

template int MPISearchUtilities::ExchangeDataAsync<char>(
    const std::vector<std::vector<char>>& rSendBuffer,
    std::vector<std::vector<char>>& rRecvBuffer,
    MPISearchData& rSearchData
    );

} // namespace Kratos