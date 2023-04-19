//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/geometrical_object.h"
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos
{
void GeometricalObjectsBinsMPI::ComputeGlobalBoundingBoxes()
{
    // TODO
}

/***********************************************************************************/
/***********************************************************************************/

inline MPI_Datatype GetMPIDatatype(const double& rValue) { return MPI_DOUBLE; }
inline MPI_Datatype GetMPIDatatype(const char& rValue)   { return MPI_CHAR; }

template< typename TDataType >
int GeometricalObjectsBinsMPI::ExchangeDataAsync(
    const std::vector<std::vector<TDataType>>& rSendBuffer,
    std::vector<std::vector<TDataType>>& rRecvBuffer
    )
{
    using SizeType = std::size_t;

    // Exchange the buffer sizes
    MPI_Alltoall(mSendSizes.data(), 1, MPI_INT, mRecvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Send Information to Candidate Partitions
    int num_comm_events     = 0;
    int num_comm_events_idx = 0;

    for(int i=0; i<mCommSize; ++i) {
        if(i != mCommRank && mRecvSizes[i]) num_comm_events++;
        if(i != mCommRank && mSendSizes[i]) num_comm_events++;
    }

    // TODO make members?
    std::vector<MPI_Request> reqs(num_comm_events);
    std::vector<MPI_Status> stats(num_comm_events);

    const MPI_Datatype mpi_datatype(GetMPIDatatype(TDataType()));

    // Exchange the data
    for (int i=0; i<mCommSize; ++i) {
        if (i != mCommRank && mRecvSizes[i]) { // TODO check what "mRecvSizes[i]" returns
            if (rRecvBuffer[i].size() != static_cast<SizeType>(mRecvSizes[i])) {
                rRecvBuffer[i].resize(mRecvSizes[i]);
            }

            MPI_Irecv(rRecvBuffer[i].data(), mRecvSizes[i],
                      mpi_datatype, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }

        if (i != mCommRank && mSendSizes[i]) {
            MPI_Isend(rSendBuffer[i].data(), mSendSizes[i],
                      mpi_datatype, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }
    }

    //wait until all communications finish
    const int err = MPI_Waitall(num_comm_events, reqs.data(), stats.data());

    return err;
}

// Explicit instantiation for function
template int GeometricalObjectsBinsMPI::ExchangeDataAsync<double>(
    const std::vector<std::vector<double>>& rSendBuffer,
    std::vector<std::vector<double>>& rRecvBuffer);

template int GeometricalObjectsBinsMPI::ExchangeDataAsync<char>(
    const std::vector<std::vector<char>>& rSendBuffer,
    std::vector<std::vector<char>>& rRecvBuffer);

}