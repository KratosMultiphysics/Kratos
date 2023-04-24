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
#include "input_output/vtk_output.h"
#include "mpi/spatial_containers/geometrical_objects_bins_mpi.h"

namespace Kratos
{
void GeometricalObjectsBinsMPI::InitializeSearch()
{
    KRATOS_TRY;

    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    // // Apply tolerance to bounding boxes
    // std::vector<double> bounding_boxes_with_tol;
    // MapperUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
    //                                                    mSearchRadius,
    //                                                    bounding_boxes_with_tol);

    // // Compute Candidate Partitions and fill the send buffer
    // MapperUtilities::FillBufferBeforeLocalSearch(mrMapperLocalSystems,
    //                                              bounding_boxes_with_tol,
    //                                              GetBufferSizeEstimate(),
    //                                              mSendBufferDouble,
    //                                              mSendSizes);

    // // copy the local information directly
    // mRecvBufferDouble[mCommRank] = mSendBufferDouble[mCommRank];

    // const int err = ExchangeDataAsync(mSendBufferDouble, mRecvBufferDouble);

    // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
    //     << "the construction of the MapperInterfaceInfos in MPI" << std::endl;

    // // Construct MapperInterfaceInfos
    // MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mRecvBufferDouble,
    //                                                       rpRefInterfaceInfo,
    //                                                       mCommRank,
    //                                                       mMapperInterfaceInfosContainer);

    MPI_Barrier(MPI_COMM_WORLD);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::FinalizeSearch()
{
    KRATOS_TRY;

    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    // FilterInterfaceInfosSuccessfulSearch();

    // MapperUtilities::FillBufferAfterLocalSearch(mMapperInterfaceInfosContainer,
    //                                             rpRefInterfaceInfo,
    //                                             mCommRank,
    //                                             mSendBufferChar,
    //                                             mSendSizes);

    // const int err = ExchangeDataAsync(mSendBufferChar, mRecvBufferChar);

    // KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the "
    //     << "serialized MapperInterfaceInfos in MPI" << std::endl;

    // MapperUtilities::DeserializeMapperInterfaceInfosFromBuffer(mRecvBufferChar,
    //                                                            rpRefInterfaceInfo,
    //                                                            mCommRank,
    //                                                            mMapperInterfaceInfosContainer);

    // AssignInterfaceInfos();

    MPI_Barrier(MPI_COMM_WORLD);

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void GeometricalObjectsBinsMPI::ComputeGlobalBoundingBoxes()
{
    // xmax, xmin, ymax, ymin, zmax, zmin
    std::array<double, 6> local_bounding_box;
    auto& r_min_point = mBoundingBox.GetMinPoint();
    auto& r_max_point = mBoundingBox.GetMaxPoint();
    for (unsigned int i=0; i<3; ++i) {
        local_bounding_box[i*2]   = r_max_point[i];
        local_bounding_box[i*2+1] = r_min_point[i];
    }

    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*mCommSize) {
        mGlobalBoundingBoxes.resize(6*mCommSize);
    }

    MPI_Allgather(local_bounding_box.data(),   6, MPI_DOUBLE,
                  mGlobalBoundingBoxes.data(), 6, MPI_DOUBLE,
                  MPI_COMM_WORLD);
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