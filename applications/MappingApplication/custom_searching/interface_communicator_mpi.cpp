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

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "interface_communicator_mpi.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{

typedef std::size_t SizeType;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
InterfaceCommunicatorMPI::InterfaceCommunicatorMPI(ModelPart& rModelPartOrigin,
                                MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                                Parameters SearchSettings) :
    InterfaceCommunicator(rModelPartOrigin,
                                 pMapperLocalSystems,
                                 SearchSettings)
{
    // set up the buffers
    MPI_Comm_rank(MPI_COMM_WORLD, &mCommRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mCommSize);

    mpMapperInterfaceInfosContainer->resize(mCommSize);

    mSendSizes.resize(mCommSize);
    mRecvSizes.resize(mCommSize);

    mSendBufferDouble.resize(mCommSize);
    mRecvBufferDouble.resize(mCommSize);

    mSendBufferChar.resize(mCommSize);
    mRecvBufferChar.resize(mCommSize);
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceCommunicatorMPI::InitializeSearch(const Kratos::Flags& rOptions,
                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
{
    InterfaceCommunicator::InitializeSearch(rOptions, rpRefInterfaceInfo, InterfaceObjectTypeOrigin);

    // Exchange Bounding Boxes => has to be done every time
    ComputeGlobalBoundingBoxes();
}

void InterfaceCommunicatorMPI::InitializeSearchIteration(const Kratos::Flags& rOptions,
                                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    // Apply tolerance to bounding boxes
    std::vector<double> bounding_boxes_with_tol;
    MapperUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
                                                       mSearchRadius,
                                                       bounding_boxes_with_tol);

    // Compute Candidate Partitions and fill the send buffer
    MapperUtilities::FillBufferBeforeLocalSearch(mpMapperLocalSystems,
                                                 bounding_boxes_with_tol,
                                                 GetBufferSizeEstimate(),
                                                 mSendBufferDouble,
                                                 mSendSizes);

    // copy the local information directly
    mRecvBufferDouble[mCommRank] = mSendBufferDouble[mCommRank];

    const int err = ExchangeDataAsync(mSendBufferDouble, mRecvBufferDouble);

    KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
        << "the construction of the MapperInterfaceInfos in MPI" << std::endl;

    // Construct MapperInterfaceInfos
    MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mRecvBufferDouble,
                                                          rpRefInterfaceInfo,
                                                          mCommRank,
                                                          mpMapperInterfaceInfosContainer);

    // TODO print info saying that ORIGIN_ONLY has no effect in MPI, the destination also has to be updated
    // => since the origin changes also the destination might be sent to other partitions!

    MPI_Barrier(MPI_COMM_WORLD);
}

void InterfaceCommunicatorMPI::FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    // Reset to zero
    std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
    std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

    FilterInterfaceInfosSuccessfulSearch();

    MapperUtilities::FillBufferAfterLocalSearch(mpMapperInterfaceInfosContainer,
                                                rpRefInterfaceInfo,
                                                mCommRank,
                                                mSendBufferChar,
                                                mSendSizes);

    const int err = ExchangeDataAsync(mSendBufferChar, mRecvBufferChar);

    KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the "
        << "serialized MapperInterfaceInfos in MPI" << std::endl;

    MapperUtilities::DeserializeMapperInterfaceInfosFromBuffer(mRecvBufferChar,
                                                               rpRefInterfaceInfo,
                                                               mCommRank,
                                                               mpMapperInterfaceInfosContainer);

    AssignInterfaceInfos();

    MPI_Barrier(MPI_COMM_WORLD);
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/
void InterfaceCommunicatorMPI::ComputeGlobalBoundingBoxes()
{
    const auto local_bounding_box = MapperUtilities::ComputeLocalBoundingBox(mrModelPartOrigin); // TODO is const ok here?

    if (static_cast<int>(mGlobalBoundingBoxes.size()) != 6*mCommSize) {
        mGlobalBoundingBoxes.resize(6*mCommSize);
    }

    MPI_Allgather(local_bounding_box.data(),   6, MPI_DOUBLE,
                  mGlobalBoundingBoxes.data(), 6, MPI_DOUBLE,
                  MPI_COMM_WORLD); // TODO working?
}

inline MPI_Datatype GetMPIDatatype(const double& rValue) { return MPI_DOUBLE; }
inline MPI_Datatype GetMPIDatatype(const char& rValue)   { return MPI_CHAR; }

template< typename TDataType >
int InterfaceCommunicatorMPI::ExchangeDataAsync(
    const std::vector<std::vector<TDataType>>& rSendBuffer,
    std::vector<std::vector<TDataType>>& rRecvBuffer)
{
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
template int InterfaceCommunicatorMPI::ExchangeDataAsync<double>(
    const std::vector<std::vector<double>>& rSendBuffer,
    std::vector<std::vector<double>>& rRecvBuffer);

template int InterfaceCommunicatorMPI::ExchangeDataAsync<char>(
    const std::vector<std::vector<char>>& rSendBuffer,
    std::vector<std::vector<char>>& rRecvBuffer);

}  // namespace Kratos.
