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
#include "interface_search_structure_mpi.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_utilities/mapper_utilities_mpi.h"

namespace Kratos
{
using SizeType = std::size_t;
using IndexType = std::size_t;

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
InterfaceSearchStructureMPI::InterfaceSearchStructureMPI(ModelPart& rModelPartOrigin,
                                MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                                Parameters SearchSettings) :
        InterfaceSearchStructureBase(rModelPartOrigin,
                                 pMapperLocalSystems,
                                 SearchSettings)
    {
        KRATOS_WATCH("MPI-Ctor of InterfaceSearchStructureMPI")

        // set up the buffers
        MPI_Comm_rank(MPI_COMM_WORLD, &mCommRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mCommSize);

        mSendSizes.resize(mCommSize);
        mRecvSizes.resize(mCommSize);

        std::fill(mSendSizes.begin(), mSendSizes.end(), 0);
        std::fill(mRecvSizes.begin(), mRecvSizes.end(), 0);

        mSendBufferDouble.resize(mCommSize);
        mRecvBufferDouble.resize(mCommSize);

        mSendBufferChar.resize(mCommSize); // TODO reserve space!
        mRecvBufferChar.resize(mCommSize);
    }

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/
void InterfaceSearchStructureMPI::PrepareSearch(const Kratos::Flags& rOptions,
                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                    InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
{
    InterfaceSearchStructureBase::PrepareSearch(rOptions, rpRefInterfaceInfo, InterfaceObjectTypeOrigin);

    // Exchange Bounding Boxes => has to be done every time
    MapperUtilitiesMPI::ComputeGlobalBoundingBoxes(mrModelPartOrigin, mGlobalBoundingBoxes);
}


void InterfaceSearchStructureMPI::InitializeSearchIteration(const Kratos::Flags& rOptions,
                                                    const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    KRATOS_INFO("InterfaceSearchStructureMPI") << "BEFORE exchanging informations in MPI" << std::endl;
    // Apply tolerance to bounding boxes
    std::vector<double> bounding_boxes_with_tol;
    MapperUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
                                                       mSearchRadius*1.2, // apply +20%
                                                       bounding_boxes_with_tol);

    // Compute Candidate Partitions and fill the send buffer
    MapperUtilities::FillBufferBeforeLocalSearch(mpMapperLocalSystems,
                                                 bounding_boxes_with_tol,
                                                 GetBufferSizeEstimate(),
                                                 mSendBufferDouble,
                                                 mSendSizes);

    // Exchange the buffer sizes
    MPI_Alltoall(mSendSizes.data(), 1, MPI_INT, mRecvSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);


    // MPI_Barrier(MPI_COMM_WORLD); // TODO remove me!
    // std::stringstream ss_send;
    // ss_send << "Rank: " << mCommRank << "; ";

    // for (const auto& i : mSendSizes)
    //     ss_send << " , " << i;

    // KRATOS_INFO("mSendSizes") << ss_send.str() << std::endl;

    // std::stringstream ss_recv;
    // ss_recv << "Rank: " << mCommRank << "; ";

    // for (const auto& i : mRecvSizes)
    //     ss_recv << " , " << i;

    // KRATOS_INFO("mRecvSizes") << ss_recv.str() << std::endl;
    // std::cout << std::endl;
    // MPI_Barrier(MPI_COMM_WORLD); // TODO remove me!


    // Send Information to Candidate Partitions

    int num_comm_events     = 0;
    int num_comm_events_idx = 0;

    for(int i=0; i<mCommSize; ++i)
    {
        if(i != mCommRank && mRecvSizes[i]) num_comm_events++;
        if(i != mCommRank && mSendSizes[i]) num_comm_events++;
    }

    std::vector<MPI_Request> reqs(num_comm_events);
    std::vector<MPI_Status> stats(num_comm_events);

    // KRATOS_INFO("num_comm_events") << "Rank: " << mCommRank << "; " << num_comm_events << std::endl;

    //Set up all receive and send events
    for (int i=0; i<mCommSize; ++i)
    {
        if(i != mCommRank && mRecvSizes[i]) // TODO check what "mRecvSizes[i]" returns
        {
            mRecvBufferDouble[i].resize(mRecvSizes[i]);
            MPI_Irecv(mRecvBufferDouble[i].data(), mRecvSizes[i],
                      MPI_DOUBLE, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }

        if(i != mCommRank && mSendSizes[i])
        {
            MPI_Isend(mSendBufferDouble[i].data(), mSendSizes[i],
                      MPI_DOUBLE, i, 0,
                      MPI_COMM_WORLD, &reqs[num_comm_events_idx++]);
        }
    }

    // KRATOS_INFO("DONE With Communication") << "Rank: " << mCommRank << std::endl;

    //wait untill all communications finish
    int err = MPI_Waitall(num_comm_events, reqs.data(), stats.data());

    KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
        << "the construction of the MapperInterfaceInfos in MPI" << std::endl;


    KRATOS_INFO("InterfaceSearchStructureMPI") << "AFTER exchanging informations in MPI" << std::endl;

    // Construct MapperInterfaceInfos
    MapperUtilities::CreateMapperInterfaceInfosFromBuffer(mRecvBufferDouble,
                                                          rpRefInterfaceInfo,
                                                          mCommRank,
                                                          mpMapperInterfaceInfosContainer);

    // TODO print info saying that ORIGIN_ONLY has no effect in MPI, the destination also has to be updated
    // => since the origin changes also the destination might be sent to other partitions!

    MPI_Barrier(MPI_COMM_WORLD);
    KRATOS_INFO("InterfaceSearchStructureMPI") << "Leaving InitializeSearchIteration" << std::endl;
}

void InterfaceSearchStructureMPI::FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo)
{
    MapperInterfaceInfoPointerVectorPointerType p_interface_infos_succ_search
        = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();

    MapperUtilities::SelectInterfaceInfosSuccessfulSearch(
        mpMapperInterfaceInfosContainer,
        GetBufferSizeEstimate(),
        p_interface_infos_succ_search);

    auto p_ref_interface_info = rpRefInterfaceInfo->Create(); // needed to "remove" the const

    MapperUtilities::MapperInterfaceInfoSerializer interface_infos_serializer_save(
        p_interface_infos_succ_search, p_ref_interface_info );

    Kratos::Serializer send_serializer;
    send_serializer.save("interface_infos", interface_infos_serializer_save);

    // Exchange Data in over MPI

    /* // TODO reenable this once the MPI-Data-Exchange is it will not work without!
    MapperUtilities::MapperInterfaceInfoSerializer interface_infos_serializer_load(
        mpMapperInterfaceInfosContainer, p_ref_interface_info );

    Kratos::Serializer load_serializer;
    // Get the stream into the serializer somehow ...
    load_serializer.load("interface_infos", interface_infos_serializer_load);

    MapperUtilities::AssignInterfaceInfosAfterRemoteSearch(
        mpMapperInterfaceInfosContainer,
        mpMapperLocalSystems);
    */
}

/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/


}  // namespace Kratos.
