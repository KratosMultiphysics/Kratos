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
    using IndexType = std::size_t;
    using SizeType = std::size_t;
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/


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



    void InterfaceSearchStructureMPI::PrepareSearchIteration(const Kratos::Flags& rOptions,
                                                       const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                                       InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    {
        KRATOS_INFO("InterfaceSearchStructureMPI") << "BEFORE exchanging informations in MPI" << std::endl;
        // Apply tolerance to bounding boxes
        std::vector<double> bounding_boxes_with_tol;
        MapperUtilities::ComputeBoundingBoxesWithTolerance(mGlobalBoundingBoxes,
                                                           mSearchRadius*1.2, // apply +20%
                                                           bounding_boxes_with_tol);

        // set up the buffers
        int comm_rank;
        int comm_size;

        MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

        std::vector<int> send_sizes(comm_size, 0);
        std::vector<int> recv_sizes(comm_size, 0);

        BufferType send_buffer(comm_size); // TODO reserve space!
        BufferType recv_buffer(comm_size); // TODO reserve space!

        // Compute Candidate Partitions and fill the send buffer
        FillSendBuffer(send_buffer, send_sizes, bounding_boxes_with_tol);

        // Exchange the buffer sizes
        MPI_Alltoall(send_sizes.data(), 1, MPI_INT, recv_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);


        MPI_Barrier(MPI_COMM_WORLD); // TODO remove me!
        std::stringstream ss_send;
        ss_send << "Rank: " << comm_rank << "; ";

        for (const auto& i : send_sizes)
            ss_send << " , " << i;

        KRATOS_INFO("send_sizes") << ss_send.str() << std::endl;

        std::stringstream ss_recv;
        ss_recv << "Rank: " << comm_rank << "; ";

        for (const auto& i : recv_sizes)
            ss_recv << " , " << i;

        KRATOS_INFO("recv_sizes") << ss_recv.str() << std::endl;
        std::cout << std::endl;
        MPI_Barrier(MPI_COMM_WORLD); // TODO remove me!


        // Send Information to Candidate Partitions

        int num_comm_events     = 0;
        int num_comm_events_idx = 0;

        for(int i=0; i<comm_size; ++i)
        {
            if(i != comm_rank && recv_sizes[i]) num_comm_events++;
            if(i != comm_rank && send_sizes[i]) num_comm_events++;
        }

        std::vector<MPI_Request> reqs(num_comm_events);
        std::vector<MPI_Status> stats(num_comm_events);

        KRATOS_INFO("num_comm_events") << "Rank: " << comm_rank << "; " << num_comm_events << std::endl;

        //Set up all receive and send events
        for(int i=0; i<comm_size; ++i)
        {
            if(i != comm_rank && recv_sizes[i])
            {
                recv_buffer[i].resize(recv_sizes[i]);
                MPI_Irecv(recv_buffer[i].data(),recv_sizes[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[num_comm_events_idx++]);
            }

            if(i != comm_rank && send_sizes[i])
            {
                MPI_Isend(send_buffer[i].data(),send_sizes[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,&reqs[num_comm_events_idx++]);
            }
        }


        // KRATOS_INFO("DONE With Communication") << "Rank: " << comm_rank << std::endl;

        //wait untill all communications finish
        int err = MPI_Waitall(num_comm_events, reqs.data(), stats.data());

        KRATOS_ERROR_IF_NOT(err == MPI_SUCCESS) << "Error in exchanging the information for "
            << "the construction of the MapperInterfaceInfos in MPI" << std::endl;


        KRATOS_INFO("InterfaceSearchStructureMPI") << "AFTER exchanging informations in MPI" << std::endl;

        // Construct MapperInterfaceInfos
        mpMapperInterfaceInfosContainer->clear();
        mpMapperInterfaceInfosContainer->reserve(comm_size);

        // IndexType local_sys_idx = 0;
        // for (const auto& r_local_sys : (*mpMapperLocalSystems))
        // {
        //     if (!r_local_sys->HasInterfaceInfo()) // Only the local_systems that have not received an InterfaceInfo create a new one
        //     {
        //         const auto& r_coords = r_local_sys->Coordinates();
        //         (*mpMapperInterfaceInfosContainer).push_back(rpRefInterfaceInfo->Create(r_coords, local_sys_idx));
        //     }
        //     ++local_sys_idx;
        // }


        // TODO print info saying that ORIGIN_ONLY has no effect in MPI, the destination also has to be updated
        // => since the origin changes also the destination might be sent to other partitions!

        /*
        1. Check which partitions have part of Interface (or send BBoxes directly...? => Don't think this is a good solution)
        2. Compute Graph for DataExchange (or do async => do I need a graph in such a case?)
        3. Send BBoxes to partitions that have part of Interface
        4. Compute CandidatePartitions
        5. Send stuff to CandidatePartitions
        6. Create Objects in the CandidatePartitions
        Afterwards do local search ...
        */
        MPI_Barrier(MPI_COMM_WORLD);
        KRATOS_INFO("InterfaceSearchStructureMPI") << "Leaving PrepareSearchIteration" << std::endl;
    }

    void InterfaceSearchStructureMPI::FinalizeSearchIteration()
    {
        /*
        1. Check with which Partitions I have to communicate
        2. Exchange this info
        3. Compute CommunicationGraph
        4. Searialize InterfaceInfos and exchange them => Big question: how to handle the fact that the interfaceinfos are pointer?
        I don't want to send pointer, so I have to think abt how to do it in a smart way...
        5. Deserialize InterfaceInfos and assign them to the LocalSystems
        */

    }

    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/
    void InterfaceSearchStructureMPI::FillSendBuffer(BufferType& rSendBuffer,
                                                     std::vector<int>& rSendSizes,
                                                     const std::vector<double>& rBoundingBoxes)
    {
        const SizeType comm_size = rSendBuffer.size();
        const SizeType num_local_sys = mpMapperLocalSystems->size();

        std::vector<double> bounding_box(6);

        // #pragma omp parallel for => "bounding_box" has to be threadprivate!
        for (IndexType i_rank=0; i_rank<comm_size; ++i_rank)
        {
            for (IndexType i=0; i<num_local_sys; ++i)
            {
                for (IndexType j=0; j<6; ++j)
                    bounding_box[j] = rBoundingBoxes[(i_rank*6) + j]; // retrieve bounding box of partition

                const auto& rp_local_sys = (*mpMapperLocalSystems)[i];

                if (!rp_local_sys->HasInterfaceInfo())
                {
                    const auto& r_coords = rp_local_sys->Coordinates();
                    if (MapperUtilities::PointIsInsideBoundingBox(bounding_box, r_coords))
                    {
                        // These push_backs are threadsafe!
                        rSendBuffer[i_rank].push_back(static_cast<double>(i)); // this it the "mSourceLocalSystemIndex" of the MapperInterfaceInfo
                        rSendBuffer[i_rank].push_back(r_coords[0]);
                        rSendBuffer[i_rank].push_back(r_coords[1]);
                        rSendBuffer[i_rank].push_back(r_coords[2]);

                        rSendSizes[i_rank] += 4;
                    }
                }
            }
        }
    }

}  // namespace Kratos.
