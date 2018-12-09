//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "mpi/includes/mpi_coloring_utilities.h"

namespace Kratos
{

std::vector<int> MPIColoringUtilities::ComputeRecvList(
    const std::vector<int>& local_destination_ids,
    MPIDataCommunicator& comm
)
{
    const int global_rank = 0;
    int current_rank = comm.Rank();

    //compute recv_list
    std::vector<std::vector<int>> recv_list;

    //gather everything on processor 0
   auto global_destination_ids = comm.Gatherv(local_destination_ids, global_rank );

    if(current_rank == global_rank)
    {
        recv_list.resize(global_destination_ids.size());

        for(unsigned int i=0; i<global_destination_ids.size(); ++i)
            for(const int& dest : global_destination_ids[i])
                recv_list[dest].push_back(i);
    }
    auto local_recv = comm.Scatterv(recv_list, global_rank );

    std::sort(local_recv.begin(), local_recv.end());
    std::unique(local_recv.begin(), local_recv.end());

    return local_recv;
}





}  // namespace Kratos.


