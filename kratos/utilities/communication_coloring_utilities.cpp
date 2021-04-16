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
#include "utilities/communication_coloring_utilities.h"

namespace Kratos
{

std::vector<int> MPIColoringUtilities::ComputeRecvList(
    const std::vector<int>& rLocalDestinationIds,
    const DataCommunicator& rComm
)
{
    const int global_rank = 0;
    int current_rank = rComm.Rank();

    //compute recv_list
    std::vector<std::vector<int>> recv_list;

    //gather everything on processor 0
    auto global_destination_ids = rComm.Gatherv(rLocalDestinationIds, global_rank );

    if(current_rank == global_rank)
    {
        recv_list.resize(global_destination_ids.size());

        for(unsigned int i=0; i<global_destination_ids.size(); ++i)
            for(const int& dest : global_destination_ids[i])
                recv_list[dest].push_back(i);
    }
    auto local_recv = rComm.Scatterv(recv_list, global_rank );

    auto last = std::unique(local_recv.begin(), local_recv.end());
    local_recv.erase(last, local_recv.end());

    return local_recv;
}

std::vector<int> MPIColoringUtilities::ComputeCommunicationScheduling(
    const std::vector<int>& rLocalDestinationIds,
    const DataCommunicator& rComm
)
{
    const int global_rank = 0;
    int current_rank = rComm.Rank();

    auto global_destination_ids = rComm.Gatherv(rLocalDestinationIds, global_rank );

    std::vector< std::vector< int >> global_colors;
    if(current_rank == global_rank)
    {
        global_colors.resize(global_destination_ids.size());

        std::map<int, std::map<int, int> > input_graph;
        for(unsigned int i=0; i<global_destination_ids.size(); ++i)
        {
            for(const auto j : global_destination_ids[i] )
            {
                (input_graph[i])[j] = 1;
                (input_graph[j])[i] = 1;
            }
        }

        int max_colors = 0;
        std::map<int, std::map<int, int> > colored_graph;
        for(unsigned int i=0; i<input_graph.size(); ++i)
        {
            for(const auto& item : input_graph[i] )
            {
                unsigned int j = item.first;
                if(j > i)
                {
                    bool found = false;
                    int color = 0;
                    while(!found)
                    {
                        if( !HasEdge(colored_graph,i,color) && !HasEdge(colored_graph,j,color) )
                        {
                            colored_graph[i][color] = j;
                            colored_graph[j][color] = i;
                            if(max_colors < static_cast<int>(color + 1))
                                max_colors = color + 1;
                            found = true;
                        }
                        color += 1;
                    }
                }
            }
        }

        for(auto& item : global_colors)
        {
            item.resize(max_colors, -1);
        }

        for(unsigned int i=0; i<colored_graph.size();  ++i)
        {
            for(const auto& item : colored_graph[i] )
            {
                global_colors[i][item.first] = item.second;
            }

            // //output, useful for debug
            // std::cout << " proc : " << i << std::endl;
            // for(auto& item : global_colors[i])
            //     std::cout << item << " ";
            // std::cout << std::endl;
        }
    }

    auto colors = rComm.Scatterv(global_colors, global_rank );

    return colors;
}

bool MPIColoringUtilities::HasEdge(std::map<int, std::map<int, int> >& rGraph,
                                   int i,
                                   int j)
{
    auto it = rGraph.find(i);
    if( it != rGraph.end())
        if(it->second.find(j) != it->second.end())
            return true;
    return false;
}


}  // namespace Kratos.


