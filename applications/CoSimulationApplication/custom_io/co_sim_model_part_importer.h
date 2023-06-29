// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#pragma once

// System includes
#include <vector>
#include <tuple>
#include <set>
#include <algorithm>
#include <memory>

// Project includes
#include "containers/model.h"
#include "includes/data_communicator.h"

// Application includes
#include "custom_utilities/co_sim_io_conversion_utilities.h"
#include "custom_external_libraries/CoSimIO/co_sim_io/includes/model_part.hpp"

namespace Kratos
{

class CoSimModelPartImporter
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static void FillCoSimIOModelPart(
        ModelPart& rKratosModelPart,
        const std::vector<std::pair<int, bool>>& rNodeIdAndPartitionIdList,
        const std::vector<std::vector<double>>& rNodeCoordinates,
        const DataCommunicator& rDataCommunicator)
    {
        KRATOS_TRY

        CoSimIO::ModelPart co_sim_model_part("dummy");

        KRATOS_ERROR_IF_NOT(rNodeIdAndPartitionIdList.size() == rNodeCoordinates.size())
            << "rNodeIdAndPartitionIdList and rNodeCoordinates size mismatch. [ rNodeIdAndPartitionIdList.size() = "
            << rNodeIdAndPartitionIdList.size() << ", rNodeCoordinates.size() = " << rNodeCoordinates.size() << " ]\n";


        // first do the communication and find out corresponding rank of each node
        std::vector<int> local_global_ids_list;
        local_global_ids_list.reserve(rNodeIdAndPartitionIdList.size());
        for (IndexType i = 0; i < rNodeIdAndPartitionIdList.size(); ++i) {
            const auto& id_parition_pair = rNodeIdAndPartitionIdList[i];

            if (std::get<1>(id_parition_pair)) {
                local_global_ids_list.push_back(std::get<0>(id_parition_pair));
            }
        }

        const auto& r_all_global_node_ids = rDataCommunicator.AllGatherv(local_global_ids_list);
        std::vector<std::set<int>> r_all_global_node_ids_set(r_all_global_node_ids.size());
        std::transform(
            r_all_global_node_ids.begin(),
            r_all_global_node_ids.end(),
            r_all_global_node_ids_set.begin(),
            [](auto& rVector) { return std::set<int>(rVector.begin(), rVector.end()); });

        for (IndexType i = 0; i < rNodeIdAndPartitionIdList.size(); ++i) {
            const auto& id_parition_pair = rNodeIdAndPartitionIdList[i];
            const auto& node_locations = rNodeCoordinates[i];

            const int global_id = std::get<0>(id_parition_pair);
            const bool is_local = std::get<1>(id_parition_pair);

            if (is_local) {
                co_sim_model_part.CreateNewNode(global_id, node_locations[0], node_locations[1], node_locations[2]);
            } else {
                // now get the rank
                IndexType rank;
                for (rank = 0; rank < r_all_global_node_ids.size(); ++rank) {
                    const auto itr = r_all_global_node_ids_set[rank].find(global_id);
                    if (itr != r_all_global_node_ids_set[rank].end()) {
                        break;
                    }
                }

                KRATOS_ERROR_IF(rank == r_all_global_node_ids.size())
                    << "Global id = " << global_id << " not found in any ranks.\n";

                co_sim_model_part.CreateNewGhostNode(global_id, node_locations[0], node_locations[1], node_locations[2], rank);
            }
        }

        CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_model_part, rKratosModelPart, rDataCommunicator);

        KRATOS_CATCH("");
    }

    ///@}
};
}


