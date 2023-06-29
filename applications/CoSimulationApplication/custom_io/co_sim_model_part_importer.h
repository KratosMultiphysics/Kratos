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

// Project includes
#include "includes/data_communicator.h"

// Application includes
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
        CoSimIO::ModelPart& rCoSimIOModelPart,
        const std::vector<std::pair<int, int>>& rNodeIdAndPartitionIdList,
        const std::vector<std::vector<double>>& rNodeCoordinates,
        const DataCommunicator& rDataCommunicator)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(rNodeIdAndPartitionIdList.size() == rNodeCoordinates.size()) 
            << "rNodeIdAndPartitionIdList and rNodeCoordinates size mismatch. [ rNodeIdAndPartitionIdList.size() = " 
            << rNodeIdAndPartitionIdList.size() << ", rNodeCoordinates.size() = " << rNodeCoordinates.size() << " ]\n";

        rDataCommunicator.

        for (IndexType i = 0; i < rNodeIdAndPartitionIdList.size(); ++i) {
            const auto& id_parition_pair = rNodeIdAndPartitionIdList[i];
            const auto& node_locations = rNodeCoordinates[i];

            const int node_id = std::get<0>(id_parition_pair);
            const int parition_id = std::get<1>(id_parition_pair);
            
            if (rDataCommunicator.Rank() == parition_id) {
                rCoSimIOModelPart.CreateNewNode(node_id, node_locations[0], node_locations[1], node_locations[2]);
            } else {
                rCoSimIOModelPart.CreateNewGhostNode(node_id, node_locations[0], node_locations[1], node_locations[2], 0);
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
};
}


