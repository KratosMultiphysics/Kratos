//     __ __      __  __________  ___                ___            __  _           
//    / //_/___ _/ / / /  _/ __ \/   |  ____  ____  / (_)________ _/ /_(_)___  ____ 
//   / ,< / __ `/ /_/ // // /_/ / /| | / __ \/ __ \/ / / ___/ __ `/ __/ / __ \/ __
//  / /| / /_/ / __  // // ____/ ___ |/ /_/ / /_/ / / / /__/ /_/ / /_/ / / /_/ / / /
// /_/ |_\__,_/_/ /_/___/_/   /_/  |_/ .___/ .___/_/_/\___/\__,_/\__/_/\____/_/ /_/ 
//                                  /_/   /_/                                       
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "custom_utilities/kahip_csr_converter.h"
#include "includes/exception.h"

namespace Kratos 
{

void KaHIPCSRConverter::ConvertToCSRFormat(
    const ConnectivitiesContainerType& rKratosFormatConnectivities,
    std::vector<kahip_idx>& rXAdj,
    std::vector<kahip_idx>& rAdjncy)
{
    const SizeType num_nodes = rKratosFormatConnectivities.size();

    // Count total entries to pre-allocate adjncy
    SizeType total_entries = 0;
    for (const auto& neighbors : rKratosFormatConnectivities) {
        total_entries += neighbors.size();
    }

    rXAdj.resize(num_nodes + 1);
    rAdjncy.resize(total_entries);

    rXAdj[0] = 0;
    SizeType adj_pos = 0;

    for (SizeType i = 0; i < num_nodes; ++i) {
        for (SizeType neighbor_id : rKratosFormatConnectivities[i]) {
            // Kratos node IDs are 1-indexed; KaHIP CSR uses 0-indexed
            KRATOS_DEBUG_ERROR_IF(neighbor_id == 0)
                << "KaHIPCSRConverter: encountered a zero node ID in connectivity (IDs must start at 1)." << std::endl;
            rAdjncy[adj_pos++] = static_cast<kahip_idx>(neighbor_id - 1);
        }
        rXAdj[i + 1] = static_cast<kahip_idx>(adj_pos);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void KaHIPCSRConverter::ConvertToCSRFormatWithWeights(
    const ConnectivitiesContainerType& rKratosFormatConnectivities,
    std::vector<kahip_idx>& rXAdj,
    std::vector<kahip_idx>& rAdjncy,
    std::vector<int>& rVWgt,
    std::vector<kahip_idx>& rAdjcWgt)
{
    ConvertToCSRFormat(rKratosFormatConnectivities, rXAdj, rAdjncy);

    const SizeType num_nodes = rKratosFormatConnectivities.size();
    const SizeType num_edges = rAdjncy.size();

    // Uniform vertex and edge weights
    rVWgt.assign(num_nodes, 1);
    rAdjcWgt.assign(num_edges, 1);
}

/***********************************************************************************/
/***********************************************************************************/

bool KaHIPCSRConverter::ValidateCSRGraph(
    int n,
    const std::vector<kahip_idx>& rXAdj,
    const std::vector<kahip_idx>& rAdjncy)
{
    bool valid = true;

    if (static_cast<int>(rXAdj.size()) != n + 1) {
        KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: xadj size "
            << rXAdj.size() << " != n+1 = " << (n + 1) << std::endl;
        valid = false;
    }

    if (rXAdj[0] != 0) {
        KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: xadj[0] = "
            << rXAdj[0] << " (expected 0)" << std::endl;
        valid = false;
    }

    // Check monotone and bounds for adjncy
    for (int i = 0; i < n; ++i) {
        if (rXAdj[i] > rXAdj[i + 1]) {
            KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: xadj not monotone at i="
                << i << " (xadj[i]=" << rXAdj[i] << " > xadj[i+1]=" << rXAdj[i + 1] << ")" << std::endl;
            valid = false;
        }

        for (kahip_idx j = rXAdj[i]; j < rXAdj[i + 1]; ++j) {
            const kahip_idx neighbor = rAdjncy[j];

            if (neighbor < 0 || neighbor >= static_cast<kahip_idx>(n)) {
                KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: adjncy["
                    << j << "] = " << neighbor << " out of bounds [0, " << n << ")" << std::endl;
                valid = false;
            }

            if (neighbor == static_cast<kahip_idx>(i)) {
                KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: self-loop at node "
                    << i << std::endl;
                valid = false;
            }
        }
    }

    // Check symmetry: for every edge (i->j) there must be an edge (j->i)
    for (int i = 0; i < n && valid; ++i) {
        for (kahip_idx j_idx = rXAdj[i]; j_idx < rXAdj[i + 1]; ++j_idx) {
            const int j = static_cast<int>(rAdjncy[j_idx]);

            bool found_reverse = false;
            for (kahip_idx k_idx = rXAdj[j]; k_idx < rXAdj[j + 1]; ++k_idx) {
                if (rAdjncy[k_idx] == static_cast<kahip_idx>(i)) {
                    found_reverse = true;
                    break;
                }
            }

            if (!found_reverse) {
                KRATOS_WARNING("KaHIPCSRConverter") << "ValidateCSRGraph: edge ("
                    << i << "->" << j << ") has no reverse edge — graph is not symmetric." << std::endl;
                valid = false;
                break;
            }
        }
    }

    return valid;
}

} // namespace Kratos
