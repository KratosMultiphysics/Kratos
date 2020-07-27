// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//

// System includes

// External includes

// Project includes
#include "omp_node_search.h"

namespace Kratos {

void OMP_NodeSearch::SearchNodesInRadius(
        NodesContainerType const& rStructureNodes,
        int const Id,
        double const Radius,
        ResultNodesContainerType& rResults ) {

    KRATOS_TRY
    int max_num_of_nodes = rStructureNodes.size();

    NodesContainerType::ContainerType& nodes_array = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());

    ResultNodesContainerType local_results(max_num_of_nodes);
    std::size_t num_of_results = 0;

    ResultNodesContainerType::iterator ResultsPointer = local_results.begin();

    num_of_results = mBins->SearchObjectsInRadius( nodes_array[Id],Radius,ResultsPointer,max_num_of_nodes);

    rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

    KRATOS_CATCH("")
}

} // namespace Kratos