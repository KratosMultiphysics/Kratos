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
#include "node_search_utility.h"

namespace Kratos {

void NodeSearchUtility::SearchNodesInRadius(
        NodeType::Pointer pNode,
        double const Radius,
        ResultNodesContainerType& rResults ) {

    KRATOS_TRY

    ResultNodesContainerType local_results(mMaxNumberOfNodes);
    std::size_t num_of_results = 0;

    auto results_iterator = local_results.begin();

    num_of_results = mpBins->SearchObjectsInRadius(pNode, Radius, results_iterator, mMaxNumberOfNodes);

    rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

    KRATOS_CATCH("")
}

} // namespace Kratos