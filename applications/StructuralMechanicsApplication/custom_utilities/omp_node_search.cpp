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

void OMP_NodeSearch::InitializeSearch(NodesContainerType const& rStructureNodes) {

    KRATOS_TRY;
    if(!mIsInitialized) {
        NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
        mBins = new NodeBinsType(nodes_ModelPart.begin(), nodes_ModelPart.end());

        mIsInitialized = true;
    }
    KRATOS_CATCH("");
}

void OMP_NodeSearch::SearchNodesInRadiusExclusiveImplementation (
        NodesContainerType const& rStructureNodes,
        int const Id,
        double const Radius,
        ResultNodesContainerType& rResults ) {

    KRATOS_TRY
    int MaxNumberOfNodes = rStructureNodes.size();

    NodesContainerType::ContainerType& nodes_array = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());

    ResultNodesContainerType localResults(MaxNumberOfNodes);
    std::size_t NumberOfResults = 0;

    ResultNodesContainerType::iterator ResultsPointer = localResults.begin();

    NumberOfResults = mBins->SearchObjectsInRadius( nodes_array[Id],Radius,ResultsPointer,MaxNumberOfNodes);

    rResults.insert(rResults.begin(), localResults.begin(), localResults.begin() + NumberOfResults);

    KRATOS_CATCH("")
}

} // namespace Kratos