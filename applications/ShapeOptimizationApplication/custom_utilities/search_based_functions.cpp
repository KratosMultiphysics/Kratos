// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "search_based_functions.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// ==============================================================================

namespace Kratos
{

SearchBasedFunctions::SearchBasedFunctions( ModelPart& modelPart )
    : mrModelPart( modelPart ),
        mBucketSize(100),
        mMaxNeighborNodes(10000)
{
    // Create search tree
    mListOfNodesInModelPart.resize(mrModelPart.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesInModelPart[counter++] = pnode;
    }

    mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesInModelPart.begin(), mListOfNodesInModelPart.end(), mBucketSize));
}

void SearchBasedFunctions::FlagNodesInRadius( ModelPart::NodesContainerType& rNodeSet, const Kratos::Flags& rFlag, double filter_radius)
{
    // Mark all nodes within one filter radius away from the fixed boundary
    bool is_max_number_too_small = true;

    while(is_max_number_too_small == true)
    {
        is_max_number_too_small = false;
        for(auto& node_i : rNodeSet)
        {
            NodeVector neighbor_nodes( mMaxNeighborNodes );
            DoubleVector resulting_squared_distances( mMaxNeighborNodes, 0.0 );
            unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( node_i,
                                                                                filter_radius,
                                                                                neighbor_nodes.begin(),
                                                                                resulting_squared_distances.begin(),
                                                                                mMaxNeighborNodes );

            if(number_of_neighbors >= mMaxNeighborNodes)
            {
                KRATOS_WARNING("ShapeOpt::FlagNodesInRadius") << "For node " << node_i.Id() << " maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) is reached! Increasing maximum number by factor 2. " << std::endl;
                is_max_number_too_small = true;
                break;
            }

            // Loop over all nodes in radius (including node on damping region itself)
            for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                neighbor_node.Set(rFlag,true);
            }
        }
        mMaxNeighborNodes *=2;
    }
}

}  // namespace Kratos.
