// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Schmoelz David, https://github.com/dschmoelz
//
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "custom_utilities/damping/thickness_damping_utilities.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/filter_function.h"
#include "shape_optimization_application.h"

// ==============================================================================

namespace Kratos
{

ThicknessDampingUtilities::ThicknessDampingUtilities( ModelPart& modelPartToDamp, Parameters DampingSettings )
    : mrModelPartToDamp( modelPartToDamp ),
      mDampingSettings( DampingSettings ),
      mMaxNeighborNodes( DampingSettings["max_neighbor_nodes"].GetInt() )
{
    Parameters default_parameters( R"(
        {
            "sub_model_part_name"   : "MODEL_PART_NAME",
            "damping_function_type" : "cosine",
            "damping_radius"        : -1.0
        }  )" );

    // Loop over all regions for which damping is to be applied
    for (auto& r_region_parameters : mDampingSettings["damping_regions"])
    {
        r_region_parameters.ValidateAndAssignDefaults(default_parameters);
        KRATOS_ERROR_IF(r_region_parameters["damping_radius"].GetDouble() < 0.0) << "ThicknessDampingUtilities: 'damping_radius' is a mandatory setting and has to be > 0.0!" << std::endl;
    }

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ThickOpt") << "Creating search tree to perform damping..." << std::endl;
    CreateListOfNodesOfModelPart();
    CreateSearchTreeWithAllNodesOfModelPart();
    KRATOS_INFO("ThickOpt") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
    InitalizeDampingFactorsToHaveNoInfluence();
    SetDampingFactorsForAllDampingRegions();
}

void ThicknessDampingUtilities::CreateListOfNodesOfModelPart()
{
    mListOfNodesOfModelPart.resize(mrModelPartToDamp.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrModelPartToDamp.NodesBegin(); node_it != mrModelPartToDamp.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesOfModelPart[counter++] = pnode;
    }
}

void ThicknessDampingUtilities::CreateSearchTreeWithAllNodesOfModelPart()
{
    mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mListOfNodesOfModelPart.begin(), mListOfNodesOfModelPart.end(), mBucketSize));
}

void ThicknessDampingUtilities::InitalizeDampingFactorsToHaveNoInfluence()
{
    for(auto& node_i : mrModelPartToDamp.Nodes())
    {
        node_i.SetValue(DAMPING_FACTOR_T,1.0);
    }
}

void ThicknessDampingUtilities::SetDampingFactorsForAllDampingRegions()
{
    KRATOS_INFO("") << std::endl;

    BuiltinTimer timer;
    KRATOS_INFO("ThickOpt") << "Starting to prepare damping..." << std::endl;

    // Loop over all regions for which damping is to be applied
    for (const auto& r_region_parameters : mDampingSettings["damping_regions"])
    {
        const std::string& sub_model_part_name = r_region_parameters["sub_model_part_name"].GetString();
        ModelPart& dampingRegion = mrModelPartToDamp.GetRootModelPart().GetSubModelPart(sub_model_part_name);

        const auto& dampingFunctionType = r_region_parameters["damping_function_type"].GetString();
        const double dampingRadius = r_region_parameters["damping_radius"].GetDouble();

        const auto p_damping_function = CreateDampingFunction(dampingFunctionType );

        // Loop over all nodes in specified damping sub-model part
        block_for_each(dampingRegion.Nodes(), [&](const ModelPart::NodeType& rNode) {
            NodeVector neighbor_nodes( mMaxNeighborNodes );
            const unsigned int number_of_neighbors = mpSearchTree->SearchInRadius( rNode,
                                                                                dampingRadius,
                                                                                neighbor_nodes.begin(),
                                                                                mMaxNeighborNodes );

            ThrowWarningIfNodeNeighborsExceedLimit( rNode, number_of_neighbors );

            // Loop over all nodes in radius (including node on damping region itself)
            for(unsigned int j_itr = 0 ; j_itr<number_of_neighbors ; j_itr++)
            {
                ModelPart::NodeType& neighbor_node = *neighbor_nodes[j_itr];
                const double dampingFactor = 1.0 - p_damping_function->ComputeWeight( rNode.Coordinates(), neighbor_node.Coordinates(), dampingRadius);

                // For every specified damping direction we check if new damping factor is smaller than the assigned one for current node.
                // In case yes, we overwrite the value. This ensures that the damping factor of a node is computed by its closest distance to the damping region
                double& damping_factor_variable = neighbor_node.GetValue(DAMPING_FACTOR_T);
                neighbor_node.SetLock();
                if(dampingFactor < damping_factor_variable)
                    damping_factor_variable = dampingFactor;
                neighbor_node.UnSetLock();
            }
        });
    }

    KRATOS_INFO("ThickOpt") << "Finished preparation of damping in: " << timer.ElapsedSeconds() << " s" << std::endl;
}

FilterFunction::Pointer ThicknessDampingUtilities::CreateDampingFunction( std::string damping_type ) const
{
    return Kratos::make_unique<FilterFunction>(damping_type);
}

void ThicknessDampingUtilities::ThrowWarningIfNodeNeighborsExceedLimit( const ModelPart::NodeType& given_node, const unsigned int number_of_neighbors ) const
{
    if(number_of_neighbors >= mMaxNeighborNodes)
        KRATOS_WARNING("ThickOpt::ThicknessDampingUtilities") << "For node " << given_node.Id() << " and specified damping radius, maximum number of neighbor nodes (=" << mMaxNeighborNodes << " nodes) reached!" << std::endl;
}

void ThicknessDampingUtilities::DampNodalVariable( const Variable<double> &rNodalVariable )
{
    block_for_each(mrModelPartToDamp.Nodes(), [&](ModelPart::NodeType& rNode) {
        const double& damping_factor = rNode.GetValue(DAMPING_FACTOR_T);
        double& nodalVariable = rNode.FastGetSolutionStepValue(rNodalVariable);

        nodalVariable *= damping_factor;
    });
}

} // namespace Kratos
