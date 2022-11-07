// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    David Schmölz
//                   Suneth Warnakulasuriya
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>
#include <unordered_map>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/model_part.h"
#include "processes/find_global_nodal_neighbours_for_entities_process.h"
#include "shape_optimization_application.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"

// ------------------------------------------------------------------------------
// Base mapper includes
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing.h"
#include "mapper_vertex_morphing_improved_integration.h"
#include "mapper_vertex_morphing_matrix_free.h"
#include "mapper_vertex_morphing_symmetric.h"

// ------------------------------------------------------------------------------
// Base class include
// ------------------------------------------------------------------------------
#include "mapper_vertex_morphing_adaptive_radius.h"

// ==============================================================================

namespace Kratos
{

template <class TBaseVertexMorphingMapper>
MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::MapperVertexMorphingAdaptiveRadius(
    ModelPart &rOriginModelPart,
    ModelPart &rDestinationModelPart,
    Parameters MapperSettings)
    : BaseType(rOriginModelPart, rDestinationModelPart, MapperSettings),
        mrOriginModelPart(rOriginModelPart),
        mrDestinationModelPart(rDestinationModelPart),
        mFilterRadiusFactor(MapperSettings["adaptive_filter_settings"]["filter_radius_factor"].GetDouble()),
        mMinimumFilterRadius(MapperSettings["adaptive_filter_settings"]["minimum_filter_radius"].GetDouble()),
        mNumberOfSmoothingIterations(MapperSettings["adaptive_filter_settings"]["filter_radius_smoothing_iterations"].GetInt()),
        mMaxNumberOfNeighbors(MapperSettings["max_nodes_in_filter_radius"].GetInt())
{
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Initialize()
{
    BaseType::Initialize();
    CalculateAdaptiveVertexMorphingRadius();

    KRATOS_INFO("ShapeOpt") << "filter_radius_factor:  " << mFilterRadiusFactor  << std::endl;
    KRATOS_INFO("ShapeOpt") << "minimum_filter_radius:  " << mMinimumFilterRadius  << std::endl;
    KRATOS_INFO("ShapeOpt") << "filter_radius_smoothing_iterations:  " << mNumberOfSmoothingIterations  << std::endl;

}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Update()
{
    BaseType::Update();
    CalculateAdaptiveVertexMorphingRadius();
}

template <class TBaseVertexMorphingMapper>
std::string MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Info() const
{
    return BaseType::Info() + "AdaptiveRadius";
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintInfo(std::ostream &rOStream) const
{
    rOStream << BaseType::Info() << "AdaptiveRadius";
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintData(std::ostream &rOStream) const
{
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateNeighbourBasedFilterRadius()
{
    KRATOS_TRY

    class GlobalPointerAdder
    {
    public:
        typedef GlobalPointersVector<NodeType> value_type;
        typedef GlobalPointersVector<NodeType> return_type;

        return_type gp_vector;
        return_type GetValue()
        {
            gp_vector.Unique();
            return gp_vector;
        }

        void LocalReduce(const value_type &rGPVector)
        {
            for (auto &r_gp : rGPVector.GetContainer())
            {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder &rOther)
        {
#pragma omp critical
            {
                for (auto &r_gp : rOther.gp_vector.GetContainer())
                {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    const auto &r_data_communicator = mrDestinationModelPart.GetCommunicator().GetDataCommunicator();

    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrDestinationModelPart.Nodes(), [&](NodeType &rNode)
                                            { return rNode.GetValue(NEIGHBOUR_NODES); });

    GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType> &rGP)
        { return rGP->Coordinates(); });

    block_for_each(mrDestinationModelPart.Nodes(), [&](NodeType &rNode) {
        double max_distance = -1.0;
        const auto& r_coordinates_origin_node = rNode.Coordinates();
        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        for (const auto& r_neighbour : r_neighbours) {
            const auto& r_coordinates_neighbour_node = coordinates_proxy.Get(r_neighbour);
            const double distance = norm_2(r_coordinates_origin_node - r_coordinates_neighbour_node);
            max_distance = std::max(max_distance, distance);
        }

        KRATOS_INFO("ShapeOpt") << "max distance of node id = " << rNode.Id() << " : " << max_distance << std::endl;
        rNode.SetValue(VERTEX_MORPHING_RADIUS, max_distance * mFilterRadiusFactor);
        KRATOS_INFO("ShapeOpt") << "VERTEX_MORPHING_RADIUS on node id = " << rNode.Id() << " : " << max_distance * mFilterRadiusFactor << std::endl;
    });

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::SmoothenNeighbourBasedFilterRadius()
{
    KRATOS_TRY

    const IndexType number_of_nodes = mrDestinationModelPart.NumberOfNodes();
    Vector raw_radius(number_of_nodes);
    Vector temp_radius(number_of_nodes);
    IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
        raw_radius[iNode] = (mrDestinationModelPart.NodesBegin() + iNode)->GetValue(VERTEX_MORPHING_RADIUS);
    });

    for (IndexType iter = 0; iter < mNumberOfSmoothingIterations; ++iter) {
        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            const auto& r_node = *(mrDestinationModelPart.NodesBegin() + iNode);
            double current_radius = r_node.GetValue(VERTEX_MORPHING_RADIUS);
            const double current_raw_radius = raw_radius[iNode];
            if (current_raw_radius > current_radius) {
                temp_radius[iNode] = current_raw_radius;
            } else {
                temp_radius[iNode] = current_radius;
            }
        });

        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            auto& r_node_i = *(mrDestinationModelPart.NodesBegin() + iNode);
            double& radius = r_node_i.GetValue(VERTEX_MORPHING_RADIUS);

            NodeVector neighbor_nodes(mMaxNumberOfNeighbors);
            std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);

            const IndexType number_of_neighbors = mpSearchTree->SearchInRadius(
                                                    r_node_i,
                                                    radius,
                                                    neighbor_nodes.begin(),
                                                    resulting_squared_distances.begin(),
                                                    mMaxNumberOfNeighbors);

            std::vector<double> list_of_weights(number_of_neighbors, 0.0);
            double sum_of_weights = 0.0;
            this->ComputeWeightForAllNeighbors(r_node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
            radius = 0.0;

            for(IndexType neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; ++neighbor_itr) {
                const NodeType& r_node_j = *neighbor_nodes[neighbor_itr];
                const double radius_j = temp_radius[r_node_j.GetValue(MAPPING_ID)];
                radius += radius_j* list_of_weights[neighbor_itr] / sum_of_weights;
            }
        });
    }

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::GetVertexMorphingRadius(const NodeType &rNode) const
{
    KRATOS_INFO("ShapeOpt") << "VERTEX_MORPHING_RADIUS on node id = " << rNode.Id() << " : " << rNode.GetValue(VERTEX_MORPHING_RADIUS) << std::endl;
    return std::max(rNode.GetValue(VERTEX_MORPHING_RADIUS), mMinimumFilterRadius);
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateAdaptiveVertexMorphingRadius()
{
    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting calculation of adaptive vertext morphing radius for  " << mrDestinationModelPart.FullName() << "..." << std::endl;

    AssignMappingIds();
    CreateListOfNodesInOriginModelPart();
    CreateSearchTreeWithAllNodesInOriginModelPart();
    CalculateNeighbourBasedFilterRadius();
    SmoothenNeighbourBasedFilterRadius();

    KRATOS_INFO("ShapeOpt") << "Finished calculation of adaptive vertext morphing radius in " << timer.ElapsedSeconds() << " s." << std::endl;
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CreateSearchTreeWithAllNodesInOriginModelPart()
{
    BuiltinTimer timer;
    KRATOS_INFO("ShapeOpt") << "Creating search tree to perform mapping..." << std::endl;
    mpSearchTree = Kratos::make_unique<KDTree>(mListOfNodesInOriginModelPart.begin(), mListOfNodesInOriginModelPart.end(), mBucketSize);
    KRATOS_INFO("ShapeOpt") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CreateListOfNodesInOriginModelPart()
{
    mListOfNodesInOriginModelPart.resize(mrOriginModelPart.Nodes().size());
    int counter = 0;
    for (ModelPart::NodesContainerType::iterator node_it = mrOriginModelPart.NodesBegin(); node_it != mrOriginModelPart.NodesEnd(); ++node_it)
    {
        NodeTypePointer pnode = *(node_it.base());
        mListOfNodesInOriginModelPart[counter++] = pnode;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::ComputeWeightForAllNeighbors(
    const ModelPart::NodeType& destination_node,
    const NodeVector& neighbor_nodes,
    const unsigned int number_of_neighbors,
    std::vector<double>& list_of_weights,
    double& sum_of_weights)
{
    for(unsigned int neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++) {
        const NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
        const double weight = this->mpFilterFunction->ComputeWeight( destination_node.Coordinates(), neighbor_node.Coordinates(), GetVertexMorphingRadius(destination_node) );

        list_of_weights[neighbor_itr] = weight;
        sum_of_weights += weight;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::AssignMappingIds()
{
    unsigned int i = 0;

    // Note: loop in the same order as in AllocateMatrix(), to avoid reallocations of the matrix.
    for(auto& node_i : mrOriginModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);

    i = 0;
    for(auto& node_i : mrDestinationModelPart.Nodes())
        node_i.SetValue(MAPPING_ID,i++);
}

// template instantiations
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>;

} // namespace Kratos.
