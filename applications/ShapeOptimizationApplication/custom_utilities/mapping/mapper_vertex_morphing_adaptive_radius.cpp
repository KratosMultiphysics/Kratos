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
#include "processes/find_global_nodal_entity_neighbours_process.h"
#include "shape_optimization_application.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "custom_utilities/geometry_utilities.h"

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
        mRadiusFunctionType(MapperSettings["adaptive_filter_settings"]["radius_function"].GetString()),
        mRadiusFunctionParameter(MapperSettings["adaptive_filter_settings"]["radius_function_parameter"].GetDouble()),
        mMinimumFilterRadius(MapperSettings["adaptive_filter_settings"]["minimum_filter_radius"].GetDouble()),
        mCurvatureLimit(MapperSettings["adaptive_filter_settings"]["curvature_limit"].GetDouble()),
        mNumberOfSmoothingIterations(MapperSettings["adaptive_filter_settings"]["filter_radius_smoothing_iterations"].GetInt()),
        mMaxNumberOfNeighbors(MapperSettings["max_nodes_in_filter_radius"].GetInt())
{
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Initialize()
{
    BaseType::Initialize();

    KRATOS_INFO("ShapeOpt") << "minimum_filter_radius: " << mMinimumFilterRadius  << std::endl;
    KRATOS_INFO("ShapeOpt") << "radius_function: " << mRadiusFunctionType  << std::endl;
    if (mRadiusFunctionType == "analytic") {
        KRATOS_INFO("ShapeOpt") << "curvature_limit: " << mCurvatureLimit  << std::endl;
    } else {
        KRATOS_INFO("ShapeOpt") << "curvature_limit: Setting is ignored for the chosen radius_function!" << std::endl;
    }
    KRATOS_INFO("ShapeOpt") << "radius_function_parameter: " << mRadiusFunctionParameter  << std::endl;
    KRATOS_INFO("ShapeOpt") << "filter_radius_smoothing_iterations: " << mNumberOfSmoothingIterations  << std::endl;
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
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CurvatureFunction(const double& rCurvature, const double& rElementSize) {

    // Using the function delta_K = a / (kappa + b)
    // two equations for parameter a and b:
    // (I) a / b = delta_K(r_min)
    // (II) a / (kappa_limit + b) = mRadiusFunctionParameter
    // This is based on the "analytic" relation (found out by numerical experiment) between curvature change and filter radius: delta_K ~ 4 / r**4
    // => r = pow(4 / delta_K, 0.25)
    if (mRadiusFunctionType == "analytic")
    {
        double delta_K_max = 4 / pow(mMinimumFilterRadius/rElementSize, 4);

        double b = (mRadiusFunctionParameter * mCurvatureLimit) / (delta_K_max - mRadiusFunctionParameter);
        double a = b * delta_K_max;
        double delta_K = a / (std::abs(rCurvature) + b);

        double filter_radius = (delta_K > 0) ? rElementSize * pow(4/delta_K, 0.25) : mMinimumFilterRadius;
        return filter_radius;
    }
    // Using linear function r = r_min + a * h * kappa
    else if (mRadiusFunctionType == "linear")
    {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * std::abs(rCurvature);
    }
    // Using square root function r = r_min + a * h * sqrt(kappa)
    else if (mRadiusFunctionType == "square_root") {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * sqrt(std::abs(rCurvature));
    }
    // Using fourth root function r = r_min + a * h * kappa**0.25
    else if (mRadiusFunctionType == "fourth_root") {
        return mMinimumFilterRadius + mRadiusFunctionParameter * rElementSize * pow(std::abs(rCurvature), 0.25);
    } else {
        KRATOS_ERROR << "ShapeOpt Adaptive Filter: Curvature function type " << mRadiusFunctionType << " not supported for adaptive filter vertex morphing method." << std::endl;
    }
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateCurvatureBasedFilterRadius()
{
    KRATOS_TRY

    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting calculation of curvature based filter radius for  " << mrDestinationModelPart.FullName() << "..." << std::endl;

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

    FindGlobalNodalEntityNeighboursProcess<ModelPart::ConditionsContainerType> find_condition_neighbours_process(
        mrDestinationModelPart,
        NEIGHBOUR_CONDITIONS);
    find_condition_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrDestinationModelPart.Nodes(), [&](NodeType &rNode)
                                            { return rNode.GetValue(NEIGHBOUR_NODES); });

    GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType> &rGP)
        { return rGP->Coordinates(); });

    GeometryUtilities geom_util(this->mrDestinationModelPart);
    geom_util.CalculateGaussianCurvature();
    block_for_each(mrDestinationModelPart.Nodes(), [&](NodeType &rNode) {
        double max_distance = -1.0;
        const auto& r_coordinates_origin_node = rNode.Coordinates();
        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        for (const auto& r_neighbour : r_neighbours) {
            const auto& r_coordinates_neighbour_node = coordinates_proxy.Get(r_neighbour);
            const double distance = norm_2(r_coordinates_origin_node - r_coordinates_neighbour_node);
            max_distance = std::max(max_distance, distance);
        }
        double gaussian_curvature = rNode.FastGetSolutionStepValue(GAUSSIAN_CURVATURE);

        double vm_radius = this->CurvatureFunction(gaussian_curvature, max_distance);
        rNode.FastGetSolutionStepValue(MAX_NEIGHBOUR_DISTANCE) = max_distance;
        rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS_RAW) = vm_radius;
        rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS) = vm_radius;
    });

    KRATOS_INFO("ShapeOpt") << "Finished calculation of curvature based filter radius " << timer.ElapsedSeconds() << " s." << std::endl;

    KRATOS_CATCH("");
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::SmoothenCurvatureBasedFilterRadius()
{
    KRATOS_TRY

    const IndexType number_of_nodes = mrDestinationModelPart.NumberOfNodes();
    Vector raw_radius(number_of_nodes);
    Vector temp_radius(number_of_nodes);
    IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
        raw_radius[iNode] = (mrDestinationModelPart.NodesBegin() + iNode)->FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS_RAW);
    });

    for (IndexType iter = 0; iter < mNumberOfSmoothingIterations; ++iter) {
        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            const auto& r_node = *(mrDestinationModelPart.NodesBegin() + iNode);
            double current_radius = r_node.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);
            const double current_raw_radius = raw_radius[iNode];
            if (current_raw_radius > current_radius) {
                temp_radius[iNode] = current_raw_radius;
            } else {
                temp_radius[iNode] = current_radius;
            }
        });

        IndexPartition<IndexType>(number_of_nodes).for_each([&](const IndexType iNode) {
            auto& r_node_i = *(mrDestinationModelPart.NodesBegin() + iNode);
            double& radius = r_node_i.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS);

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
    return std::max(rNode.FastGetSolutionStepValue(VERTEX_MORPHING_RADIUS), mMinimumFilterRadius);
}

template <class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateAdaptiveVertexMorphingRadius()
{
    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting calculation of adaptive vertex morphing radius for  " << mrDestinationModelPart.FullName() << "..." << std::endl;

    AssignMappingIds();
    CreateListOfNodesInOriginModelPart();
    CreateSearchTreeWithAllNodesInOriginModelPart();
    CalculateCurvatureBasedFilterRadius();
    SmoothenCurvatureBasedFilterRadius();

    KRATOS_INFO("ShapeOpt") << "Finished calculation of adaptive vertex morphing radius in " << timer.ElapsedSeconds() << " s." << std::endl;
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
