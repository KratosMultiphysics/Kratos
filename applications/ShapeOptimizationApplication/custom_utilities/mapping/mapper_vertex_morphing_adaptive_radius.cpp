// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>

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

template<class TBaseVertexMorphingMapper>
MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::MapperVertexMorphingAdaptiveRadius(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    Parameters MapperSettings)
    : BaseType(rOriginModelPart, rDestinationModelPart, MapperSettings),
        mrOriginModelPart(rOriginModelPart),
        mrDestinationModelPart(rDestinationModelPart),
        mFilterRadiusFactor(MapperSettings["filter_radius_factor"].GetDouble())
{
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Initialize()
{
    CalculateAdaptiveVertexMorphingRadius();
    BaseType::Initialize();
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Update()
{
    CalculateAdaptiveVertexMorphingRadius();
    BaseType::Update();
}

template<class TBaseVertexMorphingMapper>
std::string MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::Info() const
{
    return BaseType::Info() + "AdaptiveRadius";
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << BaseType::Info() << "AdaptiveRadius";
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::PrintData(std::ostream& rOStream) const
{
}

template<class TBaseVertexMorphingMapper>
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

        void LocalReduce(const value_type& rGPVector)
        {
            for (auto& r_gp : rGPVector.GetContainer()) {
                this->gp_vector.push_back(r_gp);
            }
        }
        void ThreadSafeReduce(GlobalPointerAdder& rOther)
        {
#pragma omp critical
            {
                for (auto& r_gp : rOther.gp_vector.GetContainer()) {
                    this->gp_vector.push_back(r_gp);
                }
            }
        }
    };

    const auto& r_data_communicator = mrDestinationModelPart.GetCommunicator().GetDataCommunicator();

    FindNodalNeighboursForEntitiesProcess<ModelPart::ConditionsContainerType> find_neighbours_process(
        r_data_communicator,
        mrDestinationModelPart,
        NEIGHBOUR_NODES);
    find_neighbours_process.Execute();

    GlobalPointersVector<NodeType> all_global_pointers =
        block_for_each<GlobalPointerAdder>(mrDestinationModelPart.Nodes(), [&](NodeType& rNode) {
            return rNode.GetValue(NEIGHBOUR_NODES);
        });

    GlobalPointerCommunicator<NodeType> pointer_comm(r_data_communicator, all_global_pointers);

    auto coordinates_proxy = pointer_comm.Apply(
        [](const GlobalPointer<NodeType>& rGP) { return rGP->Coordinates(); });

    block_for_each(mrDestinationModelPart.Nodes(), [&](NodeType& rNode) {
        double max_distance = -1.0;
        const auto& r_coordinates_origin_node = rNode.Coordinates();
        const auto& r_neighbours = rNode.GetValue(NEIGHBOUR_NODES).GetContainer();
        for (const auto& r_neighbour : r_neighbours) {
            const auto& r_coordinates_neighbour_node = coordinates_proxy.Get(r_neighbour);
            const double distance = norm_2(r_coordinates_origin_node - r_coordinates_neighbour_node);
            max_distance = std::max(max_distance, distance);
        }

        rNode.SetValue(VERTEX_MORPHING_RADIUS, max_distance * mFilterRadiusFactor);
    });

    KRATOS_CATCH("");
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::SmoothenNeighbourBasedFilterRadius()
{

}

template<class TBaseVertexMorphingMapper>
double MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::GetVertexMorphingRadius(const NodeType& rNode) const
{
    return rNode.GetValue(VERTEX_MORPHING_RADIUS);
}

template<class TBaseVertexMorphingMapper>
void MapperVertexMorphingAdaptiveRadius<TBaseVertexMorphingMapper>::CalculateAdaptiveVertexMorphingRadius()
{
    BuiltinTimer timer;
    KRATOS_INFO("") << std::endl;
    KRATOS_INFO("ShapeOpt") << "Starting calculation of adaptive vertext morphing radius for  " << mrDestinationModelPart.FullName() << "..." << std::endl;

    CalculateNeighbourBasedFilterRadius();
    SmoothenNeighbourBasedFilterRadius();

    KRATOS_INFO("ShapeOpt") << "Finished calculation of adaptive vertext morphing radius in " << timer.ElapsedSeconds() << " s." << std::endl;
}

// template instantiations
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphing>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingMatrixFree>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingImprovedIntegration>;
template class MapperVertexMorphingAdaptiveRadius<MapperVertexMorphingSymmetric>;

}  // namespace Kratos.
