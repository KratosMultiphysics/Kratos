//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andrea Montanino 
//

// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "write_from_sw_at_interface_process.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TDim>
const Parameters WriteFromSwAtInterfaceProcess<TDim>::GetDefaultParameters() const
{
    auto default_parameters = Parameters(R"(
    {
        "volume_model_part_name"    : "",
        "interface_model_part_name" : "",
        "store_historical_database" : false,
        "extrapolate_boundaries"    : false,
        "print_velocity_profile"    : false
    })");
    return default_parameters;
}

template<std::size_t TDim>
WriteFromSwAtInterfaceProcess<TDim>::WriteFromSwAtInterfaceProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : Process(),
        mrVolumeModelPart(rModel.GetModelPart(ThisParameters["volume_model_part_name"].GetString())),
        mrInterfaceModelPart(rModel.GetModelPart(ThisParameters["interface_model_part_name"].GetString()))
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mStoreHistorical = ThisParameters["store_historical_database"].GetBool();
    mExtrapolateBoundaries = ThisParameters["extrapolate_boundaries"].GetBool();
    mDirection = -mrVolumeModelPart.GetProcessInfo()[GRAVITY];
    mDirection /= norm_2(mDirection);
    mPrintVelocityProfile = ThisParameters["print_velocity_profile"].GetBool();

    if (!mStoreHistorical) {
        VariableUtils().SetNonHistoricalVariableToZero(MOMENTUM, mrInterfaceModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(VELOCITY, mrInterfaceModelPart.Nodes());
        VariableUtils().SetNonHistoricalVariableToZero(HEIGHT, mrInterfaceModelPart.Nodes());
    }

    // if (mExtrapolateBoundaries) {
    //     FindBoundaryNeighbors();
    // }
}

template<std::size_t TDim>
void WriteFromSwAtInterfaceProcess<TDim>::Execute()
{
    
    BinBasedFastPointLocator<TDim> locator(mrVolumeModelPart);
    locator.UpdateSearchDatabase();

    struct locator_tls {
        Vector N;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results;
        locator_tls(const int max_results = 10000) {
            N.resize(TDim+1);
            results.resize(max_results);
        }
    };
    // If you do an integral, why you are doing this for all nodes of the interface? AMontanino87
    block_for_each(mrInterfaceModelPart.Nodes(), locator_tls(), [&](NodeType& rNode, locator_tls& rTLS){
        ReadAndSetValues(rNode, locator, rTLS.results);
    });

    if (mExtrapolateBoundaries) {
        CopyValues(*mpFirstBoundaryNeighbor, *mpFirstBoundaryNode);
        CopyValues(*mpSecondBoundaryNeighbor, *mpSecondBoundaryNode);
    }
}

// template<std::size_t TDim>
// void WriteFromSwAtInterfaceProcess<TDim>::GetBoundingVolumeLimits(double& rMin, double& rMax)
// {
//     using MultipleReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

//     std::tie(rMin, rMax) = block_for_each<MultipleReduction>(mrVolumeModelPart.Nodes(), [&](NodeType& node){
//         const double distance = inner_prod(mDirection, node);
//         return std::make_tuple(distance, distance);
//     });
// }

template<std::size_t TDim>
void WriteFromSwAtInterfaceProcess<TDim>::ReadAndSetValues(
    NodeType& rNode,
    BinBasedFastPointLocator<TDim>& rLocator,
    typename BinBasedFastPointLocator<TDim>::ResultContainerType& rResults)
{
    
    // initialize the integration variables
    Element::Pointer p_elem;
    array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> momentum = ZeroVector(3);
    double height = 0.0;
    
    velocity = rNode.FastGetSolutionStepValue(VELOCITY);
    momentum = rNode.FastGetSolutionStepValue(MOMENTUM);
    height = rNode.FastGetSolutionStepValue(HEIGHT);

    SetValue(rNode, MOMENTUM, momentum);
    SetValue(rNode, VELOCITY, velocity);
    SetValue(rNode, HEIGHT, height);

}

// template<std::size_t TDim>
// array_1d<double,3> WriteFromSwAtInterfaceProcess<TDim>::InterpolateVelocity(
//     const Element::Pointer pElement,
//     const Vector& rShapeFunctionValues) const
// {
//     KRATOS_DEBUG_ERROR_IF(pElement->GetGeometry().size() != rShapeFunctionValues.size()) << "WriteFromSwAtInterfaceProcess: check the found element!" << std::endl;
//     array_1d<double,3> velocity = ZeroVector(3);
//     int n = 0;
//     for (auto& r_node : pElement->GetGeometry()) {
//         velocity += rShapeFunctionValues[n] * r_node.FastGetSolutionStepValue(VELOCITY);
//         n++;
//     }
//     return velocity;
// }

template<std::size_t TDim>
int WriteFromSwAtInterfaceProcess<TDim>::Check()
{
    const auto dimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(dimension != 2 && dimension != 3) << Info() << ": Wrong DOMAIN_SIZE equal to " << dimension << "in model part " << mrVolumeModelPart.Name() << std::endl;
    KRATOS_ERROR_IF(dimension == 2 && mExtrapolateBoundaries) << Info() << ": Is not possible to extrapolate the boundaries in a 2D simulation." << std::endl;
    KRATOS_ERROR_IF(mrVolumeModelPart.NumberOfNodes() == 0) << Info() << ": The volume model part is empty. Not possible to construct the search structure." << std::endl;
    return 0;
}

// template<std::size_t TDim>
// void WriteFromSwAtInterfaceProcess<TDim>::FindBoundaryNeighbors()
// {
//     // Step 1, find the center of the nodes
//     const int num_nodes = mrInterfaceModelPart.NumberOfNodes();
//     array_1d<double,3> center = block_for_each<SumReduction<array_1d<double,3>>>(
//         mrInterfaceModelPart.Nodes(), [&](NodeType& rNode){return rNode.Coordinates();}
//     );
//     center /= num_nodes;

//     // Step 2, compute the distances from the center
//     std::vector<double> distances(num_nodes);
//     IndexPartition<int>(num_nodes).for_each([&](int i){
//         auto it_node = mrInterfaceModelPart.NodesBegin() + i;
//         double distance = norm_2(center - *it_node);
//         distances[i] = distance;
//     });

//     // Step 3, the two further nodes are the boundaries
//     std::size_t i_first_node = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
//     double max_distance = distances[i_first_node];
//     distances[i_first_node] = 0.0;
//     std::size_t i_second_node = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
//     mpFirstBoundaryNode = &*(mrInterfaceModelPart.NodesBegin() + i_first_node);
//     mpSecondBoundaryNode = &*(mrInterfaceModelPart.NodesBegin() + i_second_node);

//     // Step 4, compute the distances from the first node, the closer is the neighbor
//     IndexPartition<int>(num_nodes).for_each([&](int i){
//         auto it_node = mrInterfaceModelPart.NodesBegin() + i;
//         double distance = norm_2(*mpFirstBoundaryNode - *it_node);
//         if (distance < 1e-6) {
//             distance = max_distance; // this is the boundary itself
//         }
//         distances[i] = distance;
//     });
//     std::size_t i_first_neigh = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
//     mpFirstBoundaryNeighbor = &*(mrInterfaceModelPart.NodesBegin() + i_first_neigh);

//     // Step 5, compute the distances from the second node, the closer is the neighbor
//     IndexPartition<int>(num_nodes).for_each([&](int i){
//         auto it_node = mrInterfaceModelPart.NodesBegin() + i;
//         double distance = norm_2(*mpSecondBoundaryNode - *it_node);
//         if (distance < 1e-6) {
//             distance = max_distance; // this is the boundary itself
//         }
//         distances[i] = distance;
//     });
//     std::size_t i_second_neigh = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
//     mpSecondBoundaryNeighbor = &*(mrInterfaceModelPart.NodesBegin() + i_second_neigh);
// }

template<std::size_t TDim>
void WriteFromSwAtInterfaceProcess<TDim>::CopyValues(const NodeType& rOriginNode, NodeType& rDestinationNode)
{
    SetValue(rDestinationNode, HEIGHT, GetValue<double>(rOriginNode, HEIGHT));
    SetValue(rDestinationNode, VELOCITY, GetValue<array_1d<double,3>>(rOriginNode, VELOCITY));
    SetValue(rDestinationNode, MOMENTUM, GetValue<array_1d<double,3>>(rOriginNode, MOMENTUM));
}

template class WriteFromSwAtInterfaceProcess<2>;
template class WriteFromSwAtInterfaceProcess<3>;

}  // namespace Kratos.
