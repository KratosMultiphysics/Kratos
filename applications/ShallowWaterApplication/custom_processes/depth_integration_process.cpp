//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/model_part.h"
#include "depth_integration_process.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

template<std::size_t TDim>
const Parameters DepthIntegrationProcess<TDim>::GetDefaultParameters() const
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
DepthIntegrationProcess<TDim>::DepthIntegrationProcess(
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

    if (mExtrapolateBoundaries) {
        FindBoundaryNeighbors();
    }
}

template<std::size_t TDim>
void DepthIntegrationProcess<TDim>::Execute()
{
    double bottom, top;
    GetBoundingVolumeLimits(bottom, top);
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
    
    block_for_each(mrInterfaceModelPart.Nodes(), locator_tls(), [&](NodeType& rNode, locator_tls& rTLS){
        Integrate(rNode, bottom, top, locator, rTLS.results, rTLS.N);
    });

    if (mExtrapolateBoundaries) {
        CopyValues(*mpFirstBoundaryNeighbor, *mpFirstBoundaryNode);
        CopyValues(*mpSecondBoundaryNeighbor, *mpSecondBoundaryNode);
    }
}

template<std::size_t TDim>
void DepthIntegrationProcess<TDim>::GetBoundingVolumeLimits(double& rMin, double& rMax)
{
    using MultipleReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

    std::tie(rMin, rMax) = block_for_each<MultipleReduction>(mrVolumeModelPart.Nodes(), [&](NodeType& node){
        const double distance = inner_prod(mDirection, node);
        return std::make_tuple(distance, distance);
    });
}

template<std::size_t TDim>
void DepthIntegrationProcess<TDim>::Integrate(
    NodeType& rNode,
    const double Bottom,
    const double Top,
    BinBasedFastPointLocator<TDim>& rLocator,
    typename BinBasedFastPointLocator<TDim>::ResultContainerType& rResults,
    Vector& rShapeFunctionsValues)
{
    // Integrate the velocity over the water column
    const int num_steps = 50;
    const double step = (Top - Bottom) / double(num_steps-1);
    const array_1d<double,3> start = rNode + mDirection * (Bottom -inner_prod(rNode, mDirection));

    // initialize the integration variables
    Element::Pointer p_elem;
    array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> momentum = ZeroVector(3);
    array_1d<double,3> elem_velocity;
    double height = 0.0;
    double bottom_depth = 0.0;
    double top_depth = 0.0;

    // perform the first search to identify the bounds
    array_1d<double,3> point(start);
    bool prev_is_found = false;
    std::size_t last_found_id = 0;
    std::size_t num_found = 0;
    for (int i = 0; i < num_steps; ++i) {
        bool is_found = rLocator.FindPointOnMesh(point, rShapeFunctionsValues, p_elem, rResults.begin());
        if (is_found) {
            if (p_elem->Id() != last_found_id) {
                last_found_id = p_elem->Id();
                num_found++;
            }
            top_depth = inner_prod(point, mDirection);
        }
        if (!prev_is_found && is_found) { // this is the first element
            bottom_depth = inner_prod(point, mDirection);
        }
        if (prev_is_found && !is_found) { // the previous element is the last
            break;
        }
        prev_is_found = is_found;
        point += step * mDirection;
    }
    height = top_depth - bottom_depth;

    // perform the custom quadrature
    double average_weight = 1.0 / double(num_found);
    array_1d<double,3> integration_step = average_weight * height * mDirection;
    array_1d<double,3> integration_point = rNode + (bottom_depth -inner_prod(rNode, mDirection)) * mDirection;
    for (std::size_t i = 0; i < num_found; ++i) { // the custom quadrature
        integration_point += integration_step;
        bool is_found = rLocator.FindPointOnMesh(integration_point, rShapeFunctionsValues, p_elem, rResults.begin());
        if (is_found) {
            array_1d<double,3> point_velocity = InterpolateVelocity(p_elem, rShapeFunctionsValues);
            momentum += average_weight * height * point_velocity;
            velocity += average_weight * point_velocity;
        }
        else {
            KRATOS_WARNING("Depth integration") << "Point not found inside the fluid domain!!!" << std::endl;
        }
    }

    SetValue(rNode, MOMENTUM, momentum);
    SetValue(rNode, VELOCITY, velocity);
    SetValue(rNode, HEIGHT, height);

    // print the debug log if specified
    if (mPrintVelocityProfile)
    {
        // initialize the file
        const double time = mrVolumeModelPart.GetProcessInfo()[TIME];
        std::ostringstream file_name;
        file_name << std::fixed;
        file_name << std::setprecision(2);
        file_name << "depth_integration/vel_" << rNode.X() << "_" << time << ".dat";
        std::ofstream log_file(file_name.str());
        KRATOS_ERROR_IF_NOT(log_file.is_open()) << "Unable to open the log file. Make sure the \"depth_integration\" folder exists" << std::endl;

        // store the mean velocity
        std::stringstream vel_mean;
        vel_mean << bottom_depth << '\t' << velocity[0] << std::endl;
        vel_mean << top_depth    << '\t' << velocity[0] << std::endl;

        // compute the beta velocity
        std::stringstream vel_beta;
        array_1d<double,3> v_beta = ZeroVector(3);
        const double beta_depth = top_depth -0.531 * height;
        const array_1d<double,3> beta_point(rNode + mDirection * (beta_depth - inner_prod(mDirection, rNode)));
        bool is_found = rLocator.FindPointOnMesh(beta_point, rShapeFunctionsValues, p_elem, rResults.begin());
        if (is_found) {
            v_beta = InterpolateVelocity(p_elem, rShapeFunctionsValues);
        }
        vel_beta << beta_depth << '\t' << v_beta[0] << std::endl;

        // compute the velocity profile
        std::stringstream vel_profile;
        point = start;
        for (int i = 1; i < num_steps; ++i) {
            is_found = rLocator.FindPointOnMesh(point, rShapeFunctionsValues, p_elem, rResults.begin());
            if (is_found) {
                elem_velocity = InterpolateVelocity(p_elem, rShapeFunctionsValues);
                vel_profile << inner_prod(point, mDirection) << '\t' << elem_velocity[0] << std::endl;
            }
            else {
                vel_profile << inner_prod(point, mDirection) << '\t' << 0.0 << std::endl;
            }
            point += step * mDirection;
        }

        // print to file
        log_file << "First two lines: mean velocity  |  Third line: beta velocity  |  Other: velocity profile" << std::endl;
        log_file << "# z \t u" << std::endl;
        log_file << vel_mean.str();
        log_file << vel_beta.str();
        log_file << vel_profile.str();
        log_file.close();
    }
}

template<std::size_t TDim>
array_1d<double,3> DepthIntegrationProcess<TDim>::InterpolateVelocity(
    const Element::Pointer pElement,
    const Vector& rShapeFunctionValues) const
{
    KRATOS_DEBUG_ERROR_IF(pElement->GetGeometry().size() != rShapeFunctionValues.size()) << "DepthIntegrationProcess: check the found element!" << std::endl;
    array_1d<double,3> velocity = ZeroVector(3);
    int n = 0;
    for (auto& r_node : pElement->GetGeometry()) {
        velocity += rShapeFunctionValues[n] * r_node.FastGetSolutionStepValue(VELOCITY);
        n++;
    }
    return velocity;
}

template<std::size_t TDim>
int DepthIntegrationProcess<TDim>::Check()
{
    const auto dimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(dimension != 2 && dimension != 3) << Info() << ": Wrong DOMAIN_SIZE equal to " << dimension << "in model part " << mrVolumeModelPart.Name() << std::endl;
    KRATOS_ERROR_IF(dimension == 2 && mExtrapolateBoundaries) << Info() << ": Is not possible to extrapolate the boundaries in a 2D simulation." << std::endl;
    KRATOS_ERROR_IF(mrVolumeModelPart.NumberOfNodes() == 0) << Info() << ": The volume model part is empty. Not possible to construct the search structure." << std::endl;
    return 0;
}

template<std::size_t TDim>
void DepthIntegrationProcess<TDim>::FindBoundaryNeighbors()
{
    // Step 1, find the center of the nodes
    const int num_nodes = mrInterfaceModelPart.NumberOfNodes();
    array_1d<double,3> center = block_for_each<SumReduction<array_1d<double,3>>>(
        mrInterfaceModelPart.Nodes(), [&](NodeType& rNode){return rNode.Coordinates();}
    );
    center /= num_nodes;

    // Step 2, compute the distances from the center
    std::vector<double> distances(num_nodes);
    IndexPartition<int>(num_nodes).for_each([&](int i){
        auto it_node = mrInterfaceModelPart.NodesBegin() + i;
        double distance = norm_2(center - *it_node);
        distances[i] = distance;
    });

    // Step 3, the two further nodes are the boundaries
    std::size_t i_first_node = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
    double max_distance = distances[i_first_node];
    distances[i_first_node] = 0.0;
    std::size_t i_second_node = std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));
    mpFirstBoundaryNode = &*(mrInterfaceModelPart.NodesBegin() + i_first_node);
    mpSecondBoundaryNode = &*(mrInterfaceModelPart.NodesBegin() + i_second_node);

    // Step 4, compute the distances from the first node, the closer is the neighbor
    IndexPartition<int>(num_nodes).for_each([&](int i){
        auto it_node = mrInterfaceModelPart.NodesBegin() + i;
        double distance = norm_2(*mpFirstBoundaryNode - *it_node);
        if (distance < 1e-6) {
            distance = max_distance; // this is the boundary itself
        }
        distances[i] = distance;
    });
    std::size_t i_first_neigh = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
    mpFirstBoundaryNeighbor = &*(mrInterfaceModelPart.NodesBegin() + i_first_neigh);

    // Step 5, compute the distances from the second node, the closer is the neighbor
    IndexPartition<int>(num_nodes).for_each([&](int i){
        auto it_node = mrInterfaceModelPart.NodesBegin() + i;
        double distance = norm_2(*mpSecondBoundaryNode - *it_node);
        if (distance < 1e-6) {
            distance = max_distance; // this is the boundary itself
        }
        distances[i] = distance;
    });
    std::size_t i_second_neigh = std::distance(distances.begin(), std::min_element(distances.begin(), distances.end()));
    mpSecondBoundaryNeighbor = &*(mrInterfaceModelPart.NodesBegin() + i_second_neigh);
}

template<std::size_t TDim>
void DepthIntegrationProcess<TDim>::CopyValues(const NodeType& rOriginNode, NodeType& rDestinationNode)
{
    SetValue(rDestinationNode, HEIGHT, GetValue<double>(rOriginNode, HEIGHT));
    SetValue(rDestinationNode, VELOCITY, GetValue<array_1d<double,3>>(rOriginNode, VELOCITY));
    SetValue(rDestinationNode, MOMENTUM, GetValue<array_1d<double,3>>(rOriginNode, MOMENTUM));
}

template class DepthIntegrationProcess<2>;
template class DepthIntegrationProcess<3>;

}  // namespace Kratos.
