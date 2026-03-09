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
        VariableUtils().SetNonHistoricalVariableToZero(VERTICAL_VELOCITY, mrInterfaceModelPart.Nodes());
    }

}

template<std::size_t TDim>
void WriteFromSwAtInterfaceProcess<TDim>::Execute()
{
    
    BinBasedFastPointLocator<TDim> locator(mrVolumeModelPart);
    locator.UpdateSearchDatabase();

    block_for_each(mrInterfaceModelPart.Nodes(), locator_tls(), [&](NodeType& rNode, locator_tls& rTLS){
        ReadAndSetValues(rNode, locator, rTLS.results);
    });

    if (mExtrapolateBoundaries) {
        CopyValues(*mpFirstBoundaryNeighbor, *mpFirstBoundaryNode);
        CopyValues(*mpSecondBoundaryNeighbor, *mpSecondBoundaryNode);
    }
}

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
    double vertical_velocity = 0.0;
    double topography = 0.0;
    
    noalias(velocity) = rNode.FastGetSolutionStepValue(VELOCITY);
    noalias(momentum) = rNode.FastGetSolutionStepValue(MOMENTUM);
    height = rNode.FastGetSolutionStepValue(HEIGHT);
    vertical_velocity = rNode.FastGetSolutionStepValue(VERTICAL_VELOCITY);
    topography = rNode.FastGetSolutionStepValue(TOPOGRAPHY);

    SetValue(rNode, MOMENTUM, momentum);
    SetValue(rNode, VELOCITY, velocity);
    SetValue(rNode, HEIGHT, height);
    SetValue(rNode, VERTICAL_VELOCITY, vertical_velocity);
    SetValue(rNode, TOPOGRAPHY, topography);

}

template<std::size_t TDim>
int WriteFromSwAtInterfaceProcess<TDim>::Check()
{
    const auto dimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(dimension != 2 && dimension != 3) << Info() << ": Wrong DOMAIN_SIZE equal to " << dimension << "in model part " << mrVolumeModelPart.Name() << std::endl;
    KRATOS_ERROR_IF(dimension == 2 && mExtrapolateBoundaries) << Info() << ": Is not possible to extrapolate the boundaries in a 2D simulation." << std::endl;
    KRATOS_ERROR_IF(mrVolumeModelPart.NumberOfNodes() == 0) << Info() << ": The volume model part is empty. Not possible to construct the search structure." << std::endl;
    return 0;
}


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
