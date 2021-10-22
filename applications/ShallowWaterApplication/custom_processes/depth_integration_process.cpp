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
#include "geometries/line_3d_2.h"
#include "depth_integration_process.h"
#include "utilities/parallel_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/find_intersected_objects_utility.h"

namespace Kratos
{

DepthIntegrationProcess::DepthIntegrationProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : Process(),
        mrVolumeModelPart(rModel.GetModelPart(ThisParameters["volume_model_part_name"].GetString())),
        mrInterfaceModelPart(rModel.GetModelPart(ThisParameters["interface_model_part_name"].GetString()))
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mDirection = ThisParameters["direction_of_integration"].GetVector();
    mDirection /= norm_2(mDirection);
}

void DepthIntegrationProcess::Execute()
{
    double bottom, top;
    GetVolumePartBounds(bottom, top);
    FindIntersectedObjectsUtility intersections(mrVolumeModelPart);
    for (auto& node : mrInterfaceModelPart.Nodes()) {
        GeometryType::Pointer integration_line = CreateIntegrationLine(node, bottom, top);
        PointerVector<GeometricalObject> intersected_objects;
        intersections.FindIntersectedObjects(integration_line, intersected_objects);
        Integrate(intersected_objects, node);
    }
}

void DepthIntegrationProcess::Integrate(PointerVector<GeometricalObject>& rObjects, NodeType& rNode)
{
    array_1d<double,3> velocity = ZeroVector(3);
    double min_elevation = 1e6;
    double max_elevation = -1e6;
    int num_nodes = 0;
    for (auto& object : rObjects) {
        array_1d<double,3> obj_velocity = ZeroVector(3);
        double obj_min_elevation = 1e6;
        double obj_max_elevation = -1e6;
        int obj_num_nodes = object.GetGeometry().size();
        for (auto& node : object.GetGeometry()) {
            velocity += node.FastGetSolutionStepValue(VELOCITY);
            obj_min_elevation = std::min(obj_min_elevation, inner_prod(mDirection, node));
            obj_max_elevation = std::max(obj_max_elevation, inner_prod(mDirection, node));
        }
        velocity += obj_velocity;
        min_elevation = std::min(min_elevation, obj_min_elevation);
        max_elevation = std::max(max_elevation, obj_max_elevation);
        num_nodes += obj_num_nodes;
    }
    velocity /= num_nodes;
    double height = max_elevation - min_elevation;
    rNode.FastGetSolutionStepValue(VELOCITY) = velocity;
    rNode.FastGetSolutionStepValue(HEIGHT) = height;
}

void DepthIntegrationProcess::GetVolumePartBounds(double& rMin, double& rMax)
{
    using MultipleReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

    std::tie(rMin, rMax) = block_for_each<MultipleReduction>(mrVolumeModelPart.Nodes(), [&](NodeType& node){
        const double distance = inner_prod(mDirection, node);
        return std::make_tuple(distance, distance);
    });
}

Geometry<Node<3>>::Pointer DepthIntegrationProcess::CreateIntegrationLine(
    const NodeType& rNode,
    const double Bottom,
    const double Top)
{
    array_1d<double,3> origin = rNode - mDirection * inner_prod(rNode, mDirection);
    array_1d<double,3> start = origin + Top * mDirection;
    array_1d<double,3> end = origin + Bottom * mDirection;
    PointerVector<Node<3>> nodes;
    nodes.push_back(NodeType::Pointer(new NodeType(0, start)));
    nodes.push_back(NodeType::Pointer(new NodeType(0, end)));
    return Kratos::make_shared<Line3D2<NodeType>>(nodes);
}

int DepthIntegrationProcess::Check()
{
    KRATOS_ERROR_IF(mDirection.size() != 3) << "DepthIntegrationProcess: The direction of integration must have three coordinates." << std::endl;
    return 0;
}

const Parameters DepthIntegrationProcess::GetDefaultParameters() const
{
    auto default_parameters = Parameters(R"(
    {
        "volume_model_part_name"    : "",
        "interface_model_part_name" : "",
        "direction_of_integration"  : [0.0, 0.0, 1.0]
    })");
    return default_parameters;
}

}  // namespace Kratos.
