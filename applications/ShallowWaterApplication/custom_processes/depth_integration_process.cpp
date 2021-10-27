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
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "depth_integration_process.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "shallow_water_application_variables.h"
#include "custom_utilities/find_intersected_objects_utility.h"

namespace Kratos
{

const Parameters DepthIntegrationProcess::GetDefaultParameters() const
{
    auto default_parameters = Parameters(R"(
    {
        "volume_model_part_name"    : "",
        "interface_model_part_name" : "",
        "direction_of_integration"  : [0.0, 0.0, 1.0],
        "store_historical_database" : false
    })");
    return default_parameters;
}

DepthIntegrationProcess::DepthIntegrationProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : Process(),
        mrVolumeModelPart(rModel.GetModelPart(ThisParameters["volume_model_part_name"].GetString())),
        mrInterfaceModelPart(rModel.GetModelPart(ThisParameters["interface_model_part_name"].GetString()))
{
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mStoreHistorical = ThisParameters["store_historical_database"].GetBool();
    mDirection = ThisParameters["direction_of_integration"].GetVector();
    mDirection /= norm_2(mDirection);

    if (rModel.HasModelPart("integration_auxiliary_model_part")) {
        mpIntegrationModelPart = &rModel.GetModelPart("integration_auxiliary_model_part");
    } else {
        mpIntegrationModelPart = &rModel.CreateModelPart("integration_auxiliary_model_part");
    }
    mDimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
}

void DepthIntegrationProcess::Execute()
{
    double bottom, top;
    GetBoundingVolumeLimits(bottom, top);
    // InitializeIntegrationModelPart();
    // FindIntersectedGeometricalObjectsProcess find_intersected_objects_process(*mpIntegrationModelPart, mrVolumeModelPart);
    // find_intersected_objects_process.ExecuteInitialize();
    // for (auto& node : mrInterfaceModelPart.Nodes()) {
    //     InitializeIntegrationLine();
    //     SetIntegrationLine(node, bottom, top);
    //     find_intersected_objects_process.FindIntersections();
    //     auto intersected_objects = find_intersected_objects_process.GetIntersections();
    //     Integrate(intersected_objects[0], node);
    // }
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
    const array_1d<double,3> momentum = height*velocity;
    SetValue(rNode, MOMENTUM, momentum);
    SetValue(rNode, VELOCITY, velocity);
    SetValue(rNode, HEIGHT, height);
}

void DepthIntegrationProcess::GetBoundingVolumeLimits(double& rMin, double& rMax)
{
    using MultipleReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

    std::tie(rMin, rMax) = block_for_each<MultipleReduction>(mrVolumeModelPart.Nodes(), [&](NodeType& node){
        const double distance = inner_prod(mDirection, node);
        return std::make_tuple(distance, distance);
    });
}

void DepthIntegrationProcess::InitializeIntegrationModelPart()
{
    // Empty the model part
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Nodes());
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Elements());
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Conditions());
    mpIntegrationModelPart->RemoveNodesFromAllLevels();
    mpIntegrationModelPart->RemoveElementsFromAllLevels();
    mpIntegrationModelPart->RemoveConditionsFromAllLevels();

    // Add a dummy node in order to allow the generation of the octree
    if (mrVolumeModelPart.NumberOfNodes() > 0) {
        auto i_node = mrVolumeModelPart.NodesBegin();
        mpIntegrationModelPart->AddNode(&*i_node);
    } else {
        KRATOS_ERROR << "DepthIntegrationProcess: The volume model part is empty." << std::endl;
    }
}

void DepthIntegrationProcess::InitializeIntegrationLine()
{
    // Remove the dummy node
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Nodes());
    mpIntegrationModelPart->RemoveNodesFromAllLevels();

    // Set the element name
    std::string element_name = (mDimension == 2) ? "Element2D2N" : "Element3D2N";

    // Create the integration lines
    ModelPart::PropertiesType::Pointer p_prop;
    if (mpIntegrationModelPart->HasProperties(1)) {
        p_prop = mpIntegrationModelPart->pGetProperties(1);
    } else {
        p_prop = mpIntegrationModelPart->CreateNewProperties(1);
    }
    mpIntegrationModelPart->CreateNewNode(1, 0.0, 0.0, 0.0);
    mpIntegrationModelPart->CreateNewNode(2, 0.0, 0.0, 1.0);
    mpIntegrationModelPart->CreateNewElement(element_name, 1, {{1, 2}}, p_prop);
}

void DepthIntegrationProcess::SetIntegrationLine(
    const NodeType& rNode,
    const double Bottom,
    const double Top)
{
    // Set the integration limits to the lines
    array_1d<double,3> base_point = rNode - mDirection * inner_prod(rNode, mDirection);
    array_1d<double,3> start = base_point + Top * mDirection;
    array_1d<double,3> end = base_point + Bottom * mDirection;
    auto& integration_line = mpIntegrationModelPart->GetElement(1);
    integration_line.GetGeometry()[0].Coordinates() = start;
    integration_line.GetGeometry()[1].Coordinates() = end;
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
    if (mDimension == 2)
        return Kratos::make_shared<Line2D2<NodeType>>(nodes);
    else
        return Kratos::make_shared<Line3D2<NodeType>>(nodes);
}

int DepthIntegrationProcess::Check()
{
    KRATOS_ERROR_IF(mDirection.size() != 3) << "DepthIntegrationProcess: The direction of integration must have three coordinates." << std::endl;
    return 0;
}

}  // namespace Kratos.
