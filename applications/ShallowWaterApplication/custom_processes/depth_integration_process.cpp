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
#include "processes/find_intersected_geometrical_objects_process.h"

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
    if (rModel.HasModelPart("integration_auxiliary_model_part")) { // This is to allow multiple instances of this process
        mpIntegrationModelPart = &rModel.GetModelPart("integration_auxiliary_model_part");
    } else {
        mpIntegrationModelPart = &rModel.CreateModelPart("integration_auxiliary_model_part");
    }
}

void DepthIntegrationProcess::Execute()
{
    double bottom, top;
    InitializeIntegrationModelPart();
    GetBoundingVolumeLimits(bottom, top);
    CreateIntegrationLines(bottom, top);
    FindIntersectedGeometricalObjectsProcess find_intersected_objects_process(*mpIntegrationModelPart, mrVolumeModelPart);
    find_intersected_objects_process.ExecuteInitialize();
    find_intersected_objects_process.FindIntersections();
    auto intersected_objects = find_intersected_objects_process.GetIntersections();
    Integrate(intersected_objects);
}

void DepthIntegrationProcess::InitializeIntegrationModelPart()
{
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Nodes());
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Elements());
    VariableUtils().SetFlag(TO_ERASE, true, mpIntegrationModelPart->Conditions());
    mpIntegrationModelPart->RemoveNodesFromAllLevels();
    mpIntegrationModelPart->RemoveElementsFromAllLevels();
    mpIntegrationModelPart->RemoveConditionsFromAllLevels();
}

void DepthIntegrationProcess::GetBoundingVolumeLimits(double& rMin, double& rMax)
{
    using MultipleReduction = CombinedReduction<MinReduction<double>,MaxReduction<double>>; 

    std::tie(rMin, rMax) = block_for_each<MultipleReduction>(mrVolumeModelPart.Nodes(), [&](NodeType& node){
        const double distance = inner_prod(mDirection, node);
        return std::make_tuple(distance, distance);
    });
}

void DepthIntegrationProcess::CreateIntegrationLines(const double Low, const double High)
{
    // Set the element name
    const auto dimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
    std::string element_name = (dimension == 2) ? "Element2D2N" : "Element3D2N";

    // Get some properties
    ModelPart::PropertiesType::Pointer p_prop;
    if (mpIntegrationModelPart->HasProperties(1)) {
        p_prop = mpIntegrationModelPart->pGetProperties(1);
    } else {
        p_prop = mpIntegrationModelPart->CreateNewProperties(1);
    }

    // Initialize the ids counter
    std::size_t node_id = 0;
    std::size_t elem_id = 0;

    // Create the integration lines
    for (auto& node : mrInterfaceModelPart.Nodes()) {
        const double distance = inner_prod(node, mDirection);
        array_1d<double,3> start = node + mDirection * (Low + distance);
        array_1d<double,3> end = node + mDirection * (High - distance);
        mpIntegrationModelPart->CreateNewNode(++node_id, start[0], start[1], start[2]);
        mpIntegrationModelPart->CreateNewNode(++node_id, end[0], end[1], end[2]);
        mpIntegrationModelPart->CreateNewElement(element_name, ++elem_id, {{node_id-1, node_id}}, p_prop);
    }
}

void DepthIntegrationProcess::Integrate(std::vector<PointerVector<GeometricalObject>>& rResults)
{
    KRATOS_ERROR_IF(rResults.size() != mrInterfaceModelPart.NumberOfNodes()) << "DepthIntegrationProcess: the number of nodes in the interface and the number of integration lines mismatch.";
    IndexPartition<int>(static_cast<int>(rResults.size())).for_each([&](int i) {
        auto& objects_in_line = rResults[i];
        auto i_node = mrInterfaceModelPart.NodesBegin() + i;
        Integrate(objects_in_line, *i_node);
    });
}

void DepthIntegrationProcess::Integrate(PointerVector<GeometricalObject>& rObjects, NodeType& rNode)
{
    array_1d<double,3> velocity = ZeroVector(3);
    array_1d<double,3> momentum = ZeroVector(3);
    double height = 0.0;
    
    if (rObjects.size() > 0)
    {
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
        height = max_elevation - min_elevation;
        momentum = height*velocity;
    }
    SetValue(rNode, MOMENTUM, momentum);
    SetValue(rNode, VELOCITY, velocity);
    SetValue(rNode, HEIGHT, height);
}

int DepthIntegrationProcess::Check()
{
    const auto dimension = mrVolumeModelPart.GetProcessInfo()[DOMAIN_SIZE];
    KRATOS_ERROR_IF(dimension != 2 && dimension != 3) << Info() << ": Wrong DOMAIN_SIZE equal to " << dimension << "in model part " << mrVolumeModelPart.Name() << std::endl;
    KRATOS_ERROR_IF(mDirection.size() != 3) << Info() << ": The direction of integration must be given with three coordinates." << std::endl;
    KRATOS_ERROR_IF(mrVolumeModelPart.NumberOfNodes() == 0) << Info() << ": The volume model part is empty. Not possible to construct the octree." << std::endl;
    return 0;
}

}  // namespace Kratos.
