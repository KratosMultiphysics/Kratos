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
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/parallel_utilities.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "custom_processes/depth_integration_process.h"
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::NodeType NodeType;

void FillModelParts2D(ModelPart& rVolumeModelPart, ModelPart& rInterfaceModelPart)
{
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions"        : 5,
        "element_name"               : "Element2D3N",
        "create_skin_sub_model_part" : false
    })");

    auto& root_model_part = rVolumeModelPart.GetRootModelPart();
    root_model_part.AddNodalSolutionStepVariable(VELOCITY);
    root_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    root_model_part.AddNodalSolutionStepVariable(HEIGHT);
    root_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    root_model_part.GetProcessInfo().SetValue(GRAVITY, array_1d<double,3>({0.0, -9.81, 0.0}));
    StructuredMeshGeneratorProcess(geometry, rVolumeModelPart, mesher_parameters).Execute();

    auto id = root_model_part.NumberOfNodes();
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 0.0, 0.0);
}

void FillModelParts3D(ModelPart& rVolumeModelPart, ModelPart& rInterfaceModelPart)
{
    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0);
    Node::Pointer p_point_5 = Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0);
    Node::Pointer p_point_6 = Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0);
    Node::Pointer p_point_7 = Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0);
    Node::Pointer p_point_8 = Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0);

    Hexahedra3D8<Node> geometry(
        p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions"        : 3,
        "element_name"               : "Element3D4N",
        "create_skin_sub_model_part" : false
    })");

    auto& root_model_part = rVolumeModelPart.GetRootModelPart();
    root_model_part.AddNodalSolutionStepVariable(VELOCITY);
    root_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    root_model_part.AddNodalSolutionStepVariable(HEIGHT);
    root_model_part.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    root_model_part.GetProcessInfo().SetValue(GRAVITY, array_1d<double,3>({0.0, 0.0, -9.81}));
    StructuredMeshGeneratorProcess(geometry, rVolumeModelPart, mesher_parameters).Execute();

    auto id = root_model_part.NumberOfNodes();
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 0.00, 0.0);
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 0.25, 0.0);
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 0.50, 0.0);
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 0.75, 0.0);
    rInterfaceModelPart.CreateNewNode(++id, 0.5, 1.00, 0.0);
}

void ApplyVelocityField(ModelPart& rModelPart)
{
    block_for_each(rModelPart.Nodes(), [](NodeType& rNode){
        array_1d<double,3>& velocity = rNode.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3>& coords = rNode.Coordinates();
        velocity[0] = coords[1] + coords[2] * coords[2];
        velocity[1] = 0.0;
        velocity[2] = 0.0;
    });
}

KRATOS_TEST_CASE_IN_SUITE(DepthIntegrationProcess2D, ShallowWaterApplicationFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("model_part");
    auto& r_volume_model_part = r_model_part.CreateSubModelPart("volume");
    auto& r_interface_model_part = r_model_part.CreateSubModelPart("interface");
    FillModelParts2D(r_volume_model_part, r_interface_model_part);
    ApplyVelocityField(r_volume_model_part);

    Parameters process_parameters(R"(
    {
        "volume_model_part_name"    : "model_part.volume",
        "interface_model_part_name" : "model_part.interface",
        "store_historical_database" : false
    })");
    DepthIntegrationProcess<2>(model, process_parameters).Execute();

    std::vector<std::vector<double>> reference;
    reference.push_back({0.55, 0.0, 0.0});

    for(std::size_t i = 0; i < r_interface_model_part.NumberOfNodes(); ++i) {
        auto i_node = r_interface_model_part.NodesBegin() + i;
        KRATOS_CHECK_VECTOR_NEAR(i_node->GetValue(VELOCITY), reference[i], 1e-6);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DepthIntegrationProcess3D, ShallowWaterApplicationFastSuite)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("model_part");
    auto& r_volume_model_part = r_model_part.CreateSubModelPart("volume");
    auto& r_interface_model_part = r_model_part.CreateSubModelPart("interface");
    FillModelParts3D(r_volume_model_part, r_interface_model_part);
    ApplyVelocityField(r_volume_model_part);

    Parameters process_parameters(R"(
    {
        "volume_model_part_name"    : "model_part.volume",
        "interface_model_part_name" : "model_part.interface",
        "store_historical_database" : false
    })");
    DepthIntegrationProcess<3>(model, process_parameters).Execute();

    std::vector<std::vector<double>> reference;
    reference.push_back({0.426304, 0.0, 0.0});
    reference.push_back({0.657407, 0.0, 0.0});
    reference.push_back({0.935185, 0.0, 0.0});
    reference.push_back({1.157407, 0.0, 0.0});
    reference.push_back({1.435185, 0.0, 0.0});

    for (std::size_t i = 0; i < r_interface_model_part.NumberOfNodes(); ++i) {
        auto i_node = r_interface_model_part.NodesBegin() + i;
        KRATOS_CHECK_VECTOR_NEAR(i_node->GetValue(VELOCITY), reference[i], 1e-6);
    }
}

} // namespace Testing

} // namespace Kratos
