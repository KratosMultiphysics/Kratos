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
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "custom_utilities/shallow_water_utilities.h"
#include "shallow_water_application_variables.h"

namespace Kratos {

namespace Testing {

typedef ModelPart::NodeType     NodeType;

KRATOS_TEST_CASE_IN_SUITE(ShallowWaterUtilitiesComputeVelocity, ShallowWaterApplicationFastSuite)
{
    Model model;

    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions":  8,
        "element_name":         "Element2D3N",
        "condition_name":       "LineCondition"
    })");

    ModelPart& r_model_part = model.CreateModelPart("model_part");
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    r_model_part.AddNodalSolutionStepVariable(MOMENTUM);
    r_model_part.AddNodalSolutionStepVariable(HEIGHT);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    auto f = [] (double x) {return 100*std::pow(x,2) * std::pow(1-x,2);};
    auto df = [] (double x) {return 200*x * std::pow(1-x,2) + 100*std::pow(x,2)*(2*x-2);};
    block_for_each(r_model_part.Nodes(), [&](NodeType& r_node){
        array_1d<double,3> q;
        q[0] = f(r_node.X()) * df(r_node.Y());
        q[1] = -df(r_node.X()) * f(r_node.Y());
        q[2] = 0.0;
        double h = df(r_node.X());
        r_node.FastGetSolutionStepValue(MOMENTUM) = q;
        r_node.FastGetSolutionStepValue(HEIGHT) = h;
    });

    bool perform_projection = true;
    ShallowWaterUtilities().ComputeVelocity(r_model_part, perform_projection);

    const double tolerance = 1e-5;
    std::vector<double> velocity1 = {0.672913, -0.672913, 0.0};
    std::vector<double> velocity10 = {0.791841, -0.457474, 0.0};
    std::vector<double> velocity20 = {3.21839, -1.39479, 0.0};
    std::vector<double> velocity30 = {11.1471, -3.59187, 0.0};
    std::vector<double> velocity40 = {7.83417, -1.68915, 0.0};
    std::vector<double> velocity50 = {0.0, 0.0, 0.0};
    std::vector<double> velocity60 = {0.0, 0.0, 0.0};
    std::vector<double> velocity70 = {0.0, 0.0, 0.0};
    std::vector<double> velocity80 = {0.0, 0.0, 0.0};

    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(1).FastGetSolutionStepValue(VELOCITY), velocity1, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(10).FastGetSolutionStepValue(VELOCITY), velocity10, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(20).FastGetSolutionStepValue(VELOCITY), velocity20, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(30).FastGetSolutionStepValue(VELOCITY), velocity30, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(40).FastGetSolutionStepValue(VELOCITY), velocity40, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(50).FastGetSolutionStepValue(VELOCITY), velocity50, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(60).FastGetSolutionStepValue(VELOCITY), velocity60, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(70).FastGetSolutionStepValue(VELOCITY), velocity70, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(
        r_model_part.GetNode(80).FastGetSolutionStepValue(VELOCITY), velocity80, tolerance)
}


KRATOS_TEST_CASE_IN_SUITE(ShallowWaterUtilitiesHydrostaticForce, ShallowWaterApplicationFastSuite)
{
    Model model;

    Node::Pointer p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    Node::Pointer p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    Node::Pointer p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    Node::Pointer p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions"        : 8,
        "element_name"               : "WaveElement2D3N",
        "condition_name"             : "WaveCondition",
        "elements_properties_id"     : 0,
        "conditions_properties_id"   : 0,
        "create_skin_sub_model_part" : true
    })");

    ModelPart& r_model_part = model.CreateModelPart("model_part");
    r_model_part.AddNodalSolutionStepVariable(HEIGHT);

    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    auto& r_process_info = r_model_part.GetProcessInfo();
    r_process_info[GRAVITY_Z] = 9.81;

    auto& r_properties = r_model_part.GetProperties(0);
    r_properties[DENSITY] = 1000;

    auto h = [] (const NodeType& rNode) {return 1.0 + 0.1*rNode.X() + 0.2*rNode.Y();};
    block_for_each(r_model_part.Nodes(), [&](NodeType& rNode){
        rNode.FastGetSolutionStepValue(HEIGHT) = h(rNode);
    });

    const auto vertical_force = ShallowWaterUtilities().ComputeHydrostaticForces(
        r_model_part.Elements(), r_process_info);
    const auto horizontal_force = ShallowWaterUtilities().ComputeHydrostaticForces(
        r_model_part.Conditions(), r_process_info);

    const std::vector<double> ref_vertical_force = {0.0, 0.0, -11281.5};
    const std::vector<double> ref_horizontal_force = {1128.15, 2256.3, 0.0};

    const double tolerance = 1e-10;
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(vertical_force, ref_vertical_force, tolerance)
    KRATOS_CHECK_VECTOR_RELATIVE_NEAR(horizontal_force, ref_horizontal_force, tolerance)
}

} // namespace Testing

} // namespace Kratos
