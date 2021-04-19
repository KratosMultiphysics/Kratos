//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt

// Project includes
#include "containers/model.h"
#include "utilities/force_and_torque_utils.h"
#include "testing/testing.h"
#include "processes/structured_mesh_generator_process.h"
#include "geometries/quadrilateral_2d_4.h"

namespace Kratos
{
namespace Testing
{


void FillModelPart(ModelPart& rModelPart, bool generate_moments=true)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT_X);
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT_Y);
    rModelPart.AddNodalSolutionStepVariable(REACTION);

    if (generate_moments) {
        rModelPart.AddNodalSolutionStepVariable(ROTATION);
        rModelPart.AddNodalSolutionStepVariable(MOMENT);
    }

    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

    // Generate mesh
    using GeometryType = Quadrilateral2D4<Node<3>>;

    GeometryType::PointsArrayType domain_corners;
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(4, 0.0, 2.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(3, 2.0, 2.0, 0.0));
    domain_corners.push_back(Kratos::make_intrusive<Node<3>>(2, 2.0, 0.0, 0.0));

    Parameters mesh_generator_parameters(R"(
    {
        "number_of_divisions" : 2,
        "element_name" : "Element2D3N",
        "create_skin_sub_model_part" : false
    })");

    StructuredMeshGeneratorProcess(
        GeometryType(domain_corners),
        rModelPart,
        mesh_generator_parameters
    ).Execute();

    // Add DoFs
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        if (generate_moments) {
            r_node.AddDof(ROTATION_Z);
        }
    }

    // Set reactions and moments
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(REACTION_X) = r_node.X() * r_node.Id();
        r_node.FastGetSolutionStepValue(REACTION_Y) = r_node.Y() * r_node.Id();
        if (generate_moments) {
            r_node.FastGetSolutionStepValue(MOMENT_Z) = r_node.Id();
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(SumForceAndTorque, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    FillModelPart(r_model_part);

    const auto force_and_moment = ForceAndTorqueUtils::SumForceAndTorque(
        r_model_part,
        array_1d<double,3> {1.0, 1.0, 0.0}
    );

    const auto& force = force_and_moment[0];
    const auto& moment = force_and_moment[1];

    KRATOS_CHECK_NEAR(force[0], 63.0, 1e-15);
    KRATOS_CHECK_NEAR(force[1], 51.0, 1e-15);
    KRATOS_CHECK_NEAR(force[2], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[0], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[1], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[2], 57.0, 1e-15);
}

KRATOS_TEST_CASE_IN_SUITE(SumForceAndTorqueNoMoment, KratosCoreFastSuite)
{
    Model model;
    ModelPart& r_model_part = model.CreateModelPart("Main", 2);
    FillModelPart(r_model_part, false);

    const auto force_and_moment = ForceAndTorqueUtils::SumForceAndTorque(
        r_model_part,
        array_1d<double,3> {1.0, 1.0, 0.0}
    );

    const auto& force = force_and_moment[0];
    const auto& moment = force_and_moment[1];

    KRATOS_CHECK_NEAR(force[0], 63.0, 1e-15);
    KRATOS_CHECK_NEAR(force[1], 51.0, 1e-15);
    KRATOS_CHECK_NEAR(force[2], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[0], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[1], 0.00, 1e-15);
    KRATOS_CHECK_NEAR(moment[2], 12.0, 1e-15);
}

} /* namespace Testing */
} /* namespace Kratos */