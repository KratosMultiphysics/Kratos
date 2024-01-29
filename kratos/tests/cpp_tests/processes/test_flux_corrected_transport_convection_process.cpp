//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
// #include "includes/gid_io.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "processes/flux_corrected_transport_convection_process.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"
#include "utilities/variable_utils.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(FluxCorrectedTransportConvectionProcess2D, KratosCoreFastSuite)
{
    // Set-up a simplicial mesh
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");
    auto p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, -0.5, 0.5, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 0.5, 0.5, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 0.5, -0.5, 0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions" : 50,
        "element_name" : "Element2D3N",
        "condition_name" : "LineCondition",
        "create_skin_sub_model_part" : true
    })");
    r_model_part.SetBufferSize(3);
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
    VariableUtils().AddDof(DISTANCE, r_model_part);

    // Calculate nodal neighbours
    // TODO: Temporary solution until we skip the neighbours calculation
    FindGlobalNodalNeighboursProcess nodal_neighs_process(r_model_part);
    nodal_neighs_process.Execute();

    // FCT convection settings
    Parameters fct_parameters(R"({
        "model_part_name" : "ModelPart",
        "echo_level" : 1,
        "max_CFL" : 0.1,
        "diffusion_constant" : 1.0
    })");


    // // Gaussian hill
    // const double a = 1.0; // height
    // const double c = 0.1; // width
    // const double b_x = 0.15; // x-coordinate of the center
    // const double b_y = 0.15; // y-coordinate of the center
    // auto dist_func = [&](Node& rNode){return a*std::exp(-(std::pow(rNode.X()-b_x,2)/2/std::pow(c,2)+(std::pow(rNode.Y()-b_y,2)/2/std::pow(c,2))));};
    // auto vel_func = [&](Node& rNode, array_1d<double,3>& rVel){rVel[0] = -rNode.Y();rVel[1] = rNode.X();rVel[2] = 0.0;};

    // // "1D" Gaussian hill
    // const double a = 1.0; // height
    // const double c = 0.05; // width
    // const double b = 0.0; // x-coordinate of the center
    // auto dist_func = [&](Node& rNode){return a*std::exp(-(std::pow(rNode.X()-b,2)/2/std::pow(c,2)));};
    // auto vel_func = [&](Node& rNode, array_1d<double,3>& rVel){rVel[0] = 1.0;rVel[1] = 0.0;rVel[2] = 0.0;};

    // "1D" hill
    const double a = 1.0; // height
    const double c = 0.075; // width
    const double b = 0.0; // x-coordinate of the center
    auto dist_func = [&](Node& rNode){return std::abs(rNode.X() - b) < c ? a : 0.0;};
    auto vel_func = [&](Node& rNode, array_1d<double,3>& rVel){rVel[0] = 1.0;rVel[1] = 0.0;rVel[2] = 0.0;};

    // Set nodal values
    array_1d<double,3> aux_v;
    for (auto& r_node : r_model_part.Nodes()) {
        vel_func(r_node, aux_v);
        r_node.FastGetSolutionStepValue(DISTANCE, 0) = dist_func(r_node);
        r_node.FastGetSolutionStepValue(DISTANCE, 1) = dist_func(r_node);
        r_node.FastGetSolutionStepValue(VELOCITY, 0) = aux_v;
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = aux_v;
        r_node.FastGetSolutionStepValue(VELOCITY, 2) = aux_v;
    }

    // Set and execute the FCT convection process (time loop)
    // GidIO<> gid_io_convection("/home/rzorrilla/Desktop/FluxCorrectedTransportProcess2D", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_convection.InitializeMesh(0);
    // gid_io_convection.WriteMesh(r_model_part.GetMesh());
    // gid_io_convection.FinalizeMesh();
    // gid_io_convection.InitializeResults(0, r_model_part.GetMesh());
    // gid_io_convection.WriteNodalResults(DISTANCE, r_model_part.Nodes(), 0, 1);
    // gid_io_convection.WriteNodalResults(VELOCITY, r_model_part.Nodes(), 0, 1);

    const double dt = 5e-3;
    const double end_time = 5e-3;
    // const double end_time = 0.5; // Set this end time to do the complete test
    r_model_part.GetProcessInfo()[TIME] = 0.0;
    r_model_part.GetProcessInfo()[DELTA_TIME] = dt;
    FluxCorrectedTransportConvectionProcess<2> fct_convection_process(current_model, fct_parameters);
    while (r_model_part.GetProcessInfo()[TIME] < end_time) {
        const double new_time = r_model_part.GetProcessInfo()[TIME] + dt;
        r_model_part.CloneTimeStep(new_time);
        fct_convection_process.Execute();
        // gid_io_convection.WriteNodalResults(DISTANCE, r_model_part.Nodes(), new_time, 0);
        // gid_io_convection.WriteNodalResults(VELOCITY, r_model_part.Nodes(), new_time, 0);
    }

    // gid_io_convection.FinalizeResults();

    // Check results
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(1451).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(1502).FastGetSolutionStepValue(DISTANCE), 0.249900031105, 1.0e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(1553).FastGetSolutionStepValue(DISTANCE), 9.99688951356e-05, 1.0e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(1604).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-8);
}

KRATOS_TEST_CASE_IN_SUITE(FluxCorrectedTransportConvectionProcessZalesak, KratosCoreFastSuite)
{
    // Set-up a simplicial mesh
    Model current_model;
    auto& r_model_part = current_model.CreateModelPart("ModelPart");
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions" : 128,
        "element_name" : "Element2D3N",
        "condition_name" : "LineCondition",
        "create_skin_sub_model_part" : true
    })");
    r_model_part.SetBufferSize(2);
    r_model_part.AddNodalSolutionStepVariable(DISTANCE);
    r_model_part.AddNodalSolutionStepVariable(VELOCITY);
    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();
    VariableUtils().AddDof(DISTANCE, r_model_part);

    // Calculate nodal neighbours
    // TODO: Temporary solution until we skip the neighbours calculation
    FindGlobalNodalNeighboursProcess nodal_neighs_process(r_model_part);
    nodal_neighs_process.Execute();

    // FCT convection settings
    Parameters fct_parameters(R"({
        "model_part_name" : "ModelPart",
        "echo_level" : 0,
        "max_CFL" : 0.1,
        "echo_level" : 2,
        "time_scheme" : "RK3-TVD"
    })");

    // Define initial distance function
    auto dist_func = [&](Node& rNode){
        double dist_disk = 0.0;
        const double aux_1 = std::sqrt(std::pow(rNode.X() - 0.5,2) + std::pow(rNode.Y() - 0.75,2));
        if (aux_1 <= 0.15 && (std::abs(rNode.X() - 0.5) >= 0.025 || rNode.Y() >= 0.85)) {
            dist_disk = 1.0;
        }

        double dist_cone = 0.0;
        const double aux_2 = std::sqrt(std::pow(rNode.X() - 0.5,2) + std::pow(rNode.Y() - 0.25,2));
        if (aux_2 <= 0.15) {
            dist_cone = 1 - aux_2/0.15;
        }

        double dist_hump = 0.0;
        const double aux_3 = std::sqrt(std::pow(rNode.X() - 0.25, 2) + std::pow(rNode.Y() - 0.5, 2));
        if (aux_3 <= 0.15) {
            dist_hump = 0.25 + 0.25*std::cos(Globals::Pi * aux_3 / 0.15);
        }
        return std::max({dist_disk, dist_cone, dist_hump});
    };

    // Define velocity function
    auto vel_func = [&](Node& rNode, array_1d<double,3>& rVel){
        rVel[0] = 0.5 - rNode.Y();
        rVel[1] = rNode.X() - 0.5;
        rVel[2] = 0.0;
    };

    // Set nodal values
    array_1d<double,3> aux_v;
    for (auto& r_node : r_model_part.Nodes()) {
        vel_func(r_node, aux_v);
        r_node.FastGetSolutionStepValue(DISTANCE, 0) = dist_func(r_node);
        r_node.FastGetSolutionStepValue(DISTANCE, 1) = dist_func(r_node);
        r_node.FastGetSolutionStepValue(VELOCITY, 0) = aux_v;
        r_node.FastGetSolutionStepValue(VELOCITY, 1) = aux_v;
    }

    // Fix DISTANCE in the square boundaries
    for (auto& r_node : r_model_part.GetSubModelPart("Skin").Nodes()) {
        r_node.Fix(DISTANCE);
    }

    // Set and execute the FCT convection process (time loop)
    // GidIO<> gid_io_convection("/home/rzorrilla/Desktop/FluxCorrectedTransportProcessZalesak", GiD_PostAscii, SingleFile, WriteDeformed, WriteConditions);
    // gid_io_convection.InitializeMesh(0);
    // gid_io_convection.WriteMesh(r_model_part.GetMesh());
    // gid_io_convection.FinalizeMesh();
    // gid_io_convection.InitializeResults(0, r_model_part.GetMesh());
    // gid_io_convection.WriteNodalResults(DISTANCE, r_model_part.Nodes(), 0, 1);
    // gid_io_convection.WriteNodalResults(VELOCITY, r_model_part.Nodes(), 0, 1);

    const double dt = 1.0e-2;
    const double t_end = 5.0e-2;
    // const double t_end = 6.28; // Set this end time to do the complete test
    r_model_part.GetProcessInfo()[TIME] = 0.0;
    r_model_part.GetProcessInfo()[DELTA_TIME] = dt;
    FluxCorrectedTransportConvectionProcess<2> fct_convection_process(current_model, fct_parameters);
    while (r_model_part.GetProcessInfo()[TIME] < t_end) {
        const double new_time = r_model_part.GetProcessInfo()[TIME] + dt;
        r_model_part.CloneTimeStep(new_time);
        fct_convection_process.Execute();
        // gid_io_convection.WriteNodalResults(DISTANCE, r_model_part.Nodes(), new_time, 0);
        // gid_io_convection.WriteNodalResults(VELOCITY, r_model_part.Nodes(), new_time, 0);
    }

    // gid_io_convection.FinalizeResults();

    // Check results
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(4191).FastGetSolutionStepValue(DISTANCE), 0.49392428659, 1.0e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(8547).FastGetSolutionStepValue(DISTANCE), 0.93486103314, 1.0e-8);
    KRATOS_EXPECT_NEAR(r_model_part.GetNode(8739).FastGetSolutionStepValue(DISTANCE), 0.950609952007, 1.0e-8);
}

} // namespace Kratos::Testing.
