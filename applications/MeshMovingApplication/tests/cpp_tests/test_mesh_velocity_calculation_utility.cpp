//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

#include <unordered_map>

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/variable_utils.h"
#include "includes/mesh_moving_variables.h"

// Application includes
#include "custom_utilities/calculate_mesh_velocity_utility.h"

namespace Kratos {
namespace Testing {

typedef std::unordered_map<std::size_t, std::vector<double>> ResultsMapType;

void CreateModelPartForMeshVelComputation(
    Model& rModel,
    const std::size_t BufferSize,
    const double DeltaTime,
    const bool RequiresMeshAcceleration)
{
    auto& r_model_part = rModel.CreateModelPart("MeshVelMP");

    Node<3>::Pointer p_point_1 = Kratos::make_shared<Node<3>>(1, 0.0, 0.0, 0.0);
    Node<3>::Pointer p_point_2 = Kratos::make_shared<Node<3>>(2, 0.0, 1.0, 0.0);
    Node<3>::Pointer p_point_3 = Kratos::make_shared<Node<3>>(3, 1.0, 1.0, 0.0);
    Node<3>::Pointer p_point_4 = Kratos::make_shared<Node<3>>(4, 1.0, 0.0, 0.0);

    Quadrilateral2D4<Node<3> > geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"(
    {
        "number_of_divisions": 3,
        "element_name": "Element2D3N"
    })");

    r_model_part.SetBufferSize(BufferSize);

    r_model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT);
    r_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    if (RequiresMeshAcceleration) {
        r_model_part.AddNodalSolutionStepVariable(MESH_ACCELERATION);
    }

    StructuredMeshGeneratorProcess(geometry, r_model_part, mesher_parameters).Execute();

    // Create a fake time loop to fill the buffer
    for (std::size_t i_step = 0; i_step < BufferSize+2; ++i_step){
        r_model_part.CloneTimeStep(i_step * DeltaTime);
        array_1d<double,3> u_val = ZeroVector(3);
        for (auto it_node : r_model_part.NodesArray()) {
            u_val(0) = 0.5*i_step * DeltaTime * it_node->X();
            u_val(1) = i_step * DeltaTime * it_node->Y();
            it_node->GetSolutionStepValue(MESH_DISPLACEMENT) = u_val;
        }
    }
}

template< class TVarType>
void CheckNodalResults(const ModelPart& rModelPart,
                       const ResultsMapType& rResultsMap,
                       const TVarType& rVariable,
                       const std::size_t StepIndex,
                       const bool SetupTest)
{
    const double tol = 1e-10;
    for (const auto& r_results : rResultsMap) {
        const std::size_t node_id = r_results.first;
        const auto& r_node = rModelPart.GetNode(node_id);
        const double result = r_node.GetSolutionStepValue(rVariable);
        const double expected_result = r_results.second[StepIndex];

        if (SetupTest) {
            std::cout << std::setprecision(15) << "Node # " << node_id << " | step: "
                      << StepIndex << " | Variable: " << rVariable.Name()
                      << " | result: " << result << std::endl;
        }
        else {
            KRATOS_CHECK_NEAR(expected_result, result, tol);
        }
    }
}

void TestBDFTimeIntegration(const std::size_t IntegrationOrder,
                            const ResultsMapType& rResultsX,
                            const ResultsMapType& rResultsY,
                            const bool SetupTest)
{
    Model current_model;
    const std::string name("bdf" + std::to_string(IntegrationOrder));
    const std::size_t buffer_size = CalculateMeshVelocityUtility::GetMinimumBufferSize(name);
    const std::size_t num_steps = 3;
    const double delta_time = 0.1;
    const bool requires_mesh_acceleration = false;

    CreateModelPartForMeshVelComputation(current_model, buffer_size, delta_time, requires_mesh_acceleration);

    auto& r_model_part = current_model.GetModelPart("MeshVelMP");

    Parameters utility_params(R"(
        {
            "time_scheme": ""
        })");
    utility_params["time_scheme"].SetString(name);

    CalculateMeshVelocityUtility utility(r_model_part, utility_params);

    // pseudo time-loop to check results
    const double current_time = r_model_part.GetProcessInfo()[TIME];
    for (std::size_t i_step = 0; i_step < num_steps; ++i_step){
        r_model_part.CloneTimeStep(current_time + (i_step+1) * delta_time);
        array_1d<double,3> u_val = ZeroVector(3);
        for (auto it_node : r_model_part.NodesArray()) {
            // apply some deformation
            u_val(0) = 2.0*std::pow(i_step, 1.82) * delta_time * it_node->X();
            u_val(1) = 1.0*std::pow(i_step, 3.951) * delta_time * it_node->Y();
            it_node->GetSolutionStepValue(MESH_DISPLACEMENT) = u_val;
        }
        utility.CalculateMeshVelocities();

        CheckNodalResults(r_model_part, rResultsX, MESH_VELOCITY_X, i_step, SetupTest);
        CheckNodalResults(r_model_part, rResultsY, MESH_VELOCITY_Y, i_step, SetupTest);
    }
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBDF1, MeshMovingApplicationFastSuite)
{
    const ResultsMapType results_x {
        {5 ,  {-0.5, 0.66666666666, 1.68720799011}},
        {13 , {-1.5, 2.0,           5.06162397033}},
    };

    const ResultsMapType results_y {
        {2 ,  {-1.0, 0.33333333333333, 4.82189918435}},
        {3 ,  {-2.0, 0.66666666666,    9.64379836869}},
    };

    TestBDFTimeIntegration(1, results_x, results_y, false);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesBDF2, MeshMovingApplicationFastSuite)
{
    const ResultsMapType results_x {
        {5 ,  {-1.08333333333, 1.3333333333333, 2.19747865183}},
        {13 , {-3.25,          4.0,             6.59243595549}},
    };

    const ResultsMapType results_y {
        {2 ,  {-2.16666666667, 1.16666666667, 7.06618210985}},
        {3 ,  {-4.33333333333, 2.33333333333, 14.1323642197}},
    };

    TestBDFTimeIntegration(2, results_x, results_y, false);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateMeshVelocitiesGeneralizedAlpha, MeshMovingApplicationFastSuite)
{
    Model current_model;
    const std::size_t buffer_size = CalculateMeshVelocityUtility::GetMinimumBufferSize("generalized_alpha");
    const std::size_t num_steps = 3;
    const double delta_time = 0.1;
    const bool requires_mesh_acceleration = true;
    const bool setup_test = false;

    CreateModelPartForMeshVelComputation(current_model, buffer_size, delta_time, requires_mesh_acceleration);

    auto& r_model_part = current_model.GetModelPart("MeshVelMP");

    Parameters utility_params(R"(
        {
            "time_scheme": "generalized_alpha",
            "alpha_m": -0.05,
            "alpha_f":  0.03
        })");

    CalculateMeshVelocityUtility utility(r_model_part, utility_params);

    const ResultsMapType results_v_x {
        {5 ,  {-0.99451303155, 2.30020830158, 1.11891669243}},
        {13 , {-2.98353909465, 6.90062490474, 3.3567500773}},
    };

    const ResultsMapType results_v_y {
        {2 ,  {-1.9890260631, 2.61139054006, 7.06529705614}},
        {3 ,  {-3.9780521262, 5.22278108012, 14.1305941123}},
    };

    const ResultsMapType results_a_x {
        {5 ,  {-17.146776406,  69.222171417,     -70.493496701}},
        {13 , {-51.4403292181, 207.666514250877, -211.480490102989}},
    };

    const ResultsMapType results_a_y {
        {2 ,  {-34.2935528121, 104.150790021846, 1.37195405457}},
        {3 ,  {-68.5871056241, 208.301580043693, 2.74390810914}},
    };

    // pseudo time-loop to check results
    const double current_time = r_model_part.GetProcessInfo()[TIME];
    for (std::size_t i_step = 0; i_step < num_steps; ++i_step){
        r_model_part.CloneTimeStep(current_time + (i_step+1) * delta_time);
        array_1d<double,3> u_val = ZeroVector(3);
        for (auto it_node : r_model_part.NodesArray()) {
            // apply some deformation
            u_val(0) = 2.0*std::pow(i_step, 1.82) * delta_time * it_node->X();
            u_val(1) = 1.0*std::pow(i_step, 3.951) * delta_time * it_node->Y();
            it_node->GetSolutionStepValue(MESH_DISPLACEMENT) = u_val;
        }
        utility.CalculateMeshVelocities();

        CheckNodalResults(r_model_part, results_v_x, MESH_VELOCITY_X, i_step, setup_test);
        CheckNodalResults(r_model_part, results_v_y, MESH_VELOCITY_Y, i_step, setup_test);
        CheckNodalResults(r_model_part, results_a_x, MESH_ACCELERATION_X, i_step, setup_test);
        CheckNodalResults(r_model_part, results_a_y, MESH_ACCELERATION_Y, i_step, setup_test);
    }
}

}
}  // namespace Kratos.
