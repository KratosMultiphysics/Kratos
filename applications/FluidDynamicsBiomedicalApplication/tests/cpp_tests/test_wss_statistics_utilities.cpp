//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:   Ruben Zorrilla
//                  Eduardo Soudah
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"
#include "processes/structured_mesh_generator_process.h"
#include "testing/testing.h"
#include "utilities/normal_calculation_utils.h"

// Application includes
#include "fluid_dynamics_biomedical_application_variables.h"
#include "custom_utilities/wss_statistics_utilities.h"

namespace Kratos {
namespace Testing {

void TetrahedraModelPartForWSSTests(ModelPart& rModelPart)
{
    // Model part data
    rModelPart.AddNodalSolutionStepVariable(REACTION);
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 3);
    rModelPart.GetProcessInfo().SetValue(TIME, 2.0); // Required for the OSI-related indicators
    rModelPart.GetProcessInfo().SetValue(DELTA_TIME, 2.0); // Required for the OSI-related indicators

    // Geometry creation
    auto p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0);
    auto p_point_5 = Kratos::make_intrusive<Node>(5, 0.0, 0.0, 1.0);
    auto p_point_6 = Kratos::make_intrusive<Node>(6, 1.0, 0.0, 1.0);
    auto p_point_7 = Kratos::make_intrusive<Node>(7, 1.0, 1.0, 1.0);
    auto p_point_8 = Kratos::make_intrusive<Node>(8, 0.0, 1.0, 1.0);

    Hexahedra3D8<Node> geometry(
        p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

    Parameters mesher_parameters(R"({
        "number_of_divisions"        : 2,
        "element_name"               : "Element3D4N",
        "condition_name"             : "SurfaceCondition",
        "create_skin_sub_model_part" : true
    })");
    StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();
}

KRATOS_TEST_CASE_IN_SUITE(WSSStatisticsUtilitiesWSS, FluidDynamicsBiomedicalApplicationFastSuite)
{
    // Set up the test model part
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    TetrahedraModelPartForWSSTests(r_test_model_part);
    auto& r_test_skin_model_part = r_test_model_part.GetSubModelPart("Skin");

    // Calculate nodal normals
    const std::size_t domain_size = r_test_model_part.GetProcessInfo()[DOMAIN_SIZE];
    NormalCalculationUtils().CalculateOnSimplexNonHistorical(r_test_skin_model_part, domain_size, NORMAL);

    // Set the skin nodal data emulating the CFD reaction
    const auto& r_skin_mp = r_test_model_part.GetSubModelPart("Skin");
    for (auto& r_node : r_skin_mp.Nodes()) {
        r_node.FastGetSolutionStepValue(REACTION_X) = r_node.X();
        r_node.FastGetSolutionStepValue(REACTION_Y) = 2.0*r_node.Y();
        r_node.FastGetSolutionStepValue(REACTION_Z) = 4.0*r_node.Z();
    }

    // Calculate WSS
    const bool is_normal_historical = false;
    WssStatisticsUtilities::InitializeWSSVariables(r_test_skin_model_part);
    WssStatisticsUtilities::CalculateWSS(r_test_skin_model_part, NORMAL, is_normal_historical);

    // Check results
    const double tolerance = 1.0e-8;
    const auto &r_node = *(r_test_skin_model_part.NodesEnd() - 1);
    std::vector<double> expected_wss_norm_stress = {16.1658075373, 16.1658075373, 16.1658075373};
    std::vector<double> expected_wss_tang_stress = {-9.23760430703, -2.30940107676, 11.5470053838};
    KRATOS_EXPECT_NEAR(r_node.GetValue(WSS), 14.9666295471, tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(WSS_NORMAL_STRESS), expected_wss_norm_stress, tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(WSS_TANGENTIAL_STRESS), expected_wss_tang_stress, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(WSSStatisticsUtilitiesOSI, FluidDynamicsBiomedicalApplicationFastSuite)
{
    // Set up the test model part
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    TetrahedraModelPartForWSSTests(r_test_model_part);
    auto& r_test_skin_model_part = r_test_model_part.GetSubModelPart("Skin");

    // Set the skin nodal data emulating the already computed tangential WSS and old OSI
    WssStatisticsUtilities::InitializeWSSVariables(r_test_skin_model_part);
    array_1d<double,3> temporal_osi;
    array_1d<double,3> wss_tang_stress;
    for (auto& r_node : r_test_skin_model_part.Nodes()) {
        temporal_osi = r_node.Coordinates() / 2.0;
        wss_tang_stress = r_node.Coordinates() * -2.0;
        r_node.SetValue(TEMPORAL_OSI, temporal_osi);
        r_node.SetValue(WSS_TANGENTIAL_STRESS, wss_tang_stress);
    }

    // Calculate OSI
    WssStatisticsUtilities::CalculateOSI(r_test_skin_model_part);


    // Check results
    const double tolerance = 1.0e-8;
    const auto &r_node = *(r_test_skin_model_part.NodesEnd() - 1);
    std::vector<double> expected_temporal_osi = {-3.5,-3.5,-3.5};
    KRATOS_EXPECT_NEAR(r_node.GetValue(OSI), 0.0625, tolerance);
    KRATOS_EXPECT_NEAR(r_node.GetValue(RRT), 0.329914439537, tolerance);
    KRATOS_EXPECT_NEAR(r_node.GetValue(TWSS), 6.92820323028, tolerance);
    KRATOS_EXPECT_NEAR(r_node.GetValue(ECAP), 0.0180421959122, tolerance);
    KRATOS_EXPECT_NEAR(r_node.GetValue(TAWSS), 3.46410161514, tolerance);
    KRATOS_EXPECT_VECTOR_NEAR(r_node.GetValue(TEMPORAL_OSI), expected_temporal_osi, tolerance);
}

}
}