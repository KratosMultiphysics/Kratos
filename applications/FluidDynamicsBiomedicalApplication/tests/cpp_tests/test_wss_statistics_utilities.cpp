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
    rModelPart.GetProcessInfo().SetValue(STEP, 2); // Required as WSS checks step > buffer

    // Geometry creation
    auto p_point_1 = Kratos::make_intrusive<Node<3>>(1, 0.0, 0.0, 0.0);
    auto p_point_2 = Kratos::make_intrusive<Node<3>>(2, 1.0, 0.0, 0.0);
    auto p_point_3 = Kratos::make_intrusive<Node<3>>(3, 1.0, 1.0, 0.0);
    auto p_point_4 = Kratos::make_intrusive<Node<3>>(4, 0.0, 1.0, 0.0);
    auto p_point_5 = Kratos::make_intrusive<Node<3>>(5, 0.0, 0.0, 1.0);
    auto p_point_6 = Kratos::make_intrusive<Node<3>>(6, 1.0, 0.0, 1.0);
    auto p_point_7 = Kratos::make_intrusive<Node<3>>(7, 1.0, 1.0, 1.0);
    auto p_point_8 = Kratos::make_intrusive<Node<3>>(8, 0.0, 1.0, 1.0);

    Hexahedra3D8<Node<3>> geometry(
        p_point_1, p_point_2, p_point_3, p_point_4, p_point_5, p_point_6, p_point_7, p_point_8);

    Parameters mesher_parameters(R"({
        "number_of_divisions"        : 2,
        "element_name"               : "Element3D4N",
        "condition_name"             : "SurfaceCondition",
        "create_skin_sub_model_part" : true
    })");
    StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();

    // Set the skin nodal data emulating the CFD reaction
    const auto& r_skin_mp = rModelPart.GetSubModelPart("Skin");
    for (auto& r_node : r_skin_mp.Nodes()) {
        r_node.FastGetSolutionStepValue(REACTION_X) = r_node.X();
        r_node.FastGetSolutionStepValue(REACTION_Y) = r_node.Y();
        r_node.FastGetSolutionStepValue(REACTION_Z) = r_node.Z();
    }
}

KRATOS_TEST_CASE_IN_SUITE(WSSUtilities3DWSS, FluidDynamicsBiomedicalApplicationFastSuite)
{
    // Set up the test model part
    Model model;
    auto& r_test_model_part = model.CreateModelPart("TestModelPart");
    TetrahedraModelPartForWSSTests(r_test_model_part);

    // Calculate nodal normals
    const std::size_t domain_size = r_test_model_part.GetProcessInfo()[DOMAIN_SIZE];
    NormalCalculationUtils().CalculateOnSimplexNonHistorical(r_test_model_part, domain_size, NORMAL);

    // Calculate WSS
    const bool is_normal_historical = false;
    WssStatisticsUtilities::InitializeWSSVariables(r_test_model_part);
    WssStatisticsUtilities::CalculateWSS(r_test_model_part, NORMAL, is_normal_historical);

    // Check results
    const double tolerance = 1.0e-8;
    const auto& r_node = *(r_test_model_part.GetSubModelPart("Skin").NodesEnd() - 1);
    KRATOS_WATCH(r_node)
    KRATOS_WATCH(r_node.GetValue(NORMAL))
    KRATOS_WATCH(r_node.GetValue(FACE_LOAD))
    KRATOS_WATCH(r_node.FastGetSolutionStepValue(REACTION))
    KRATOS_WATCH(r_node.GetValue(REACTION))
    std::cout << std::setprecision(12) << r_node.GetValue(WSS) << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_NORMAL_STRESS)[0] << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_NORMAL_STRESS)[1] << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_NORMAL_STRESS)[2] << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_TANGENTIAL_STRESS)[0] << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_TANGENTIAL_STRESS)[1] << std::endl;
    std::cout << std::setprecision(12) << r_node.GetValue(WSS_TANGENTIAL_STRESS)[2] << std::endl;
    KRATOS_CHECK_NEAR(r_node.GetValue(WSS), 0.0, 1e-8);
    KRATOS_CHECK_VECTOR_NEAR(r_node.GetValue(WSS_NORMAL_STRESS), 2.0, 1e-8);
    KRATOS_CHECK_VECTOR_NEAR(r_node.GetValue(WSS_TANGENTIAL_STRESS), 0.0, 1e-8);
}

}
}