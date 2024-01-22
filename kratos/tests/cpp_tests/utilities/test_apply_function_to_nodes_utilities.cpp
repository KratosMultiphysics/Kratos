//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/apply_function_to_nodes_utility.h"
#include "tests/test_utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(ApplyFunctionToNodesUtility, KratosCoreFastSuite)
{
    Parameters parameters(R"input(
    {
        "origin" : [0,0,0],
        "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

    }
    )input");

    auto p_function = Kratos::make_shared<GenericFunctionUtility>("x+2*y", parameters);
    KRATOS_EXPECT_TRUE(p_function->DependsOnSpace());
    KRATOS_EXPECT_TRUE(p_function->UseLocalSystem());
    KRATOS_EXPECT_EQ(p_function->FunctionBody(), "x+2*y");
    KRATOS_EXPECT_DOUBLE_EQ(p_function->CallFunction(4.0,3.0,0.0,0.0), 10);
    KRATOS_EXPECT_DOUBLE_EQ(p_function->RotateAndCallFunction(4.0,3.0,0.0,0.0), 11);

    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 2;

    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N", false);

    auto utility = ApplyFunctionToNodesUtility(r_model_part.Nodes(), p_function);
    utility.ApplyFunction(TEMPERATURE, 1.0);

    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(TEMPERATURE) - (r_node.Y() + 2.0 * r_node.X()), 0.0);
    }
}

}   // namespace Testing
}  // namespace Kratos.

