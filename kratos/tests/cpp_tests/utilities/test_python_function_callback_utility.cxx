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
#include "utilities/python_function_callback_utility.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(PythonGenericFunctionUtility1, KratosCoreFastSuite)
{
    auto function1 = PythonGenericFunctionUtility("x**2+y**2");
    KRATOS_CHECK(function1.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function1.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function1.FunctionBody(), "x**2+y**2");
    KRATOS_CHECK_DOUBLE_EQUAL(function1.CallFunction(4.0,3.0,0.0,0.0), 25);

    auto function2 = PythonGenericFunctionUtility("3*t");
    KRATOS_CHECK_IS_FALSE(function2.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function2.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function2.FunctionBody(), "3*t");
    KRATOS_CHECK_DOUBLE_EQUAL(function2.CallFunction(0.0,0.0,0.0,5.0), 15);

    auto function3 = PythonGenericFunctionUtility("X**2+Y**2");
    KRATOS_CHECK(function3.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function3.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function3.FunctionBody(), "X**2+Y**2");
    KRATOS_CHECK_DOUBLE_EQUAL(function3.CallFunction(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25);

    auto function4 = PythonGenericFunctionUtility("(cos(x)+sin(y))*t");
    KRATOS_CHECK(function4.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function4.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function4.FunctionBody(), "(cos(x)+sin(y))*t");
    KRATOS_CHECK_DOUBLE_EQUAL(function4.CallFunction(0.25,0.15,0.0,1.5), 1.5*(std::cos(0.25) + std::sin(0.15)));
}

KRATOS_TEST_CASE_IN_SUITE(PythonGenericFunctionUtility2, KratosCoreFastSuite)
{
    Parameters parameters(R"input(
    {
        "origin" : [0,0,0],
        "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

    }
    )input");

    auto function = PythonGenericFunctionUtility("x+2*y", parameters);
    KRATOS_CHECK(function.DependsOnSpace());
    KRATOS_CHECK(function.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function.FunctionBody(), "x+2*y");
    KRATOS_CHECK_DOUBLE_EQUAL(function.CallFunction(4.0,3.0,0.0,0.0), 10);
    KRATOS_CHECK_DOUBLE_EQUAL(function.RotateAndCallFunction(4.0,3.0,0.0,0.0), 11);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyFunctionToNodesUtility, KratosCoreFastSuite)
{
    Parameters parameters(R"input(
    {
        "origin" : [0,0,0],
        "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

    }
    )input");

    auto p_function = Kratos::make_shared<PythonGenericFunctionUtility>("x+2*y", parameters);
    KRATOS_CHECK(p_function->DependsOnSpace());
    KRATOS_CHECK(p_function->UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(p_function->FunctionBody(), "x+2*y");
    KRATOS_CHECK_DOUBLE_EQUAL(p_function->CallFunction(4.0,3.0,0.0,0.0), 10);
    KRATOS_CHECK_DOUBLE_EQUAL(p_function->RotateAndCallFunction(4.0,3.0,0.0,0.0), 11);

    Model this_model;
    ModelPart& r_model_part = this_model.CreateModelPart("Main", 2);
    ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
    r_current_process_info[DOMAIN_SIZE] = 2;

    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N", false);

    auto utility = ApplyFunctionToNodesUtility(r_model_part.Nodes(), p_function);
    utility.ApplyFunction<Variable<double>>(TEMPERATURE, 1.0);

    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_CHECK_DOUBLE_EQUAL(r_node.FastGetSolutionStepValue(TEMPERATURE) - (r_node.Y() + 2.0 * r_node.X()), 0.0);
    }
}

}   // namespace Testing
}  // namespace Kratos.
