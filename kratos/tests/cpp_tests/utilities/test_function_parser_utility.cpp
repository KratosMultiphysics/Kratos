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
#include "includes/global_variables.h"
#include "utilities/function_parser_utility.h"

namespace Kratos {
namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(GenericFunctionUtility1, KratosCoreFastSuite)
{
    auto function0 = GenericFunctionUtility("2*x");
    KRATOS_EXPECT_TRUE(function0.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function0.UseLocalSystem());
    KRATOS_EXPECT_EQ(function0.FunctionBody(), "2*x");
    KRATOS_EXPECT_DOUBLE_EQ(function0.CallFunction(4.0,3.0,0.0,0.0), 8);

    auto function1 = GenericFunctionUtility("x**2+y**2");
    KRATOS_EXPECT_TRUE(function1.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function1.UseLocalSystem());
    KRATOS_EXPECT_EQ(function1.FunctionBody(), "x**2+y**2");
    KRATOS_EXPECT_DOUBLE_EQ(function1.CallFunction(4.0,3.0,0.0,0.0), 25);

    auto function2 = GenericFunctionUtility("3*t");
    KRATOS_EXPECT_FALSE(function2.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function2.UseLocalSystem());
    KRATOS_EXPECT_EQ(function2.FunctionBody(), "3*t");
    KRATOS_EXPECT_DOUBLE_EQ(function2.CallFunction(0.0,0.0,0.0,5.0), 15);

    auto function3 = GenericFunctionUtility("X**2+Y**2");
    KRATOS_EXPECT_TRUE(function3.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function3.UseLocalSystem());
    KRATOS_EXPECT_EQ(function3.FunctionBody(), "X**2+Y**2");
    KRATOS_EXPECT_DOUBLE_EQ(function3.CallFunction(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25);

    auto function4 = GenericFunctionUtility("(cos(x*pi)+sin(y*pi))*t");
    KRATOS_EXPECT_TRUE(function4.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function4.UseLocalSystem());
    KRATOS_EXPECT_EQ(function4.FunctionBody(), "(cos(x*pi)+sin(y*pi))*t");
    KRATOS_EXPECT_DOUBLE_EQ(function4.CallFunction(0.25,0.15,0.0,1.5), 1.5*(std::cos(0.25*Globals::Pi) + std::sin(0.15*Globals::Pi)));

    auto function5 = GenericFunctionUtility("(1.0)*(50*(cos(t)-2))");
    KRATOS_EXPECT_FALSE(function5.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function5.UseLocalSystem());
    KRATOS_EXPECT_EQ(function5.FunctionBody(), "(1.0)*(50*(cos(t)-2))");
    KRATOS_EXPECT_DOUBLE_EQ(function5.CallFunction(0.0,0.0,0.0,0.0), -50);

    auto function6 = GenericFunctionUtility("(1.0)*(50*(exp(t)-2))");
    KRATOS_EXPECT_FALSE(function6.DependsOnSpace());
    KRATOS_EXPECT_FALSE(function6.UseLocalSystem());
    KRATOS_EXPECT_EQ(function6.FunctionBody(), "(1.0)*(50*(exp(t)-2))");
    KRATOS_EXPECT_DOUBLE_EQ(function6.CallFunction(0.0,0.0,0.0,0.0), -50);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(GenericFunctionUtility("NotAValidExpression"), "Error: \nParsing error in function: NotAValidExpression\nError occurred near here :                   ^ (char [18])\nCheck your locale (e.g. if \".\" or \",\" is used as decimal point)");
}

KRATOS_TEST_CASE_IN_SUITE(GenericFunctionUtility2, KratosCoreFastSuite)
{
    Parameters parameters(R"input(
    {
        "origin" : [0,0,0],
        "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

    }
    )input");

    auto function = GenericFunctionUtility("x+2*y", parameters);
    KRATOS_EXPECT_TRUE(function.DependsOnSpace());
    KRATOS_EXPECT_TRUE(function.UseLocalSystem());
    KRATOS_EXPECT_EQ(function.FunctionBody(), "x+2*y");
    KRATOS_EXPECT_DOUBLE_EQ(function.CallFunction(4.0,3.0,0.0,0.0), 10);
    KRATOS_EXPECT_DOUBLE_EQ(function.RotateAndCallFunction(4.0,3.0,0.0,0.0), 11);
}

KRATOS_TEST_CASE_IN_SUITE(FunctionParser, KratosCoreFastSuite)
{
    auto parser_function1 = FunctionParser("x**2+y**2");
    auto function1 = parser_function1.GetFunctionSpace();
    KRATOS_EXPECT_DOUBLE_EQ(function1(4.0,3.0,0.0), 25);

    auto parser_function2 = FunctionParser("3*t");
    auto function2 = parser_function2.GetFunction();
    KRATOS_EXPECT_DOUBLE_EQ(function2(0.0,0.0,0.0,5.0), 15);

    auto parser_function3 = FunctionParser("X**2+Y**2");
    auto function3 = parser_function3.GetFunctionInitialCoordinates();
    KRATOS_EXPECT_DOUBLE_EQ(function3(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25);

    auto parser_function4 = FunctionParser("(cos(x*pi)+sin(y*pi))*t");
    auto function4 = parser_function4.GetFunction();
    KRATOS_EXPECT_DOUBLE_EQ(function4(0.25,0.15,0.0,1.5), 1.5*(std::cos(0.25*Globals::Pi) + std::sin(0.15*Globals::Pi)));
}

KRATOS_TEST_CASE_IN_SUITE(FunctionParserHexNumbers, KratosCoreFastSuite)
{
    auto fun_0 = GenericFunctionUtility("2 + 0xB0B0");
    KRATOS_EXPECT_FALSE(fun_0.DependsOnSpace());

    auto fun_1 = GenericFunctionUtility("3*y");
    KRATOS_EXPECT_TRUE(fun_1.DependsOnSpace());

    auto fun_2 = GenericFunctionUtility("Z");
    KRATOS_EXPECT_TRUE(fun_2.DependsOnSpace());

    auto fun_3 = GenericFunctionUtility("x + 0x34321");
    KRATOS_EXPECT_TRUE(fun_3.DependsOnSpace());
}

}   // namespace Testing
}  // namespace Kratos.
