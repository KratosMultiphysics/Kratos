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
    auto function1 = GenericFunctionUtility("x**2+y**2");
    KRATOS_CHECK(function1.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function1.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function1.FunctionBody(), "x**2+y**2");
    KRATOS_CHECK_DOUBLE_EQUAL(function1.CallFunction(4.0,3.0,0.0,0.0), 25);

    auto function2 = GenericFunctionUtility("3*t");
    KRATOS_CHECK_IS_FALSE(function2.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function2.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function2.FunctionBody(), "3*t");
    KRATOS_CHECK_DOUBLE_EQUAL(function2.CallFunction(0.0,0.0,0.0,5.0), 15);

    auto function3 = GenericFunctionUtility("X**2+Y**2");
    KRATOS_CHECK(function3.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function3.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function3.FunctionBody(), "X**2+Y**2");
    KRATOS_CHECK_DOUBLE_EQUAL(function3.CallFunction(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25);

    auto function4 = GenericFunctionUtility("(cos(x*pi)+sin(y*pi))*t");
    KRATOS_CHECK(function4.DependsOnSpace());
    KRATOS_CHECK_IS_FALSE(function4.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function4.FunctionBody(), "(cos(x*pi)+sin(y*pi))*t");
    KRATOS_CHECK_DOUBLE_EQUAL(function4.CallFunction(0.25,0.15,0.0,1.5), 1.5*(std::cos(0.25*Globals::Pi) + std::sin(0.15*Globals::Pi)));

    KRATOS_CHECK_EXCEPTION_IS_THROWN(GenericFunctionUtility("A quien le importa lo que yo haga A quien le importa lo que yo diga  Yo soy asi, y asi seguire, nunca cambiare"), "Error: Parsing error in function: A quien le importa lo que yo haga A quien le importa lo que yo diga  Yo soy asi, y asi seguire, nunca cambiare");
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
    KRATOS_CHECK(function.DependsOnSpace());
    KRATOS_CHECK(function.UseLocalSystem());
    KRATOS_CHECK_STRING_EQUAL(function.FunctionBody(), "x+2*y");
    KRATOS_CHECK_DOUBLE_EQUAL(function.CallFunction(4.0,3.0,0.0,0.0), 10);
    KRATOS_CHECK_DOUBLE_EQUAL(function.RotateAndCallFunction(4.0,3.0,0.0,0.0), 11);
}

KRATOS_TEST_CASE_IN_SUITE(FunctionParser, KratosCoreFastSuite)
{
    auto parser_function1 = FunctionParser("x**2+y**2");
    auto function1 = parser_function1.GetFunctionSpace();
    KRATOS_CHECK_DOUBLE_EQUAL(function1(4.0,3.0,0.0), 25);

    auto parser_function2 = FunctionParser("3*t");
    auto function2 = parser_function2.GetFunction();
    KRATOS_CHECK_DOUBLE_EQUAL(function2(0.0,0.0,0.0,5.0), 15);

    auto parser_function3 = FunctionParser("X**2+Y**2");
    auto function3 = parser_function3.GetFunctionInitialCoordinates();
    KRATOS_CHECK_DOUBLE_EQUAL(function3(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25);

    auto parser_function4 = FunctionParser("(cos(x*pi)+sin(y*pi))*t");
    auto function4 = parser_function4.GetFunction();
    KRATOS_CHECK_DOUBLE_EQUAL(function4(0.25,0.15,0.0,1.5), 1.5*(std::cos(0.25*Globals::Pi) + std::sin(0.15*Globals::Pi)));
}

}   // namespace Testing
}  // namespace Kratos.
