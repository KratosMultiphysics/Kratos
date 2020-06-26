//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes


// External includes


// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"

namespace Kratos {
    namespace Testing {

        KRATOS_TEST_CASE_IN_SUITE(ParametersIterator, KratosCoreFastSuite)
        {
            Parameters parameters(R"input(
            {
               "int_value" : 10,   "double_value": 2.0,   "bool_value" : true,   "string_value" : "hello",
               "level1":
               {
                 "list_value":[ 3, "hi", false],
                 "tmp" : 5.0
               }
            }
            )input");

            auto i_parameter = parameters.begin();
            KRATOS_CHECK((*i_parameter).IsBool());
            KRATOS_CHECK_EQUAL((*i_parameter).GetBool(), true);
            ++i_parameter;
            KRATOS_CHECK(i_parameter->IsDouble());
            KRATOS_CHECK_EQUAL(i_parameter->GetDouble(), 2.0);
            ++i_parameter;
            KRATOS_CHECK_EQUAL(i_parameter.name(), std::string("int_value"));
            KRATOS_CHECK(i_parameter->IsInt());
            KRATOS_CHECK_EQUAL(i_parameter->GetInt(), 10);

            unsigned int size = 0;

            for(auto i = parameters.begin() ; i != parameters.end() ; i++)
                size++;

            KRATOS_CHECK_EQUAL(size, 5);


        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Constructor1, KratosCoreFastSuite)
        {
            Parameters p;
            KRATOS_CHECK_EQUAL(p.WriteJsonString(), "{}");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Swap1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": 1.0})"};
            Parameters p2{R"({"bar": false})"};
            p1.swap(p2);
            KRATOS_CHECK(p1.Has("bar"));
            KRATOS_CHECK(p1["bar"].GetBool() == false);
            KRATOS_CHECK(p2.Has("foo"));
            KRATOS_CHECK(p2["foo"].GetDouble() == 1.0);
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Reset1, KratosCoreFastSuite)
        {
            Parameters p{R"({"foo": {"bar": 1.0}})"};
            p.Reset();
            KRATOS_CHECK_EQUAL(p.WriteJsonString(), "{}");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_MoveConstructor1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": {"bar": 1.0}})"};
            Parameters p2{std::move(p1)};
            KRATOS_CHECK_EQUAL(p1.WriteJsonString(), "{}");
            KRATOS_CHECK_EQUAL(p2.WriteJsonString(), R"({"foo":{"bar":1.0}})");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_MoveAssignment1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": {"bar": 1.0}})"};
            Parameters p2{};
            p2 = std::move(p1);
            KRATOS_CHECK_EQUAL(p1.WriteJsonString(), "{}");
            KRATOS_CHECK_EQUAL(p2.WriteJsonString(), R"({"foo":{"bar":1.0}})");
        }
    }
}  // namespace Kratos.
