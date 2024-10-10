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
            KRATOS_EXPECT_TRUE((*i_parameter).IsBool());
            KRATOS_EXPECT_EQ((*i_parameter).GetBool(), true);
            ++i_parameter;
            KRATOS_EXPECT_TRUE(i_parameter->IsDouble());
            KRATOS_EXPECT_EQ(i_parameter->GetDouble(), 2.0);
            ++i_parameter;
            KRATOS_EXPECT_EQ(i_parameter.name(), std::string("int_value"));
            KRATOS_EXPECT_TRUE(i_parameter->IsInt());
            KRATOS_EXPECT_EQ(i_parameter->GetInt(), 10);

            unsigned int size = 0;

            for(auto i = parameters.begin() ; i != parameters.end() ; i++)
                size++;

            KRATOS_EXPECT_EQ(size, 5);


        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Constructor1, KratosCoreFastSuite)
        {
            Parameters p;
            KRATOS_EXPECT_EQ(p.WriteJsonString(), "{}");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Swap1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": 1.0})"};
            Parameters p2{R"({"bar": false})"};
            p1.swap(p2);
            KRATOS_EXPECT_TRUE(p1.Has("bar"));
            KRATOS_EXPECT_TRUE(p1["bar"].GetBool() == false);
            KRATOS_EXPECT_TRUE(p2.Has("foo"));
            KRATOS_EXPECT_TRUE(p2["foo"].GetDouble() == 1.0);
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_Reset1, KratosCoreFastSuite)
        {
            Parameters p{R"({"foo": {"bar": 1.0}})"};
            p.Reset();
            KRATOS_EXPECT_EQ(p.WriteJsonString(), "{}");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_MoveConstructor1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": {"bar": 1.0}})"};
            Parameters p2{std::move(p1)};
            KRATOS_EXPECT_EQ(p1.WriteJsonString(), "{}");
            KRATOS_EXPECT_EQ(p2.WriteJsonString(), R"({"foo":{"bar":1.0}})");
        }

        KRATOS_TEST_CASE_IN_SUITE(Parameters_MoveAssignment1, KratosCoreFastSuite)
        {
            Parameters p1{R"({"foo": {"bar": 1.0}})"};
            Parameters p2{};
            p2 = std::move(p1);
            KRATOS_EXPECT_EQ(p1.WriteJsonString(), "{}");
            KRATOS_EXPECT_EQ(p2.WriteJsonString(), R"({"foo":{"bar":1.0}})");
        }
    }
}  // namespace Kratos.
