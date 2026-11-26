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
//                   Riccardo Rossi
//

// System includes
#include <filesystem>

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"
#include "tests/test_utilities/scoped_file.h"

namespace Kratos {
namespace Testing {

// input string with ugly formatting
std::string GetJSONString()
{
    return R"(
    {
      "bool_value" : true, "double_value": 2.0, "int_value" : 10,
      "level1":
      {
        "list_value":[ 3, "hi", false],
        "tmp" : 5.0
      },
      "string_value" : "hello"
    })";
}

std::string GetJSONStringPrettyOut()
{
    return R"({
    "bool_value": true,
    "double_value": 2.0,
    "int_value": 10,
    "level1": {
        "list_value": [
            3,
            "hi",
            false
        ],
        "tmp": 5.0
    },
    "string_value": "hello"
})";
}

std::string GetJSONStringPrettyOutAfterChange()
{
    return R"({
    "bool_value": true,
    "double_value": 2.0,
    "int_value": 10,
    "level1": {
        "list_value": [
            "changed",
            "hi",
            false
        ],
        "tmp": 5.0
    },
    "string_value": "hello"
})";
}

// here the level1 var is set to a double so that a validation error should be thrown
std::string GetJSONStringWrongType()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": 0.0,
        "string_value": "hello"
    })";
}

// int value is badly spelt
std::string GetJSONStringWrongSpelling()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_values": 10,
        "level1": 0.0,
        "string_value": "hello"
    })";
}

// wrong on the first level
// error shall be only detective by recursive validation
std::string GetJSONStringWrongLevel2()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": { "a":0.0 },
        "string_value": "hello"
    })";
}

std::string GetJSONStringDefaults()
{
    return R"({
        "bool_value": false,
        "double_value": 2.0,
        "int_value": 10,
        "level1": {
            "list_value": [
                3,
                "hi",
                false
            ],
            "tmp": "here we expect a string"
        },
        "new_default_obj": {
            "aaa": "string",
            "bbb": false,
            "ccc": 22
        },
        "new_default_value": -123.0,
        "string_value": "hello"
    })";
}

std::string GetJSONStringIncomplete()
{
    return R"({
        "level1": {
        },
        "new_default_obj": {
            "aaa": "string",
            "bbb": false,
            "ccc": 22
        },
        "new_default_value": -123.0,
        "string_value": "hello"
    })";
}

std::string GetJSONStringIncompleteWithExtraParameter()
{
    return R"({
        "level1": {
            "new_sublevel": "this should only be assigned in recursive"
        },
        "new_default_obj": {
            "aaa": "string",
            "bbb": false,
            "ccc": 22
        },
        "new_default_value": -123.0,
        "string_value": "hello"
    })";
}

std::string GetJSONStringExpectedValidationOutput()
{
    return R"({
    "bool_value": true,
    "double_value": 2.0,
    "int_value": 10,
    "level1": {
        "list_value": [
            3,
            "hi",
            false
        ],
        "tmp": 5.0
    },
    "new_default_obj": {
        "aaa": "string",
        "bbb": false,
        "ccc": 22
    },
    "new_default_value": -123.0,
    "string_value": "hello"
})";
}

std::string GetJSONStringFourLevels()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": {
            "level2": {
                "level3": {
                    "level4": {
                    }
                }
            }
        },
        "string_value": "hello"
    })";
}

std::string GetJSONStringForLevelsVariation()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": {
            "a":11.0,
            "level2": {
                "level3": {
                    "level4": {
                    }
                }
            }
        },
        "string_value": "hello"
    })";
}

std::string GetJSONStringForLevelsWrongVariation()
{
    return R"({
        "int_value": 10,
        "double_value": "hi",
        "bool_value": true,
        "string_value": "hello",
        "level1": {
            "a":11.0,
            "level2": {
                "level3": {
                    "level4": {
                    }
                }
            }
        }
    })";
}

std::string GetJSONStringForLevelsDefaults()
{
    return R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": {
            "a":1.0,
            "level2": {
                "b":2.0,
                "level3": {
                    "c":3.0,
                    "level4": {
                        "d":4.0
                    }
                }
            }
        },
        "string_value": "hello"
    })";
}

std::string GetJSONStringWithIncludes()
{
    return R"(
    {
      "bool_value" : true, "double_value": 2.0, "int_value" : 10,
      "@include_json" : "test_included_parameters.json"
    })";
}

std::string GetIncludedJSONString()
{
    return R"({
        "level1":
        {
        "list_value":[ 3, "hi", false],
        "tmp" : 5.0
        },
        "@include_json" : "test_included_parameters_level2.json"
    })";
}

std::string GetIncludedJSONLevel2String()
{
    return R"({"string_value":"hello"})";
}

std::string GetCircularIncludeJSONString(int FileIndex, int IncludeIndex)
{
    std::stringstream stream;
    stream << R"({"@include_json":"test_cyclic_)" << FileIndex << "_" << IncludeIndex << R"(.json"})"; // could be nicer with fmtlib or C++20
    return stream.str();
}

Parameters MakeParameters()
{
    return Parameters(R"({
        "int_value" : 10,               /* Block comments. */
        "double_value": 2.0,            // Line comments.
        "bool_value" : true,            /* Nihilistic comment for historical reasons */
        "string_value" : "hello",
        "string_array" : ["hello", "world"],
        "vector_value" : [5.2,-3.1,4.33],
        "matrix_value" : [[1,2],[3,4],[5,6]]
    })");
}

Parameters MakeParametersVectorInterface()
{
    return Parameters(R"({
        "valid_vectors" : [
            []
        ],
        "false_vectors" : [ 
            [[]],
            [[2,3],2],
            [2,3,[2]],
            [2,3,[]],
            [{"key":3},2],
            [2,3,{"key":3}],
            [true,2],
            [2,3,true],
            [5,"string",2]
        ]
    })");
}

Parameters MakeMatrixShapeParameters()
{
    return Parameters(R"({
        "valid_matrices" : [ 
            [[]],
            [[],[]],
            [[-9.81,8, 5.47]]
        ],
        "false_matrices" : [ 
            [],
            [[[]]],
            [[3.3] , [1,2]],
            [[2,1.5,3.3] , [3,{"key":3},2]],
            [[2,1.5,3.3] , [5,false,2]],
            [[2,1.5,3.3] , [[2,3],1,2]],
            [[2,1.5,3.3] , ["string",2,9]]
        ]
    })");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParameters, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONString());
    KRATOS_EXPECT_EQ(
        kp.WriteJsonString(),
        R"({"bool_value":true,"double_value":2.0,"int_value":10,"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"})"
    );

    KRATOS_EXPECT_TRUE(kp.Has("int_value"));
    KRATOS_EXPECT_FALSE(kp.Has("unextisting_value"));

    KRATOS_EXPECT_EQ(kp["int_value"].GetInt(), 10);
    KRATOS_EXPECT_EQ(kp["double_value"].GetDouble(), 2.0);
    KRATOS_EXPECT_EQ(kp["bool_value"].GetBool(), true);
    KRATOS_EXPECT_EQ(kp["string_value"].GetString(), "hello");

    KRATOS_EXPECT_EQ(kp.PrettyPrintJsonString(), GetJSONStringPrettyOut());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersChangeParameters, KratosCoreFastSuite)
{
    // Now change one item in the sublist
    Parameters kp = Parameters(GetJSONString());
    Parameters subparams = kp["level1"];

    Parameters my_list = subparams["list_value"];

    for (auto& r_param : my_list) {
        if (r_param.IsBool()) {
            KRATOS_EXPECT_FALSE(r_param.GetBool())
        }
    }

    // my_list = subparams["list_value"]
    subparams["list_value"][0].SetString("changed");

    KRATOS_EXPECT_EQ(
        kp.PrettyPrintJsonString(),
        GetJSONStringPrettyOutAfterChange()
    );
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersCopy, KratosCoreFastSuite)
{
    // Try to make a copy
    Parameters kp = Parameters(GetJSONString());
    auto original_out = kp.PrettyPrintJsonString();
    auto other_copy = kp.Clone();

    KRATOS_EXPECT_EQ(
        other_copy.PrettyPrintJsonString(),
        original_out
    );

    other_copy["int_value"].SetInt(-1);
    KRATOS_EXPECT_EQ(kp["int_value"].GetInt(), 10);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWrongParameters, KratosCoreFastSuite)
{
    // Should check which errors are thrown!!
    Parameters kp = Parameters(GetJSONString());
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(kp["no_value"].GetInt(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsDueToWrongTypes, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongType());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(kp.ValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsDueToWrongSpelling, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongSpelling());
    Parameters  defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(kp.ValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsErrorsOnFirstLevel, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongLevel2());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(kp.RecursivelyValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursiveValidation4Levels, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringFourLevels());
    Parameters kp_variation = Parameters(GetJSONStringForLevelsVariation());
    Parameters kp_wrong_wariation = Parameters(GetJSONStringForLevelsWrongVariation());
    Parameters defaults_params = Parameters(GetJSONStringForLevelsDefaults());

    kp.RecursivelyValidateAndAssignDefaults(defaults_params);
    kp_variation.RecursivelyValidateAndAssignDefaults(defaults_params);

    KRATOS_EXPECT_TRUE( kp.IsEquivalentTo(defaults_params) );
    KRATOS_EXPECT_FALSE( kp_variation.IsEquivalentTo(defaults_params) );

    KRATOS_EXPECT_TRUE( kp.HasSameKeysAndTypeOfValuesAs(defaults_params) );
    KRATOS_EXPECT_TRUE( kp_variation.HasSameKeysAndTypeOfValuesAs(defaults_params) );
    KRATOS_EXPECT_FALSE( kp_wrong_wariation.HasSameKeysAndTypeOfValuesAs(defaults_params) );
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationSuccedsErroronFirstLevel, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongLevel2());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());

    // Here no error shall be thrown since validation is only done on level0
    kp.ValidateAndAssignDefaults(defaults_params);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationSucceeds, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONString());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());
    defaults_params["level1"]["tmp"].SetDouble(2.0);  // this does not coincide with the value in kp, but is of the same type

    kp.ValidateAndAssignDefaults(defaults_params);
    KRATOS_EXPECT_EQ(kp.PrettyPrintJsonString(), GetJSONStringExpectedValidationOutput());

    KRATOS_EXPECT_DOUBLE_EQ(kp["level1"]["tmp"].GetDouble(), 5.0);  // not 2, since kp overwrites the defaults
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddMissingParameters, KratosCoreFastSuite)
{
    // Only missing parameters are added, no complaints if there already exist more than in the defaults
    Parameters kp = Parameters(GetJSONString());
    Parameters tmp = Parameters(GetJSONStringIncompleteWithExtraParameter());

    kp.AddMissingParameters(tmp);

    KRATOS_EXPECT_EQ(kp["new_default_obj"]["aaa"].GetString(), "string");
    KRATOS_EXPECT_EQ(kp["string_value"].GetString(), "hello");
    KRATOS_EXPECT_FALSE(kp["level1"].Has("new_sublevel"));
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursivelyAddMissingParameters, KratosCoreFastSuite)
{
    // Only missing parameters are added, no complaints if there already exist more than in the defaults
    Parameters kp = Parameters(GetJSONString());
    Parameters tmp = Parameters(GetJSONStringIncompleteWithExtraParameter());

    kp.RecursivelyAddMissingParameters(tmp);

    KRATOS_EXPECT_TRUE(kp["level1"].Has("new_sublevel"));
    KRATOS_EXPECT_EQ(kp["level1"]["new_sublevel"].GetString(), "this should only be assigned in recursive");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidateDefaults, KratosCoreFastSuite)
{
    // Only parameters from defaults are validated, no new values are added
    Parameters kp = Parameters(GetJSONStringIncompleteWithExtraParameter());
    Parameters tmp = Parameters(GetJSONStringDefaults());

    kp.ValidateDefaults(tmp);

    KRATOS_EXPECT_FALSE(kp.Has("bool_value"));
    KRATOS_EXPECT_FALSE(kp.Has("double_value"));
    KRATOS_EXPECT_TRUE(kp.Has("level1"));
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursivelyValidateDefaults, KratosCoreFastSuite)
{
    // Only parameters from defaults are validated, no new values are added
    Parameters kp = Parameters(GetJSONStringIncomplete());
    Parameters tmp = Parameters(GetJSONStringDefaults());

    kp.RecursivelyValidateDefaults(tmp);

    KRATOS_EXPECT_FALSE(kp.Has("bool_value"));
    KRATOS_EXPECT_FALSE(kp.Has("double_value"));
    KRATOS_EXPECT_TRUE(kp.Has("level1"));
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursivelyValidateDefaultsFail, KratosCoreFastSuite)
{
    // only parameters from defaults are validated, no new values are added
    Parameters kp = Parameters(GetJSONStringIncompleteWithExtraParameter());
    Parameters tmp = Parameters(GetJSONStringDefaults());

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(kp.RecursivelyValidateDefaults(tmp), "");

    // Sub_level
    KRATOS_EXPECT_FALSE(kp["level1"].Has("tmp"));
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddValue, KratosCoreFastSuite)
{
    Parameters kp = Parameters(R"({})");
    kp.AddEmptyValue("new_double").SetDouble(1.0);

    KRATOS_EXPECT_TRUE(kp.Has("new_double"));
    KRATOS_EXPECT_EQ(kp["new_double"].GetDouble(), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddEmptyArray, KratosCoreFastSuite)
{
    Parameters kp = Parameters(R"({})");
    kp.AddEmptyArray("new_array");

    KRATOS_EXPECT_TRUE(kp.Has("new_array"));
    KRATOS_EXPECT_EQ(kp["new_array"].size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRemoveValue, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONString());
    KRATOS_EXPECT_TRUE(kp.Has("int_value"));
    KRATOS_EXPECT_TRUE(kp.Has("level1"));

    kp.RemoveValue("int_value");
    kp.RemoveValue("level1");

    KRATOS_EXPECT_FALSE(kp.Has("int_value"));
    KRATOS_EXPECT_FALSE(kp.Has("level1"));
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIterators, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONString());

    // Iteration by range
    int nitems = 0;
    for(auto it=kp.begin(); it!=kp.end(); ++it) {
        ++nitems;
    }
    KRATOS_EXPECT_EQ(nitems, 5);

    // Iteration by items
    for(auto it=kp.begin(); it!=kp.end(); ++it) {
        KRATOS_EXPECT_EQ(kp[it.name()].PrettyPrintJsonString(), it->PrettyPrintJsonString());
    }

    // Testing values
    std::vector<std::string> expected_keys ({"bool_value", "double_value", "int_value", "level1", "string_value"});
    int counter = 0;
    for(auto it=kp.begin(); it!=kp.end(); ++it) {
        KRATOS_EXPECT_EQ(it.name(), expected_keys[counter]);
        ++counter;
    }
}

//////////////////////////////////////////////////////////////////////
// Query                                                            //
//////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsInt, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_TRUE(parameters["int_value"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsInt());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsInt());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsDouble, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_FALSE(parameters["int_value"].IsDouble());
    KRATOS_EXPECT_TRUE(parameters["double_value"].IsDouble());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsDouble());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsDouble());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsDouble());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsDouble());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsDouble());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsBool, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsBool());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsBool());
    KRATOS_EXPECT_TRUE(parameters["bool_value"].IsBool());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsBool());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsBool());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsBool());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsBool());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsString, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsString());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsString());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsString());
    KRATOS_EXPECT_TRUE(parameters["string_value"].IsString());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsString());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsString());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsString());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsPath, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsPath());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsPath());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsPath());
    KRATOS_EXPECT_TRUE(parameters["string_value"].IsPath());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsPath());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsPath());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsPath());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsStringArray, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsStringArray());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsStringArray());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsStringArray());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsStringArray());
    KRATOS_EXPECT_TRUE(parameters["string_array"].IsStringArray());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsStringArray());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsStringArray());

    // Extended Cases of Valid String Arrays
    Parameters tmp = Parameters(R"({
        "valid_string_arrays" : [ 
            ["hello", "world"],
            ["array", "string"],
            ["true","2","string"]
        ],
        "false_string_arrays" : [ 
            [[]],
            [[2,3],2],
            [2,3,[2]],
            ["string",[]],
            [{"key":3},"2"],
            [2,3,{"key":3}],
            [true,"string"],
            ["5","string",2]
        ]
    })");

    for (std::size_t i = 0;  i < tmp["valid_string_arrays"].size(); ++i) {
        const auto& valid_vector = tmp["valid_string_arrays"][i];
        KRATOS_EXPECT_TRUE(valid_vector.IsStringArray());
    }

    for (std::size_t i = 0;  i < tmp["false_string_arrays"].size(); ++i) {
        const auto& false_vector = tmp["false_string_arrays"][i];
        KRATOS_EXPECT_FALSE(false_vector.IsStringArray());
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsPathArray, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsPathArray());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsPathArray());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsPathArray());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsPathArray());
    KRATOS_EXPECT_TRUE(parameters["string_array"].IsPathArray());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsPathArray());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsPathArray());

    // Extended Cases of Valid String Arrays
    Parameters tmp = Parameters(R"({
        "valid_string_arrays" : [ 
            ["hello", "world"],
            ["array", "string"],
            ["true","2","string"]
        ],
        "false_string_arrays" : [ 
            [[]],
            [[2,3],2],
            [2,3,[2]],
            ["string",[]],
            [{"key":3},"2"],
            [2,3,{"key":3}],
            [true,"string"],
            ["5","string",2]
        ]
    })");

    for (std::size_t i = 0;  i < tmp["valid_string_arrays"].size(); ++i) {
        const auto& valid_vector = tmp["valid_string_arrays"][i];
        KRATOS_EXPECT_TRUE(valid_vector.IsPathArray());
    }

    for (std::size_t i = 0;  i < tmp["false_string_arrays"].size(); ++i) {
        const auto& false_vector = tmp["false_string_arrays"][i];
        KRATOS_EXPECT_FALSE(false_vector.IsPathArray());
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsVector, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsVector());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsVector());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsVector());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsVector());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsVector());
    KRATOS_EXPECT_TRUE(parameters["vector_value"].IsVector());
    KRATOS_EXPECT_FALSE(parameters["matrix_value"].IsVector());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsMatrix, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    KRATOS_EXPECT_FALSE(parameters["int_value"].IsMatrix());
    KRATOS_EXPECT_FALSE(parameters["double_value"].IsMatrix());
    KRATOS_EXPECT_FALSE(parameters["bool_value"].IsMatrix());
    KRATOS_EXPECT_FALSE(parameters["string_value"].IsMatrix());
    KRATOS_EXPECT_FALSE(parameters["string_array"].IsMatrix());
    KRATOS_EXPECT_FALSE(parameters["vector_value"].IsMatrix());
    KRATOS_EXPECT_TRUE(parameters["matrix_value"].IsMatrix());
}

//////////////////////////////////////////////////////////////////////
// Getter                                                           //
//////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetInt, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_EQ(parameters["int_value"].GetInt(), 10);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetInt(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetInt(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetInt(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetInt(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetInt(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetDouble, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_DOUBLE_EQ(parameters["double_value"].GetDouble(), 2.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetDouble(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetDouble(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetDouble(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetDouble(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetDouble(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetBool, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_TRUE(parameters["bool_value"].GetBool());

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetBool(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetBool(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetBool(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetBool(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetBool(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetBool(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetString, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_EQ(parameters["string_value"].GetString(), "hello");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetString(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetString(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetString(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetString(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetString(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetString(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetPath, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    KRATOS_EXPECT_EQ(parameters["string_value"].GetPath(), "hello");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetPath(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetPath(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetPath(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetPath(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetPath(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetPath(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetStringArray, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    auto values = parameters["string_array"].GetStringArray();

    KRATOS_EXPECT_EQ(values[0], "hello");
    KRATOS_EXPECT_EQ(values[1], "world");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetStringArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetStringArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetStringArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetStringArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetStringArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetStringArray(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetPathArray, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();

    auto values = parameters["string_array"].GetPathArray();

    KRATOS_EXPECT_EQ(values[0], "hello");
    KRATOS_EXPECT_EQ(values[1], "world");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetPathArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetPathArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetPathArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetPathArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetPathArray(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetPathArray(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetVector, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    const auto vector = parameters["vector_value"].GetVector();
    KRATOS_EXPECT_EQ(vector.size(), 3);
    KRATOS_EXPECT_DOUBLE_EQ(vector[0], 5.2);
    KRATOS_EXPECT_DOUBLE_EQ(vector[1], -3.1);
    KRATOS_EXPECT_DOUBLE_EQ(vector[2], 4.33);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetVector(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetVector(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetVector(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetVector(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetVector(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["matrix_value"].GetVector(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetMatrix, KratosCoreFastSuite)
{
    const auto parameters = MakeParameters();
    const auto matrix = parameters["matrix_value"].GetMatrix();
    KRATOS_EXPECT_EQ(matrix.size1(), 3);
    KRATOS_EXPECT_EQ(matrix.size2(), 2);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(0, 1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(1, 0), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(1, 1), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(2, 0), 5.0);
    KRATOS_EXPECT_DOUBLE_EQ(matrix(2, 1), 6.0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["int_value"].GetMatrix(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["double_value"].GetMatrix(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["bool_value"].GetMatrix(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_value"].GetMatrix(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["string_array"].GetMatrix(), "");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(parameters["vector_value"].GetMatrix(), "");
}

//////////////////////////////////////////////////////////////////////
// Setter                                                           //
//////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetInt, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    parameters["value"].SetInt(10);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":10})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetDouble, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    parameters["value"].SetDouble(2.0);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":2.0})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetBool, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    parameters["value"].SetBool(true);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":true})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetString, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    parameters["value"].SetString("hello");

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":"hello"})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetPath, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    const auto path = std::filesystem::path("mesh.mdpa");
    parameters["value"].SetPath(path);

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":"mesh.mdpa"})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetVector, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    Vector vector(3);
    vector[0] = 5.2;
    vector[1] = -3.1;
    vector[2] = 4.33;
    parameters["value"].SetVector(vector);

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":[5.2,-3.1,4.33]})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetMatrix, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    
    Matrix matrix(3, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(1, 0) = 3.0;
    matrix(1, 1) = 4.0;
    matrix(2, 0) = 5.0;
    matrix(2, 1) = 6.0;
    parameters["value"].SetMatrix(matrix);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":[[1.0,2.0],[3.0,4.0],[5.0,6.0]]})");
}

//////////////////////////////////////////////////////////////////////
// Setter Override                                                  //
//////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideInt, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":5})");
    
    parameters["value"].SetInt(10);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":10})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideDouble, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":1.0})");
    
    parameters["value"].SetDouble(2.0);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":2.0})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideBool, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":false})");
    
    parameters["value"].SetBool(true);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":true})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideString, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":"bye"})");
    
    parameters["value"].SetString("hello");

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":"hello"})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverridePath, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":"mosh.mpda"})");
    
    const auto path = std::filesystem::path("mesh.mdpa");
    parameters["value"].SetPath(path);

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":"mesh.mdpa"})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideVector, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":["-2.5","1.3","-33.4"]})");
    
    Vector vector(3);
    vector[0] = 5.2;
    vector[1] = -3.1;
    vector[2] = 4.33;
    parameters["value"].SetVector(vector);

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":[5.2,-3.1,4.33]})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersOverrideMatrix, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":[[0.1,0.2],[0.3,0.4],[0.5,0.6]]})");
    
    Matrix matrix(3, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(1, 0) = 3.0;
    matrix(1, 1) = 4.0;
    matrix(2, 0) = 5.0;
    matrix(2, 1) = 6.0;
    parameters["value"].SetMatrix(matrix);
    
    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":[[1.0,2.0],[3.0,4.0],[5.0,6.0]]})");
}

//////////////////////////////////////////////////////////////////////
// Adder                                                            //
//////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddInt, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");
    
    std::string key = "int";
    tmp.AddInt(key, 10);
    
    KRATOS_EXPECT_EQ(tmp[key].GetInt(),10);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddDuble, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");
    
    std::string key = "double";
    tmp.AddDouble(key, 2.0);

    KRATOS_EXPECT_DOUBLE_EQ(tmp[key].GetDouble(),2.0);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddBool, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");
    
    std::string key = "bool";
    tmp.AddBool(key, true);

    KRATOS_EXPECT_EQ(tmp[key].GetBool(),true);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddString, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");

    std::string key = "string";
    tmp.AddString(key, "hello");

    KRATOS_EXPECT_EQ(tmp[key].GetString(),"hello");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddVector, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");
    
    std::string key = "vector";
    Vector vector = ZeroVector(3);
    vector[0] = 5.2;
    vector[1] = -3.1;
    vector[2] = 4.33;
    tmp.AddVector(key, vector);
    const auto& V = tmp[key].GetVector();

    KRATOS_EXPECT_VECTOR_EQ(V,vector);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddMatrix, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({})");
    
    std::string key = "matrix";
    Matrix matrix = ZeroMatrix(3,2);
    matrix(0,0) = 1.0;
    matrix(0,1) = 2.0;
    matrix(1,0) = 3.0;
    matrix(1,1) = 4.0;
    matrix(2,0) = 5.0;
    matrix(2,1) = 6.0;
    tmp.AddMatrix(key, matrix);
    const auto& A = tmp[key].GetMatrix();
    
    KRATOS_EXPECT_MATRIX_EQ(A,matrix);
}

//////////////////////////////////////////////////////////////////////
// Corner Cases                                                     //
//////////////////////////////////////////////////////////////////////

// Unicode specific tests

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetUnicodeString, KratosCoreFastSuite)
{
    const Parameters parameters(R"({"value":"\u548c\u8c10"})"); // 和谐
    KRATOS_EXPECT_EQ(parameters["value"].GetString(), "\xE5\x92\x8C\xE8\xB0\x90");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetUnicodeString, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    const std::string value = "\xE4\xBC\xA0\xE7\xBB\x9F"; // 传统
    parameters["value"].SetString(value);

    const Parameters serialized(parameters.WriteJsonString());
    KRATOS_EXPECT_EQ(serialized["value"].GetString(), value);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsUnicodeString, KratosCoreFastSuite)
{
    const Parameters parameters(R"({"value":"\u548c\u8c10"})"); // 和谐
    KRATOS_EXPECT_TRUE(parameters["value"].IsString());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetUnicodePath, KratosCoreFastSuite)
{
    const Parameters parameters(R"({"input_filename":"\u548c\u8c10/\u4f20\u7edf.mdpa"})"); // 和谐/传统.mdpa
    const auto expected = std::filesystem::path(u"\u548c\u8c10") / u"\u4f20\u7edf.mdpa";

    KRATOS_EXPECT_EQ(parameters["input_filename"].GetPath(), expected);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetUnicodePath, KratosCoreFastSuite)
{
    Parameters parameters(R"({"input_filename":null})");
    const auto filename = std::filesystem::path(u"\u4f20\u7edf.mdpa"); // 传统.mdpa
    parameters["input_filename"].SetPath(filename);

    const Parameters serialized(parameters.WriteJsonString());
    KRATOS_EXPECT_EQ(serialized["input_filename"].GetString(), "\xE4\xBC\xA0\xE7\xBB\x9F.mdpa");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsUnicodePath, KratosCoreFastSuite)
{
    const Parameters parameters(R"({"input_filename":"\u548c\u8c10/\u4f20\u7edf.mdpa"})"); // 和谐/传统.mdpa
    KRATOS_EXPECT_TRUE(parameters["input_filename"].IsPath());
}

// Vectros with different shapes

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsVectorWithDifferentShapes, KratosCoreFastSuite)
{
    const auto tmp = MakeParametersVectorInterface();

    // Check the IsVector Method
    for (std::size_t i = 0;  i < tmp["valid_vectors"].size(); ++i) {
        const auto& valid_vector = tmp["valid_vectors"][i];
        KRATOS_EXPECT_TRUE(valid_vector.IsVector());
    }

    for (std::size_t i = 0;  i < tmp["false_vectors"].size(); ++i) {
        const auto& false_vector = tmp["false_vectors"][i];
        KRATOS_EXPECT_FALSE(false_vector.IsVector());
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetVectorWithDifferentShapes, KratosCoreFastSuite)
{
    const auto tmp = MakeParametersVectorInterface();

    // Check the GetVector Method on the valid Matrices
    for (std::size_t i = 0;  i < tmp["valid_vectors"].size(); ++i) {
        const auto& valid_vector = tmp["valid_vectors"][i];
        valid_vector.GetVector();
    }

    // Check that the errors of the GetVector method are thrown correctly
    for (std::size_t i = 0;  i < tmp["false_vectors"].size(); ++i) {
        const auto& false_vector = tmp["false_vectors"][i];
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(false_vector.GetVector(), "");
    }
}

// Set on empty

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetVectorOnEmptyValue, KratosCoreFastSuite)
{
    // Manually assign and check a Vector
    Parameters tmp(R"({})");

    Vector vec = ZeroVector(3);
    vec[0] = 1.32;
    vec[1] = -2.22;
    vec[2] = 5.5;

    tmp.AddEmptyValue("vector_value");
    tmp["vector_value"].SetVector(vec);

    const auto V2 = tmp["vector_value"].GetVector();
    KRATOS_EXPECT_VECTOR_EQ(V2,vec);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetVectorArrayInterface, KratosCoreFastSuite)
{
    // Valid input: 3D vectors
    Parameters tmp(R"({
       "not_array":{},
       "empty_vectors": [],
        "not_vectors": 5,
        "first_not_array": [
            1.0,
            [2.0, 3.0]
        ],
        "vectors3D": [
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0]
        ],
        "vectors2D": [
            [1.0, 2.0],
            [3.0, 4.0]
        ],
        "bad_size": [
            [1.0, 2.0, 3.0],
            [4.0, 5.0]
        ],
        "inner_not_array": [
            [1.0, 2.0, 3.0],
            5.0
        ],
        "not_number": [
            [1.0, 2.0, "foo"],
            [4.0, 5.0, 6.0]
        ]
    })");

    const auto vectors3D = tmp["vectors3D"].GetVectorArray();
    KRATOS_EXPECT_EQ(vectors3D.size(), 3);
    for (std::size_t i = 0; i < vectors3D.size(); ++i) {
        KRATOS_EXPECT_EQ(vectors3D[i].size(), 3);
        for (std::size_t j = 0; j < 3; ++j) {
            KRATOS_EXPECT_DOUBLE_EQ(vectors3D[i][j], 3.0 * i + j + 1.0);
        }
    }

    const auto vectors2D = tmp["vectors2D"].GetVectorArray();
    KRATOS_EXPECT_EQ(vectors2D.size(), 2);
    for (std::size_t i = 0; i < vectors2D.size(); ++i) {
        KRATOS_EXPECT_EQ(vectors2D[i].size(), 2);
    }

    KRATOS_EXPECT_EQ(tmp["empty_vectors"].GetVectorArray().size(), 0);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["not_array"].GetVectorArray(),
        "Argument must be an array of arrays"
    );

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["not_vectors"].GetVectorArray(),
        "Argument must be an array of arrays"
    );

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["first_not_array"].GetVectorArray(),
        "Entry 0 is not an array"
    );

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["bad_size"].GetVectorArray(),
        "All inner arrays must have the same size"
    );

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["inner_not_array"].GetVectorArray(),
        "Entry 1 is not an array"
    );

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        tmp["not_number"].GetVectorArray(),
        "Entry (0,2) is not a number"
    );
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsMatrixWithDifferentShapes, KratosCoreFastSuite)
{
    const auto tmp = MakeMatrixShapeParameters();

    for (std::size_t i = 0;  i < tmp["valid_matrices"].size(); ++i) {
        const auto& valid_matrix = tmp["valid_matrices"][i];
        KRATOS_EXPECT_TRUE(valid_matrix.IsMatrix());
    }

    for (std::size_t i = 0;  i < tmp["false_matrices"].size(); ++i) {
        const auto& false_matrix = tmp["false_matrices"][i];
        KRATOS_EXPECT_FALSE(false_matrix.IsMatrix());
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetMatrixWithDifferentShapes, KratosCoreFastSuite)
{
    const auto tmp = MakeMatrixShapeParameters();

    for (std::size_t i = 0;  i < tmp["valid_matrices"].size(); ++i) {
        const auto& valid_matrix = tmp["valid_matrices"][i];
        valid_matrix.GetMatrix();
    }

    // Check that the errors of the GetMatrix method are thrown correctly
    for (std::size_t i = 0;  i < tmp["false_matrices"].size(); ++i) {
        const auto& false_matrix = tmp["false_matrices"][i];
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(false_matrix.GetMatrix(), "");
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetMatrixOnEmptyValue, KratosCoreFastSuite)
{
    Parameters tmp(R"({})");

    Matrix mat = ZeroMatrix(3,2);
    mat(0,0) = 1.0;
    mat(0,1) = 2.0;
    mat(1,0) = 3.0;
    mat(1,1) = 4.0;
    mat(2,0) = 5.0;
    mat(2,1) = 6.0;

    tmp.AddEmptyValue("matrix_value");
    tmp["matrix_value"].SetMatrix(mat);

    const auto& A2 = tmp["matrix_value"].GetMatrix();
    KRATOS_EXPECT_MATRIX_EQ(A2, mat);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersNullvsNullValidation, KratosCoreFastSuite)
{
    // Supplied settings
    Parameters null_custom = Parameters(R"({
        "parameter": null
    })");

    // Default settings
    Parameters null_default = Parameters(R"({
        "parameter": null
    })");

    // This should NOT raise, hence making the test to pass
    null_custom.ValidateAndAssignDefaults(null_default);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersDoublevsNullValidation, KratosCoreFastSuite)
{
    // Supplied settings
    Parameters double_custom = Parameters(R"({
        "parameter": 0.0
    })");

    // Default settings
    Parameters null_default = Parameters(R"({
        "parameter": null
    })");

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(double_custom.ValidateAndAssignDefaults(null_default), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetStringArrayValid, KratosCoreFastSuite)
{
    Parameters tmp = Parameters(R"({
        "parameter": ["foo", "bar"]
    })");
    auto v = tmp["parameter"].GetStringArray();
    KRATOS_EXPECT_EQ(v.size(), 2);
    KRATOS_EXPECT_EQ(v[0], "foo");
    KRATOS_EXPECT_EQ(v[1], "bar");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetStringArrayValid, KratosCoreFastSuite)
{
    Parameters parameters(R"({"value":null})");
    parameters["value"].SetStringArray({"foo", "bar"});

    KRATOS_EXPECT_EQ(parameters.WriteJsonString(), R"({"value":["foo","bar"]})");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWithIncludes, KratosCoreFastSuite)
{
    ScopedFile included_json("test_included_parameters.json");
    ScopedFile included_json_level2("test_included_parameters_level2.json");

    included_json << GetIncludedJSONString();
    included_json_level2 << GetIncludedJSONLevel2String();

    Parameters kp = Parameters(GetJSONStringWithIncludes());
    KRATOS_EXPECT_EQ(
        kp.WriteJsonString(),
        R"({"bool_value":true,"double_value":2.0,"int_value":10,"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"})"
    );
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWithRepeatedIncludes, KratosCoreFastSuite)
{
    ScopedFile included_json("test_included_parameters.json");
    ScopedFile included_json_level2("test_included_parameters_level2.json");

    included_json << GetIncludedJSONString();
    included_json_level2 << GetIncludedJSONLevel2String();

    Parameters parameters(R"({
        "another_include" : {
            "@include_json" : "test_included_parameters.json"
        },
        "@include_json" : "test_included_parameters.json"
    })");
    KRATOS_EXPECT_EQ(
        parameters.WriteJsonString(),
        R"({"another_include":{"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"},"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"})"
    );
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWithSelfInclude, KratosCoreFastSuite)
{
    ScopedFile file_0_includes_0("test_cyclic_0_0.json");
    file_0_includes_0 << GetCircularIncludeJSONString(0, 0);

    try {
        Parameters(R"({"@include_json" : "test_cyclic_0_0.json"})");
    } catch (Exception& rException) { // std::exceptions are not caught and indicate parsing errors
        KRATOS_EXPECT_NE(std::string(rException.what()).find("cycle in json"), std::string::npos);
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWithCyclicInclude, KratosCoreFastSuite)
{
    ScopedFile file_0_includes_1("test_cyclic_0_1.json");
    ScopedFile file_1_includes_2("test_cyclic_1_2.json");
    ScopedFile file_2_includes_0("test_cyclic_2_0.json");

    file_0_includes_1 << GetCircularIncludeJSONString(0, 1);
    file_1_includes_2 << GetCircularIncludeJSONString(1, 2);
    file_2_includes_0 << GetCircularIncludeJSONString(2, 0);

    try {
        Parameters(R"({"@include_json" : "test_cyclic_0_1.json"})");
    } catch (Exception& rException) { // std::exceptions are not caught and indicate parsing errors
        KRATOS_EXPECT_NE(std::string(rException.what()).find("cycle in json"), std::string::npos);
    }
}

}  // namespace Testing.
}  // namespace Kratos.
