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

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_parameters.h"
#include "testing/scoped_file.h"

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

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsMethods, KratosCoreFastSuite)
{
    // This method checks all the "IsXXX" Methods
    Parameters tmp = Parameters(R"({
        "int_value" : 10, /* This is comment to check that comments work */
        "double_value": 2.0, // This is comment too, but using another comment
        "bool_value" : true, // This is another comment being meta as realizing that all the possibilities are already check
        "string_value" : "hello",/* This is a nihilist comment about the futile existence of the previous comment as a metacomment */
        "s_array_value" : ["hello", "world"],
        "vector_value" : [5,3,4],
        "matrix_value" : [[1,2],[3,6]]
    })"); // if you add more values to this, make sure to add the corresponding in the loop

    for(auto it=tmp.begin(); it!=tmp.end(); ++it) {
        const std::string key = it.name();

        if (key.find("int") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsInt());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsInt());
        }

        if (key.find("double") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsDouble());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsDouble());
        }

        if (key.find("bool") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsBool());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsBool());
        }

        if (key.find("string") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsString());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsString());
        }

        if (key.find("s_array") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsStringArray());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsStringArray());
        }

        if (key.find("vector") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsVector());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsVector());
        }

        if (key.find("matrix") != std::string::npos) {
            KRATOS_EXPECT_TRUE(tmp[key].IsMatrix());
        } else {
            KRATOS_EXPECT_FALSE(tmp[key].IsMatrix());
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetMethods, KratosCoreFastSuite)
{
    // This method checks all the "GetXXX" Methods if they throw an error
    Parameters tmp = Parameters(R"({
        "int_value" : 10,
        "double_value": 2.0,
        "bool_value" : true,
        "string_value" : "hello",
        "vector_value" : [5.2,-3.1,4.33],
        "matrix_value" : [[1,2],[3,4],[5,6]]
    })"); // if you add more values to this, make sure to add the corresponding in the loop

    for(auto it=tmp.begin(); it!=tmp.end(); ++it) {
        const std::string key = it.name();

        // Int and Double are checked tgth bcs both internally call "IsNumber"
        if (key.find("double") != std::string::npos || key.find("int") != std::string::npos) {
            if (key.find("int") != std::string::npos) {
                KRATOS_EXPECT_EQ(tmp[key].GetInt(),10);
            }
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetInt(), "");
        }

        if (key.find("double") != std::string::npos || key.find("int") != std::string::npos) {
            if (key.find("double") != std::string::npos) {
                KRATOS_EXPECT_DOUBLE_EQ(tmp[key].GetDouble(),2.0);
            }
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetDouble(), "");
        }

        if (key.find("bool") != std::string::npos) {
            KRATOS_EXPECT_EQ(tmp[key].GetBool(), true);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetBool(), "");
        }

        if (key.find("string") != std::string::npos) {
            KRATOS_EXPECT_EQ(tmp[key].GetString(),"hello");
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetString(), "");
        }

        if (key.find("vector") != std::string::npos) {
            const auto& V = tmp[key].GetVector();
            KRATOS_EXPECT_DOUBLE_EQ(V[0],5.2);
            KRATOS_EXPECT_DOUBLE_EQ(V[1],-3.1);
            KRATOS_EXPECT_DOUBLE_EQ(V[2],4.33);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetVector(), "");
        }

        if (key.find("matrix") != std::string::npos) {
            const auto& A = tmp[key].GetMatrix();
            KRATOS_EXPECT_DOUBLE_EQ(A(0,0), 1.0);
            KRATOS_EXPECT_DOUBLE_EQ(A(0,1), 2.0);
            KRATOS_EXPECT_DOUBLE_EQ(A(1,0), 3.0);
            KRATOS_EXPECT_DOUBLE_EQ(A(1,1), 4.0);
            KRATOS_EXPECT_DOUBLE_EQ(A(2,0), 5.0);
            KRATOS_EXPECT_DOUBLE_EQ(A(2,1), 6.0);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetMatrix(), "");
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetMethods, KratosCoreFastSuite)
{
    // This method checks all the "GetXXX" Methods if they throw an error
    Parameters tmp = Parameters(R"({
        "int_value" : 0,
        "double_value": 0.0,
        "bool_value" : false,
        "string_value" : "",
        "vector_value" : [],
        "matrix_value" : [[0]]
    })"); // if you add more values to this, make sure to add the corresponding in the loop

    for(auto it=tmp.begin(); it!=tmp.end(); ++it) {
        const std::string key = it.name();

        // Int and Double are checked tgth bcs both internally call "IsNumber"
        if (key.find("double") != std::string::npos || key.find("int") != std::string::npos) {
            if (key.find("int") != std::string::npos) {
                tmp[key].SetInt(10);
                KRATOS_EXPECT_EQ(tmp[key].GetInt(),10);
            }
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetInt(), "");
        }

        if (key.find("double") != std::string::npos || key.find("int") != std::string::npos) {
            if (key.find("double") != std::string::npos) {
                tmp[key].SetDouble(2.0);
                KRATOS_EXPECT_EQ(tmp[key].GetDouble(),2.0);
            }
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetDouble(), "");
        }

        if (key.find("bool") != std::string::npos) {
            tmp[key].SetBool(true);
            KRATOS_EXPECT_EQ(tmp[key].GetBool(),true);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetBool(), "");
        }

        if (key.find("string") != std::string::npos) {
            tmp[key].SetString("hello");
            KRATOS_EXPECT_EQ(tmp[key].GetString(),"hello");
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetString(), "");
        }

        if (key.find("vector") != std::string::npos) {
            Vector vector = ZeroVector(3);
            vector[0] = 5.2;
            vector[1] = -3.1;
            vector[2] = 4.33;
            tmp[key].SetVector(vector);
            const auto& V = tmp[key].GetVector();
            KRATOS_EXPECT_VECTOR_EQ(V,vector);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetVector(), "");
        }

        if (key.find("matrix") != std::string::npos) {
            Matrix matrix = ZeroMatrix(3,2);
            matrix(0,0) = 1.0;
            matrix(0,1) = 2.0;
            matrix(1,0) = 3.0;
            matrix(1,1) = 4.0;
            matrix(2,0) = 5.0;
            matrix(2,1) = 6.0;
            tmp[key].SetMatrix(matrix);
            const auto& A = tmp[key].GetMatrix();
            KRATOS_EXPECT_MATRIX_EQ(A,matrix);
        } else {
            KRATOS_EXPECT_EXCEPTION_IS_THROWN(tmp[key].GetMatrix(), "");
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddMethods, KratosCoreFastSuite)
{
    // This method checks all the "GetXXX" Methods if they throw an error
    Parameters tmp = Parameters(R"({})");

    std::string key = "int";
    tmp.AddInt(key, 10);
    KRATOS_EXPECT_EQ(tmp[key].GetInt(),10);

    key = "double";
    tmp.AddDouble(key, 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(tmp[key].GetDouble(),2.0);

    key = "bool";
    tmp.AddBool(key, true);
    KRATOS_EXPECT_EQ(tmp[key].GetBool(),true);

    key = "string";
    tmp.AddString(key, "hello");
    KRATOS_EXPECT_EQ(tmp[key].GetString(),"hello");

    key = "vector";
    Vector vector = ZeroVector(3);
    vector[0] = 5.2;
    vector[1] = -3.1;
    vector[2] = 4.33;
    tmp.AddVector(key, vector);
    const auto& V = tmp[key].GetVector();
    KRATOS_EXPECT_VECTOR_EQ(V,vector);

    key = "matrix";
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

KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsStringArray, KratosCoreFastSuite)
{
    // Read and check string arrays from a Parameters-Object
    Parameters tmp = Parameters(R"({
        "valid_string_arrays" : [ ["hello", "world"],
                                ["array", "string"],
                                ["true","2","string"]
        ],
        "false_string_arrays" : [ [[]],
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

KRATOS_TEST_CASE_IN_SUITE(KratosParametersVectorInterface, KratosCoreFastSuite)
{
    // Read and check Vectors from a Parameters-Object
    Parameters tmp = Parameters(R"({
        "valid_vectors" : [ []
        ],
        "false_vectors" : [ [[]],
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

    // Check the IsVector Method
    for (std::size_t i = 0;  i < tmp["valid_vectors"].size(); ++i) {
        const auto& valid_vector = tmp["valid_vectors"][i];
        KRATOS_EXPECT_TRUE(valid_vector.IsVector());
    }

    for (std::size_t i = 0;  i < tmp["false_vectors"].size(); ++i) {
        const auto& false_vector = tmp["false_vectors"][i];
        KRATOS_EXPECT_FALSE(false_vector.IsVector());
    }

    // Check the GetVector Method also on the valid Matrices
    for (std::size_t i = 0;  i < tmp["valid_vectors"].size(); ++i) {
        const auto& valid_vector = tmp["valid_vectors"][i];
        valid_vector.GetVector();
    }

    // Check that the errors of the GetVector method are thrown correctly
    for (std::size_t i = 0;  i < tmp["false_vectors"].size(); ++i) {
        const auto& false_vector = tmp["false_vectors"][i];
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(false_vector.GetVector(), "");
    }

    // Manually assign and check a Vector
    Vector vec = ZeroVector(3);
    vec[0] = 1.32;
    vec[1] = -2.22;
    vec[2] = 5.5;

    tmp.AddEmptyValue("vector_value");
    tmp["vector_value"].SetVector(vec);

    KRATOS_EXPECT_TRUE(tmp["vector_value"].IsVector());

    const auto V2 = tmp["vector_value"].GetVector();
    KRATOS_EXPECT_VECTOR_EQ(V2,vec);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersMatrixInterface, KratosCoreFastSuite)
{
    // Read and check Matrices from a Parameters-Object
    Parameters tmp = Parameters(R"({
        "valid_matrices" : [ [[]],
                              [[],[]],
                              [[-9.81,8, 5.47]]
        ],
        "false_matrices" : [ [],
                              [[[]]],
                              [[3.3] , [1,2]],
                              [[2,1.5,3.3] , [3,{"key":3},2]],
                              [[2,1.5,3.3] , [5,false,2]],
                              [[2,1.5,3.3] , [[2,3],1,2]],
                              [[2,1.5,3.3] , ["string",2,9]]
        ]
    })");

    // Check the IsMatrix Method
    for (std::size_t i = 0;  i < tmp["valid_matrices"].size(); ++i) {
        const auto& valid_matrix = tmp["valid_matrices"][i];
        KRATOS_EXPECT_TRUE(valid_matrix.IsMatrix());
    }

    for (std::size_t i = 0;  i < tmp["false_matrices"].size(); ++i) {
        const auto& false_matrix = tmp["false_matrices"][i];
        KRATOS_EXPECT_FALSE(false_matrix.IsMatrix());
    }

    // Check the GetMatrix Method also on the valid Matrices
    for (std::size_t i = 0;  i < tmp["valid_matrices"].size(); ++i) {
        const auto& valid_matrix = tmp["valid_matrices"][i];
        valid_matrix.GetMatrix();
    }

    // Check that the errors of the GetMatrix method are thrown correctly
    for (std::size_t i = 0;  i < tmp["false_matrices"].size(); ++i) {
        const auto& false_matrix = tmp["false_matrices"][i];
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(false_matrix.GetMatrix(), "");
    }

    // Manually assign and check a Matrix
    Matrix mat = ZeroMatrix(3,2);
    mat(0,0) = 1.0;
    mat(0,1) = 2.0;
    mat(1,0) = 3.0;
    mat(1,1) = 4.0;
    mat(2,0) = 5.0;
    mat(2,1) = 6.0;

    tmp.AddEmptyValue("matrix_value");
    tmp["matrix_value"].SetMatrix(mat);

    KRATOS_EXPECT_TRUE(tmp["matrix_value"].IsMatrix());

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

KRATOS_TEST_CASE_IN_SUITE(KratosParametersGeStringArrayValid, KratosCoreFastSuite)
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
    Parameters initial = Parameters(R"({
        "parameter": ["foo", "bar"]
    })");
    auto string_array = initial["parameter"].GetStringArray();

    Parameters new_param = Parameters();
    new_param.AddEmptyValue("new_parameter");
    new_param["new_parameter"].SetStringArray(string_array);

    auto new_string_array = initial["parameter"].GetStringArray();

    int counter = 0;
    for (auto& r_string : string_array) {
        KRATOS_EXPECT_EQ(new_string_array[counter], r_string);
        ++counter;
    }
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
