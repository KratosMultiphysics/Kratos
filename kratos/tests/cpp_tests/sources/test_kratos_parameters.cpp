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
#include "includes/kratos_parameters.h"

namespace Kratos {
namespace Testing {

// input string with ugly formatting
std::string GetJSONString()
{
    const std::string json_string = R"(
    {
      "bool_value" : true, "double_value": 2.0, "int_value" : 10,
      "level1":
      {
        "list_value":[ 3, "hi", false],
        "tmp" : 5.0
      },
      "string_value" : "hello"
    })";
    return json_string;
}

std::string GetJSONStringPrettyOut()
{
    const std::string pretty_out =  R"({
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
    return pretty_out;
}

std::string GetJSONStringPrettyOutAfterChange()
{
    const std::string pretty_out_after_change = R"({
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

    return pretty_out_after_change;
}

// here the level1 var is set to a double so that a validation error should be thrown
std::string GetJSONStringWrongType()
{
    const std::string wrong_type = R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": 0.0,
        "string_value": "hello"
    })";

    return wrong_type;
}

// int value is badly spelt
std::string GetJSONStringWrongSpelling()
{
    const std::string wrong_spelling = R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_values": 10,
        "level1": 0.0,
        "string_value": "hello"
    })";

    return wrong_spelling;
}

// wrong on the first level
// error shall be only detective by recursive validation
std::string GetJSONStringWrongLevel2()
{
    const std::string wrong_lev2 = R"({
        "bool_value": true,
        "double_value": 2.0,
        "int_value": 10,
        "level1": { "a":0.0 },
        "string_value": "hello"
    })";
    return wrong_lev2;
}

std::string GetJSONStringDefaults()
{
    const std::string defaults = R"({
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
    return defaults;
}

std::string GetJSONStringIncomplete()
{
    const std::string incomplete = R"({
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
    return incomplete;
}

std::string GetJSONStringIncompleteWithExtraParameter()
{
    const std::string incomplete_with_extra_parameter = R"({
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
    return incomplete_with_extra_parameter;
}

std::string GetJSONStringExpectedValidationOutput()
{
    const std::string expected_validation_output = R"({
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
    return expected_validation_output;
}

std::string GetJSONStringFourLevels()
{
    const std::string four_levels = R"({
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
    return four_levels;
}

std::string GetJSONStringForLevelsVariation()
{
    const std::string four_levels_variation = R"({
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
    return four_levels_variation;
}

std::string GetJSONStringForLevelsWrongVariation()
{
    const std::string four_levels_wrong_variation = R"({
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
    return four_levels_wrong_variation;
}

std::string GetJSONStringForLevelsDefaults()
{
    const std::string four_levels_defaults = R"({
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
    return four_levels_defaults;
}

KRATOS_TEST_CASE_IN_SUITE(KratosParameters, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONString());
    KRATOS_CHECK_STRING_EQUAL(
        kp.WriteJsonString(),
        R"({"bool_value":true,"double_value":2.0,"int_value":10,"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"})"
    );

    KRATOS_CHECK(kp.Has("int_value"));
    KRATOS_CHECK_IS_FALSE(kp.Has("unextisting_value"));

    KRATOS_CHECK_EQUAL(kp["int_value"].GetInt(), 10);
    KRATOS_CHECK_EQUAL(kp["double_value"].GetDouble(), 2.0);
    KRATOS_CHECK_EQUAL(kp["bool_value"].GetBool(), true);
    KRATOS_CHECK_EQUAL(kp["string_value"].GetString(), "hello");

    KRATOS_CHECK_STRING_EQUAL(kp.PrettyPrintJsonString(), GetJSONStringPrettyOut());
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersChangeParameters, KratosCoreFastSuite)
{
    // Now change one item in the sublist
    Parameters kp = Parameters(GetJSONString());
    Parameters subparams = kp["level1"];

    Parameters my_list = subparams["list_value"];

    for (auto& r_param : my_list) {
        if (r_param.IsBool()) {
            KRATOS_CHECK_IS_FALSE(r_param.GetBool())
        }
    }

    // my_list = subparams["list_value"]
    subparams["list_value"][0].SetString("changed");

    KRATOS_CHECK_STRING_EQUAL(
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

    KRATOS_CHECK_STRING_EQUAL(
        other_copy.PrettyPrintJsonString(),
        original_out
    );

    other_copy["int_value"].SetInt(-1);
    KRATOS_CHECK_EQUAL(kp["int_value"].GetInt(), 10);
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersWrongParameters, KratosCoreFastSuite)
{
    // Should check which errors are thrown!!
    Parameters kp = Parameters(GetJSONString());
    KRATOS_CHECK_EXCEPTION_IS_THROWN(kp["no_value"].GetInt(), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsDueToWrongTypes, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongType());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_CHECK_EXCEPTION_IS_THROWN(kp.ValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsDueToWrongSpelling, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongSpelling());
    Parameters  defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_CHECK_EXCEPTION_IS_THROWN(kp.ValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidationFailsErrorsOnFirstLevel, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringWrongLevel2());
    Parameters defaults_params = Parameters(GetJSONStringDefaults());

    // Should check which errors are thrown!!
    KRATOS_CHECK_EXCEPTION_IS_THROWN(kp.RecursivelyValidateAndAssignDefaults(defaults_params), "");
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursiveValidation4Levels, KratosCoreFastSuite)
{
    Parameters kp = Parameters(GetJSONStringFourLevels());
    Parameters kp_variation = Parameters(GetJSONStringForLevelsVariation());
    Parameters kp_wrong_wariation = Parameters(GetJSONStringForLevelsWrongVariation());
    Parameters defaults_params = Parameters(GetJSONStringForLevelsDefaults());

    kp.RecursivelyValidateAndAssignDefaults(defaults_params);
    kp_variation.RecursivelyValidateAndAssignDefaults(defaults_params);

    KRATOS_CHECK( kp.IsEquivalentTo(defaults_params) );
    KRATOS_CHECK_IS_FALSE( kp_variation.IsEquivalentTo(defaults_params) );

    KRATOS_CHECK( kp.HasSameKeysAndTypeOfValuesAs(defaults_params) );
    KRATOS_CHECK( kp_variation.HasSameKeysAndTypeOfValuesAs(defaults_params) );
    KRATOS_CHECK_IS_FALSE( kp_wrong_wariation.HasSameKeysAndTypeOfValuesAs(defaults_params) );
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
    KRATOS_CHECK_STRING_EQUAL(kp.PrettyPrintJsonString(), GetJSONStringExpectedValidationOutput());

    KRATOS_CHECK_DOUBLE_EQUAL(kp["level1"]["tmp"].GetDouble(), 5.0);  // not 2, since kp overwrites the defaults
}

KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddMissingParameters, KratosCoreFastSuite)
{
    // Only missing parameters are added, no complaints if there already exist more than in the defaults
    Parameters kp = Parameters(GetJSONString());
    Parameters tmp = Parameters(GetJSONStringIncompleteWithExtraParameter());

    kp.AddMissingParameters(tmp);

    KRATOS_CHECK_STRING_EQUAL(kp["new_default_obj"]["aaa"].GetString(), "string");
    KRATOS_CHECK_STRING_EQUAL(kp["string_value"].GetString(), "hello");
    KRATOS_CHECK_IS_FALSE(kp["level1"].Has("new_sublevel"));
}

// KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursivelyAddMissingParameters, KratosCoreFastSuite)
// {
//     // only missing parameters are added, no complaints if there already exist more than in the defaults
//     Parameters kp = Parameters(GetJSONString());
//     Parameters tmp = Parameters(GetJSONStringIncompleteWithExtraParameter());
//
//     kp.RecursivelyAddMissingParameters(tmp)
//
//     KRATOS_CHECK(kp["level1"].Has("new_sublevel"))
//     KRATOS_CHECK_EQUAL(kp["level1"]["new_sublevel"].GetString(), "this should only be assigned in recursive")
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersValidateDefaults, KratosCoreFastSuite)
// {
//     // only parameters from defaults are validated, no new values are added
//     Parameters kp = Parameters(incomplete_with_extra_parameter)
//     tmp = Parameters(defaults)
//
//     kp.ValidateDefaults(tmp)
//
//     KRATOS_CHECK_IS_FALSE(kp.Has("bool_value"))
//     KRATOS_CHECK_IS_FALSE(kp.Has("double_value"))
//     KRATOS_CHECK(kp.Has("level1"))
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersRecursivelyValidateDefaults, KratosCoreFastSuite)
// {
//     // only parameters from defaults are validated, no new values are added
//     Parameters kp = Parameters(incomplete)
//     tmp = Parameters(defaults)
//
//     kp.RecursivelyValidateDefaults(tmp)
//
//     KRATOS_CHECK_IS_FALSE(kp.Has("bool_value"))
//     KRATOS_CHECK_IS_FALSE(kp.Has("double_value"))
//     KRATOS_CHECK(kp.Has("level1"))
// }
//
// def test_recursively_validate_defaults_fails(self):
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetStringArrayValid, KratosCoreFastSuite)
// {
//     // only parameters from defaults are validated, no new values are added
//     Parameters kp = Parameters(incomplete_with_extra_parameter)
//     tmp = Parameters(defaults)
//
//     with self.assertRaises(RuntimeError):
//         kp.RecursivelyValidateDefaults(tmp)
//
//     // sub_level
//     KRATOS_CHECK_IS_FALSE(kp["level1"].Has("tmp"))
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddValue, KratosCoreFastSuite)
// {
//     Parameters kp = Parameters("{}")
//     kp.AddEmptyValue("new_double").SetDouble(1.0)
//
//     KRATOS_CHECK(kp.Has("new_double"))
//     KRATOS_CHECK_EQUAL(kp["new_double"].GetDouble(), 1.0)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddEmptyArray, KratosCoreFastSuite)
// {
//     Parameters kp = Parameters("{}")
//     kp.AddEmptyArray("new_array")
//
//     KRATOS_CHECK(kp.Has("new_array"))
//     KRATOS_CHECK_EQUAL(kp["new_array"].size(), 0)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersIterators, KratosCoreFastSuite)
// {
//     Parameters kp = Parameters(GetJSONString());
//
//     //iteration by range
//     nitems = 0
//     for iterator in kp:
//         nitems = nitems + 1
//     KRATOS_CHECK_EQUAL(nitems, 5)
//
//     //iteration by items
//     for key,value in kp.items():
//         //print(value.PrettyPrintJsonString())
//         KRATOS_CHECK_EQUAL(kp[key].PrettyPrintJsonString(), value.PrettyPrintJsonString())
//         //print(key,value)
//
//     //testing values
//     expected_values = ['true', '2.0', '10', '{"list_value":[3,"hi",false],"tmp":5.0}','"hello"']
//     counter = 0
//
//     for value in kp.values():
//         KRATOS_CHECK_EQUAL(value.WriteJsonString(), expected_values[counter])
//         counter += 1
//
//     //testing values
//     expected_keys = ['bool_value', 'double_value', 'int_value', 'level1', 'string_value']
//     counter = 0
//     for key in kp.keys():
//         KRATOS_CHECK_EQUAL(key, expected_keys[counter])
//         counter += 1
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersRemoveValue, KratosCoreFastSuite)
// {
//     Parameters kp = Parameters(GetJSONString());
//     KRATOS_CHECK(kp.Has("int_value"))
//     KRATOS_CHECK(kp.Has("level1"))
//
//     kp.RemoveValue("int_value")
//     kp.RemoveValue("level1")
//
//     KRATOS_CHECK_IS_FALSE(kp.Has("int_value"))
//     KRATOS_CHECK_IS_FALSE(kp.Has("level1"))
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersIsMethods, KratosCoreFastSuite)
// {
//     // This method checks all the "IsXXX" Methods
//     Parameters tmp = Parameters(R"({
//         "int_value" : 10, /* This is comment to check that comments work */
//         "double_value": 2.0, // This is comment too, but using another comment
//         "bool_value" : true, // This is another comment being meta as realizing that all the possibilities are already check
//         "string_value" : "hello",/* This is a nihilist comment about the futile existence of the previous comment as a metacomment */
//         "vector_value" : [5,3,4],
//         "matrix_value" : [[1,2],[3,6]]
//     })"); // if you add more values to this, make sure to add the corresponding in the loop
//
//     for key in tmp.keys():
//         val_type = key[:-6] // removing "_value"
//
//         if val_type == "int":
//             KRATOS_CHECK(tmp[key].IsInt())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsInt())
//
//         if val_type == "double":
//             KRATOS_CHECK(tmp[key].IsDouble())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsDouble())
//
//         if val_type == "bool":
//             KRATOS_CHECK(tmp[key].IsBool())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsBool())
//
//         if val_type == "string":
//             KRATOS_CHECK(tmp[key].IsString())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsString())
//
//         if val_type == "vector":
//             KRATOS_CHECK(tmp[key].IsVector())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsVector())
//
//         if val_type == "matrix":
//             KRATOS_CHECK(tmp[key].IsMatrix())
//         else:
//             KRATOS_CHECK_IS_FALSE(tmp[key].IsMatrix())
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersGetMethods, KratosCoreFastSuite)
// {
//     // This method checks all the "GetXXX" Methods if they throw an error
//     Parameters tmp = Parameters(R"({
//         "int_value" : 10,
//         "double_value": 2.0,
//         "bool_value" : true,
//         "string_value" : "hello",
//         "vector_value" : [5.2,-3.1,4.33],
//         "matrix_value" : [[1,2],[3,4],[5,6]]
//     })"); // if you add more values to this, make sure to add the corresponding in the loop
//
//     for key in tmp.keys():
//         val_type = key[:-6] // removing "_value"
//
//         // Int and Double are checked tgth bcs both internally call "IsNumber"
//         if val_type == "int" or val_type == "double":
//             if val_type == "int":
//                 KRATOS_CHECK_EQUAL(tmp[key].GetInt(),10)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetInt()
//
//         if val_type == "double" or val_type == "int":
//             if val_type == "double":
//                 KRATOS_CHECK_EQUAL(tmp[key].GetDouble(),2.0)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetDouble()
//
//         if val_type == "bool":
//             KRATOS_CHECK_EQUAL(tmp[key].GetBool(),True)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetBool()
//
//         if val_type == "string":
//             KRATOS_CHECK_EQUAL(tmp[key].GetString(),"hello")
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetString()
//
//         if val_type == "vector":
//             V = tmp[key].GetVector()
//             KRATOS_CHECK_EQUAL(V[0],5.2)
//             KRATOS_CHECK_EQUAL(V[1],-3.1)
//             KRATOS_CHECK_EQUAL(V[2],4.33)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetVector()
//
//         if val_type == "matrix":
//             A = tmp[key].GetMatrix()
//             KRATOS_CHECK_EQUAL(A[0,0],1.0)
//             KRATOS_CHECK_EQUAL(A[0,1],2.0)
//             KRATOS_CHECK_EQUAL(A[1,0],3.0)
//             KRATOS_CHECK_EQUAL(A[1,1],4.0)
//             KRATOS_CHECK_EQUAL(A[2,0],5.0)
//             KRATOS_CHECK_EQUAL(A[2,1],6.0)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetMatrix()
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetMethods KratosCoreFastSuite)
// {
//     // This method checks all the "GetXXX" Methods if they throw an error
//     Parameters tmp = Parameters(R"({
//         "int_value" : 0,
//         "double_value": 0.0,
//         "bool_value" : false,
//         "string_value" : "",
//         "vector_value" : [],
//         "matrix_value" : [[0]]
//     })"); // if you add more values to this, make sure to add the corresponding in the loop
//
//     for key in tmp.keys():
//         val_type = key[:-6] // removing "_value"
//
//         // Int and Double are checked tgth bcs both internally call "IsNumber"
//         if val_type == "int" or val_type == "double":
//             if val_type == "int":
//                 tmp[key].SetInt(10)
//                 KRATOS_CHECK_EQUAL(tmp[key].GetInt(),10)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetInt()
//
//         if val_type == "double" or val_type == "int":
//             if val_type == "double":
//                 tmp[key].SetDouble(2.0)
//                 KRATOS_CHECK_EQUAL(tmp[key].GetDouble(),2.0)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetDouble()
//
//         if val_type == "bool":
//             tmp[key].SetBool(True)
//             KRATOS_CHECK_EQUAL(tmp[key].GetBool(),True)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetBool()
//
//         if val_type == "string":
//             tmp[key].SetString("hello")
//             KRATOS_CHECK_EQUAL(tmp[key].GetString(),"hello")
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetString()
//
//         if val_type == "vector":
//             vector = Vector(3)
//             vector[0] = 5.2
//             vector[1] = -3.1
//             vector[2] = 4.33
//             tmp[key].SetVector(vector)
//             V = tmp[key].GetVector()
//             KRATOS_CHECK_EQUAL(V[0],5.2)
//             KRATOS_CHECK_EQUAL(V[1],-3.1)
//             KRATOS_CHECK_EQUAL(V[2],4.33)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetVector()
//
//         if val_type == "matrix":
//             matrix = Matrix(3,2)
//             matrix[0,0] = 1.0
//             matrix[0,1] = 2.0
//             matrix[1,0] = 3.0
//             matrix[1,1] = 4.0
//             matrix[2,0] = 5.0
//             matrix[2,1] = 6.0
//             tmp[key].SetMatrix(matrix)
//             A = tmp[key].GetMatrix()
//             KRATOS_CHECK_EQUAL(A[0,0],1.0)
//             KRATOS_CHECK_EQUAL(A[0,1],2.0)
//             KRATOS_CHECK_EQUAL(A[1,0],3.0)
//             KRATOS_CHECK_EQUAL(A[1,1],4.0)
//             KRATOS_CHECK_EQUAL(A[2,0],5.0)
//             KRATOS_CHECK_EQUAL(A[2,1],6.0)
//         else:
//             with self.assertRaises(RuntimeError):
//                 tmp[key].GetMatrix()
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersAddMethods, KratosCoreFastSuite)
// {
//     // This method checks all the "GetXXX" Methods if they throw an error
//     Parameters tmp = Parameters(R"({})");
//
//     key = "int"
//     tmp.AddInt(key, 10)
//     KRATOS_CHECK_EQUAL(tmp[key].GetInt(),10)
//
//     key = "double"
//     tmp.AddDouble(key, 2.0)
//     KRATOS_CHECK_EQUAL(tmp[key].GetDouble(),2.0)
//
//     key = "bool"
//     tmp.AddBool(key, True)
//     KRATOS_CHECK_EQUAL(tmp[key].GetBool(),True)
//
//     key = "string"
//     tmp.AddString(key, "hello")
//     KRATOS_CHECK_EQUAL(tmp[key].GetString(),"hello")
//
//     key = "vector"
//     vector = Vector(3)
//     vector[0] = 5.2
//     vector[1] = -3.1
//     vector[2] = 4.33
//     tmp.AddVector(key, vector)
//     V = tmp[key].GetVector()
//     KRATOS_CHECK_EQUAL(V[0],5.2)
//     KRATOS_CHECK_EQUAL(V[1],-3.1)
//     KRATOS_CHECK_EQUAL(V[2],4.33)
//
//     key = "matrix"
//     matrix = Matrix(3,2)
//     matrix[0,0] = 1.0
//     matrix[0,1] = 2.0
//     matrix[1,0] = 3.0
//     matrix[1,1] = 4.0
//     matrix[2,0] = 5.0
//     matrix[2,1] = 6.0
//     tmp.AddMatrix(key, matrix)
//     A = tmp[key].GetMatrix()
//     KRATOS_CHECK_EQUAL(A[0,0],1.0)
//     KRATOS_CHECK_EQUAL(A[0,1],2.0)
//     KRATOS_CHECK_EQUAL(A[1,0],3.0)
//     KRATOS_CHECK_EQUAL(A[1,1],4.0)
//     KRATOS_CHECK_EQUAL(A[2,0],5.0)
//     KRATOS_CHECK_EQUAL(A[2,1],6.0)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersVectorInterface, KratosCoreFastSuite)
// {
//     // Read and check Vectors from a Parameters-Object
//     Parameters tmp = Parameters(R"({
//         "valid_vectors" : [ []
//         ],
//         "false_vectors" : [ [[]],
//                             [[2,3],2],
//                             [2,3,[2]],
//                             [2,3,[]],
//                             [{"key":3},2],
//                             [2,3,{"key":3}],
//                             [true,2],
//                             [2,3,true],
//                             [5,"string",2]
//         ]
//     })");
//
//     // Check the IsVector Method
//     for i in range(tmp["valid_vectors"].size()):
//         valid_vector = tmp["valid_vectors"][i]
//         KRATOS_CHECK(valid_vector.IsVector())
//
//     for i in range(tmp["false_vectors"].size()):
//         false_vector = tmp["false_vectors"][i]
//         KRATOS_CHECK_IS_FALSE(false_vector.IsVector())
//
//     // Check the GetVector Method also on the valid Matrices
//     for i in range(tmp["valid_vectors"].size()):
//         valid_vector = tmp["valid_vectors"][i]
//         valid_vector.GetVector()
//
//     // Check that the errors of the GetVector method are thrown correctly
//     for i in range(tmp["false_vectors"].size()):
//         false_vector = tmp["false_vectors"][i]
//         with self.assertRaises(RuntimeError):
//             false_vector.GetVector()
//
//     // Manually assign and check a Vector
//     vec = Vector(3)
//     vec[0] = 1.32
//     vec[1] = -2.22
//     vec[2] = 5.5
//
//     tmp.AddEmptyValue("vector_value")
//     tmp["vector_value"].SetVector(vec)
//
//     KRATOS_CHECK(tmp["vector_value"].IsVector())
//
//     V2 = tmp["vector_value"].GetVector()
//     KRATOS_CHECK_EQUAL(V2[0],1.32)
//     KRATOS_CHECK_EQUAL(V2[1],-2.22)
//     KRATOS_CHECK_EQUAL(V2[2],5.5)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersMatrixInterface, KratosCoreFastSuite)
// {
//     // Read and check Matrices from a Parameters-Object
//     Parameters tmp = Parameters(R"({
//         "valid_matrices" : [ [[]],
//                               [[],[]],
//                               [[-9.81,8, 5.47]]
//         ],
//         "false_matrices" : [ [],
//                               [[[]]],
//                               [[3.3] , [1,2]],
//                               [[2,1.5,3.3] , [3,{"key":3},2]],
//                               [[2,1.5,3.3] , [5,false,2]],
//                               [[2,1.5,3.3] , [[2,3],1,2]],
//                               [[2,1.5,3.3] , ["string",2,9]]
//         ]
//     })");
//
//     // Check the IsMatrix Method
//     for i in range(tmp["valid_matrices"].size()):
//         valid_matrix = tmp["valid_matrices"][i]
//         KRATOS_CHECK(valid_matrix.IsMatrix())
//
//     for i in range(tmp["false_matrices"].size()):
//         false_matrix = tmp["false_matrices"][i]
//         KRATOS_CHECK_IS_FALSE(false_matrix.IsMatrix())
//
//     // Check the GetMatrix Method also on the valid Matrices
//     for i in range(tmp["valid_matrices"].size()):
//         valid_matrix = tmp["valid_matrices"][i]
//         valid_matrix.GetMatrix()
//
//     // Check that the errors of the GetMatrix method are thrown correctly
//     for i in range(tmp["false_matrices"].size()):
//         false_matrix = tmp["false_matrices"][i]
//         with self.assertRaises(RuntimeError):
//             false_matrix.GetMatrix()
//
//     // Manually assign and check a Matrix
//     mat = Matrix(3,2)
//     mat[0,0] = 1.0
//     mat[0,1] = 2.0
//     mat[1,0] = 3.0
//     mat[1,1] = 4.0
//     mat[2,0] = 5.0
//     mat[2,1] = 6.0
//
//     tmp.AddEmptyValue("matrix_value")
//     tmp["matrix_value"].SetMatrix(mat)
//
//     KRATOS_CHECK(tmp["matrix_value"].IsMatrix())
//
//     A2 = tmp["matrix_value"].GetMatrix()
//     KRATOS_CHECK_EQUAL(A2[0,0],1.0);
//     KRATOS_CHECK_EQUAL(A2[0,1],2.0);
//     KRATOS_CHECK_EQUAL(A2[1,0],3.0);
//     KRATOS_CHECK_EQUAL(A2[1,1],4.0);
//     KRATOS_CHECK_EQUAL(A2[2,0],5.0);
//     KRATOS_CHECK_EQUAL(A2[2,1],6.0);
// }
//
// def test_null_vs_null_validation(self):
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersNullvsNullValidation, KratosCoreFastSuite)
// {
//     // Supplied settings
//     Parameters null_custom = Parameters(R"({
//         "parameter": null
//     })");
//
//     // Default settings
//     Parameters null_default = Parameters(R"({
//         "parameter": null
//     })");
//
//     // This should NOT raise, hence making the test to pass
//     null_custom.ValidateAndAssignDefaults(null_default)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersDoublevsNullValidation, KratosCoreFastSuite)
// {
//     // Supplied settings
//     Parameters double_custom = Parameters(R"({
//         "parameter": 0.0
//     })");
//
//     // Default settings
//     Parameters null_default = Parameters(R"({
//         "parameter": null
//     })");
//
//     with self.assertRaises(RuntimeError):
//         double_custom.ValidateAndAssignDefaults(null_default)
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersGeStringArrayValid, KratosCoreFastSuite)
// {
//     Parameters tmp = Parameters(R"({
//         "parameter": ["foo", "bar"]
//     })");
//     v = tmp["parameter"].GetStringArray()
//     KRATOS_CHECK_EQUAL(len(v), 2);
//     KRATOS_CHECK_EQUAL(v[0], "foo");
//     KRATOS_CHECK_EQUAL(v[1], "bar");
// }
//
// KRATOS_TEST_CASE_IN_SUITE(KratosParametersSetStringArrayValid, KratosCoreFastSuite)
// {
//     Parameters initial = Parameters(R"({
//         "parameter": ["foo", "bar"]
//     })");
//     string_array = initial["parameter"].GetStringArray();
//
//     Parameters new_param = Parameters();
//     new_param.AddEmptyValue("new_parameter");
//     new_param["new_parameter"].SetStringArray(string_array);
//
//     new_string_array = initial["parameter"].GetStringArray();
//
//     self.assertListEqual(new_string_array, string_array);
// }


}  // namespace Testing.
}  // namespace Kratos.
