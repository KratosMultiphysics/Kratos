from KratosMultiphysics import Parameters
from KratosMultiphysics import Vector
from KratosMultiphysics import Matrix
from KratosMultiphysics import FileSerializer, StreamSerializer, SerializerTraceType

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

# Use cPickle on Python 2.7 (Note that only the cPickle module is supported on Python 2.7)
# Source: https://pybind11.readthedocs.io/en/stable/advanced/classes.html
pickle_message = ""
try:
    import cPickle as pickle
    have_pickle_module = True
except ImportError:
    try:
        import pickle
        have_pickle_module = True
    except ImportError:
        have_pickle_module = False
        pickle_message = "No pickle module found"

# Import copy
import copy

# input string with ugly formatting
json_string = """
{
   "bool_value" : true, "double_value": 2.0, "int_value" : 10,
   "level1":
   {
     "list_value":[ 3, "hi", false],
     "tmp" : 5.0
   },
   "string_value" : "hello"
}
"""

pretty_out = """{
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
}"""

pretty_out_after_change = """{
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
}"""

# here the level1 var is set to a double so that a validation error should be thrown
wrong_type = """{
    "bool_value": true,
    "double_value": 2.0,
    "int_value": 10,
    "level1": 0.0,
    "string_value": "hello"
}"""

# int value is badly spelt
wrong_spelling = """{
    "bool_value": true,
    "double_value": 2.0,
    "int_values": 10,
    "level1": 0.0,
    "string_value": "hello"
}"""

# wrong on the first level
# error shall be only detective by recursive validation
wrong_lev2 = """{
    "bool_value": true,
    "double_value": 2.0,
    "int_value": 10,
    "level1": { "a":0.0 },
    "string_value": "hello"
}"""

defaults = """
{
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
}
"""

incomplete = """
{
    "level1": {
    },
    "new_default_obj": {
        "aaa": "string",
        "bbb": false,
        "ccc": 22
    },
    "new_default_value": -123.0,
    "string_value": "hello"
}"""

incomplete_with_extra_parameter = """
{
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
}"""

expected_validation_output = """{
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
}"""

four_levels = """{
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
}"""

four_levels_variation = """{
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
}"""

four_levels_wrong_variation = """{
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
}"""

four_levels_defaults = """{
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
}"""

json_with_includes = """
{
   "bool_value" : true, "double_value": 2.0, "int_value" : 10,
   "@include_json" : "cpp_tests/auxiliar_files_for_cpp_unnitest/test_included_parameters.json"
}
"""

more_levels_json_with_includes = """
{
   "bool_value" : true, "double_value": 2.0, "int_value" : 10,
   "level1":
   {
     "@include_json" : "cpp_tests/auxiliar_files_for_cpp_unnitest/more_levels_test_included_parameters.json"
   },
   "string_value" : "hello"
}
"""

vector_json_with_includes = """
{
   "bool_value" : true, "double_value": 2.0, "int_value" : 10,
   "level1": [
    {
         "@include_json" : "cpp_tests/auxiliar_files_for_cpp_unnitest/more_levels_test_included_parameters.json",
         "vector": [
            {
                "@include_json" : "cpp_tests/auxiliar_files_for_cpp_unnitest/more_levels_test_included_parameters.json"
            }
         ]
    },
    {
        "hello": 0
    },
    {
         "@include_json" : "cpp_tests/auxiliar_files_for_cpp_unnitest/more_levels_test_included_parameters.json"
    }
   ],
   "string_value" : "hello"
}
"""

class TestParameters(KratosUnittest.TestCase):

    def setUp(self):
        self.kp = Parameters(json_string)
        self.compact_expected_output = """{"bool_value":true,"double_value":2.0,"int_value":10,"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"}"""

    def test_kratos_parameters(self):
        self.assertEqual(
            self.kp.WriteJsonString(),
            self.compact_expected_output
        )

        self.assertTrue(self.kp.Has("int_value"))
        self.assertFalse(self.kp.Has("unextisting_value"))

        self.assertEqual(self.kp["int_value"].GetInt(), 10)
        self.assertEqual(self.kp["double_value"].GetDouble(), 2.0)
        self.assertEqual(self.kp["bool_value"].GetBool(), True)
        self.assertEqual(self.kp["string_value"].GetString(), "hello")

        self.assertEqual(self.kp.PrettyPrintJsonString(), pretty_out)

    def test_kratos_include_parameters(self):
        param = Parameters(json_with_includes)
        self.assertEqual(
            param.WriteJsonString(),
            self.compact_expected_output
        )

    def test_kratos_more_levels_include_parameters(self):
        param = Parameters(more_levels_json_with_includes)
        self.assertEqual(
            param.WriteJsonString(),
            self.compact_expected_output
        )

    def test_vector_json_with_include_parameters(self):
        param = Parameters(vector_json_with_includes)
        self.assertEqual(
            param.WriteJsonString(),
            '{"bool_value":true,"double_value":2.0,"int_value":10,"level1":[{"list_value":[3,"hi",false],"tmp":5.0,"vector":[{"list_value":[3,"hi",false],"tmp":5.0}]},{"hello":0},{"list_value":[3,"hi",false],"tmp":5.0}],"string_value":"hello"}'
        )

    def test_kratos_change_parameters(self):
        # now change one item in the sublist
        subparams = self.kp["level1"]

        my_list = subparams["list_value"]

        for i in range(my_list.size()):
            if my_list[i].IsBool():
                self.assertEqual(my_list[i].GetBool(), False)

        # my_list = subparams["list_value"]
        subparams["list_value"][0].SetString("changed")

        self.assertEqual(
            self.kp.PrettyPrintJsonString(),
            pretty_out_after_change
        )

    def test_kratos_copy_parameters(self):
        # try to make a copy
        original_out = self.kp.PrettyPrintJsonString()
        other_copy = self.kp.Clone()

        self.assertEqual(
            other_copy.PrettyPrintJsonString(),
            original_out
        )

        other_copy["int_value"].SetInt(-1)
        self.assertEqual(self.kp["int_value"].GetInt(), 10)
        # self.assertEqual(other_copy["int_value").GetString(),-1)

    def test_set_value(self):
        kp = Parameters(json_string)
        kp1 = Parameters(pretty_out_after_change)

        kp["bool_value"] = kp1["level1"]
        kp["bool_value"].PrettyPrintJsonString()
        self.assertEqual(kp["bool_value"].PrettyPrintJsonString(), kp1["level1"].PrettyPrintJsonString())

    def test_kratos_wrong_parameters(self):
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "no_value"):
            self.kp["no_value"].GetInt()

    def test_validation_fails_due_to_wrong_type(self):
        kp = Parameters(wrong_type)
        defaults_params = Parameters(defaults)

        # should check which errors are thrown!!
        with self.assertRaises(RuntimeError):
            kp.ValidateAndAssignDefaults(defaults_params)

    def test_validation_fails_due_to_wrong_spelling(self):
        kp = Parameters(wrong_spelling)
        defaults_params = Parameters(defaults)

        # should check which errors are thrown!!
        with self.assertRaises(RuntimeError):
            kp.ValidateAndAssignDefaults(defaults_params)

    def test_recursive_validation_fails_error_on_first_level(self):
        kp = Parameters(wrong_lev2)
        defaults_params = Parameters(defaults)

        # should check which errors are thrown!!
        with self.assertRaises(RuntimeError):
            kp.RecursivelyValidateAndAssignDefaults(defaults_params)

    def test_recursive_validation_4_levels(self):
        kp = Parameters(four_levels)
        kp_variation = Parameters(four_levels_variation)
        kp_wrong_wariation = Parameters(four_levels_wrong_variation)
        defaults_params = Parameters(four_levels_defaults)

        kp.RecursivelyValidateAndAssignDefaults(defaults_params)
        kp_variation.RecursivelyValidateAndAssignDefaults(defaults_params)

        self.assertTrue( kp.IsEquivalentTo(defaults_params) )
        self.assertFalse( kp_variation.IsEquivalentTo(defaults_params) )

        self.assertTrue( kp.HasSameKeysAndTypeOfValuesAs(defaults_params) )
        self.assertTrue( kp_variation.HasSameKeysAndTypeOfValuesAs(defaults_params) )
        self.assertFalse( kp_wrong_wariation.HasSameKeysAndTypeOfValuesAs(defaults_params) )

    def test_validation_succeds_error_on_first_level(self):
        kp = Parameters(wrong_lev2)
        defaults_params = Parameters(defaults)

        # here no error shall be thrown since validation is only done on level0
        kp.ValidateAndAssignDefaults(defaults_params)

    def test_validation_succeeds(self):
        kp = Parameters(json_string)
        defaults_params = Parameters(defaults)
        defaults_params["level1"]["tmp"].SetDouble(2.0)  # this does not coincide with the value in kp, but is of the same type

        kp.ValidateAndAssignDefaults(defaults_params)
        self.assertEqual(kp.PrettyPrintJsonString(), expected_validation_output)

        self.assertEqual(kp["level1"]["tmp"].GetDouble(), 5.0)  # not 2, since kp overwrites the defaults

    def test_add_missing_parameters(self):
        # only missing parameters are added, no complaints if there already exist more than in the defaults
        kp = Parameters(json_string)
        tmp = Parameters(incomplete_with_extra_parameter)

        kp.AddMissingParameters(tmp)

        self.assertEqual(kp["new_default_obj"]["aaa"].GetString(), "string")
        self.assertEqual(kp["string_value"].GetString(), "hello")
        self.assertFalse(kp["level1"].Has("new_sublevel"))

    def test_recursively_add_missing_parameters(self):
        # only missing parameters are added, no complaints if there already exist more than in the defaults
        kp = Parameters(json_string)
        tmp = Parameters(incomplete_with_extra_parameter)

        kp.RecursivelyAddMissingParameters(tmp)

        self.assertTrue(kp["level1"].Has("new_sublevel"))
        self.assertEqual(kp["level1"]["new_sublevel"].GetString(), "this should only be assigned in recursive")

    def test_validate_defaults(self):
        # only parameters from defaults are validated, no new values are added
        kp = Parameters(incomplete_with_extra_parameter)
        tmp = Parameters(defaults)

        kp.ValidateDefaults(tmp)

        self.assertFalse(kp.Has("bool_value"))
        self.assertFalse(kp.Has("double_value"))
        self.assertTrue(kp.Has("level1"))

    def test_recursively_validate_defaults(self):
        # only parameters from defaults are validated, no new values are added
        kp = Parameters(incomplete)
        tmp = Parameters(defaults)

        kp.RecursivelyValidateDefaults(tmp)

        self.assertFalse(kp.Has("bool_value"))
        self.assertFalse(kp.Has("double_value"))
        self.assertTrue(kp.Has("level1"))


    def test_recursively_validate_defaults_fails(self):
        # only parameters from defaults are validated, no new values are added
        kp = Parameters(incomplete_with_extra_parameter)
        tmp = Parameters(defaults)

        with self.assertRaises(RuntimeError):
            kp.RecursivelyValidateDefaults(tmp)

        # sub_level
        self.assertFalse(kp["level1"].Has("tmp"))

    def test_add_value(self):
        kp = Parameters("{}")
        kp.AddEmptyValue("new_double").SetDouble(1.0)

        self.assertTrue(kp.Has("new_double"))
        self.assertEqual(kp["new_double"].GetDouble(), 1.0)

    def test_add_empty_array(self):
        kp = Parameters("{}")
        kp.AddEmptyArray("new_array")

        self.assertTrue(kp.Has("new_array"))
        self.assertEqual(kp["new_array"].size(), 0)

    def test_iterators(self):
        kp = Parameters(json_string)

        #iteration by items
        for key,value in kp.items():
            #print(value.PrettyPrintJsonString())
            self.assertEqual(kp[key].PrettyPrintJsonString(), value.PrettyPrintJsonString())
            #print(key,value)

        #testing values
        expected_values = ['true', '2.0', '10', '{"list_value":[3,"hi",false],"tmp":5.0}','"hello"']
        counter = 0

        for value in kp.values():
            self.assertEqual(value.WriteJsonString(), expected_values[counter])
            counter += 1

        #testing values
        expected_keys = ['bool_value', 'double_value', 'int_value', 'level1', 'string_value']
        counter = 0
        for key in kp.keys():
            self.assertEqual(key, expected_keys[counter])
            counter += 1

    def test_remove_value(self):
        kp = Parameters(json_string)
        self.assertTrue(kp.Has("int_value"))
        self.assertTrue(kp.Has("level1"))

        kp.RemoveValue("int_value")
        kp.RemoveValue("level1")

        self.assertFalse(kp.Has("int_value"))
        self.assertFalse(kp.Has("level1"))

    def test_remove_values(self):
        kp = Parameters(json_string)
        self.assertTrue(kp.Has("int_value"))
        self.assertTrue(kp.Has("level1"))

        list_remove = ["int_value", "level1", "You_ll_never_take_me_alive"]
        success = kp.RemoveValues(list_remove)
        self.assertFalse(success)
        self.assertTrue(kp.Has("int_value"))
        self.assertTrue(kp.Has("level1"))

        list_remove = ["int_value", "level1"]
        kp.RemoveValues(list_remove)

        self.assertFalse(kp.Has("int_value"))
        self.assertFalse(kp.Has("level1"))

    def test_copy_deepcopy(self):
        kp = Parameters(json_string)

        # Copy values
        kp_copy1 = kp.__copy__()
        kp_copy2 = copy.copy(kp)
        kp_deepcopy1 = kp.__deepcopy__()
        kp_deepcopy2 = copy.deepcopy(kp)

        # Check is the same
        self.assertTrue(kp.Has("int_value"))
        self.assertTrue(kp.Has("level1"))
        self.assertTrue(kp_copy1.Has("int_value"))
        self.assertTrue(kp_copy1.Has("level1"))
        self.assertTrue(kp_copy2.Has("int_value"))
        self.assertTrue(kp_copy2.Has("level1"))
        self.assertTrue(kp_deepcopy1.Has("int_value"))
        self.assertTrue(kp_deepcopy1.Has("level1"))
        self.assertTrue(kp_deepcopy2.Has("int_value"))
        self.assertTrue(kp_deepcopy2.Has("level1"))

        # Remove values
        kp.RemoveValue("int_value")
        kp.RemoveValue("level1")

        # Check the deep copies is the same
        self.assertFalse(kp.Has("int_value"))
        self.assertFalse(kp.Has("level1"))
        self.assertFalse(kp_copy1.Has("int_value"))
        self.assertFalse(kp_copy1.Has("level1"))
        self.assertFalse(kp_copy2.Has("int_value"))
        self.assertFalse(kp_copy2.Has("level1"))
        self.assertTrue(kp_deepcopy1.Has("int_value"))
        self.assertTrue(kp_deepcopy1.Has("level1"))
        self.assertTrue(kp_deepcopy2.Has("int_value"))
        self.assertTrue(kp_deepcopy2.Has("level1"))

    def test_is_methods(self):
        # This method checks all the "IsXXX" Methods
        tmp = Parameters("""{
            "int_value" : 10, /* This is comment to check that comments work */
            "double_value": 2.0, // This is comment too, but using another comment
            "bool_value" : true, // This is another comment being meta as realizing that all the possibilities are already check
            "string_value" : "hello",/* This is a nihilist comment about the futile existence of the previous comment as a metacomment */
            "string_array_value" : ["hello", "world"],
            "vector_value" : [5,3,4],
            "matrix_value" : [[1,2],[3,6]]
        }""") # if you add more values to this, make sure to add the corresponding in the loop

        for key in tmp.keys():
            val_type = key[:-6] # removing "_value"

            if val_type == "int":
                self.assertTrue(tmp[key].IsInt())
            else:
                self.assertFalse(tmp[key].IsInt())

            if val_type == "double":
                self.assertTrue(tmp[key].IsDouble())
            else:
                self.assertFalse(tmp[key].IsDouble())

            if val_type == "bool":
                self.assertTrue(tmp[key].IsBool())
            else:
                self.assertFalse(tmp[key].IsBool())

            if val_type == "string":
                self.assertTrue(tmp[key].IsString())
            else:
                self.assertFalse(tmp[key].IsString())

            if val_type == "string_array":
                self.assertTrue(tmp[key].IsStringArray())
            else:
                self.assertFalse(tmp[key].IsStringArray())

            if val_type == "vector":
                self.assertTrue(tmp[key].IsVector())
            else:
                self.assertFalse(tmp[key].IsVector())

            if val_type == "matrix":
                self.assertTrue(tmp[key].IsMatrix())
            else:
                self.assertFalse(tmp[key].IsMatrix())

    def test_get_methods(self):
        # This method checks all the "GetXXX" Methods if they throw an error
        tmp = Parameters("""{
            "int_value" : 10,
            "double_value": 2.0,
            "bool_value" : true,
            "string_value" : "hello",
            "vector_value" : [5.2,-3.1,4.33],
            "matrix_value" : [[1,2],[3,4],[5,6]]
        }""") # if you add more values to this, make sure to add the corresponding in the loop

        for key in tmp.keys():
            val_type = key[:-6] # removing "_value"

            # Int and Double are checked tgth bcs both internally call "IsNumber"
            if val_type == "int" or val_type == "double":
                if val_type == "int":
                    self.assertEqual(tmp[key].GetInt(),10)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetInt()

            if val_type == "double" or val_type == "int":
                if val_type == "double":
                    self.assertEqual(tmp[key].GetDouble(),2.0)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetDouble()

            if val_type == "bool":
                self.assertEqual(tmp[key].GetBool(),True)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetBool()

            if val_type == "string":
                self.assertEqual(tmp[key].GetString(),"hello")
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetString()

            if val_type == "vector":
                V = tmp[key].GetVector()
                self.assertEqual(V[0],5.2)
                self.assertEqual(V[1],-3.1)
                self.assertEqual(V[2],4.33)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetVector()

            if val_type == "matrix":
                A = tmp[key].GetMatrix()
                self.assertEqual(A[0,0],1.0)
                self.assertEqual(A[0,1],2.0)
                self.assertEqual(A[1,0],3.0)
                self.assertEqual(A[1,1],4.0)
                self.assertEqual(A[2,0],5.0)
                self.assertEqual(A[2,1],6.0)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetMatrix()

    def test_set_methods(self):
        # This method checks all the "GetXXX" Methods if they throw an error
        tmp = Parameters("""{
            "int_value" : 0,
            "double_value": 0.0,
            "bool_value" : false,
            "string_value" : "",
            "vector_value" : [],
            "matrix_value" : [[0]]
        }""") # if you add more values to this, make sure to add the corresponding in the loop

        for key in tmp.keys():
            val_type = key[:-6] # removing "_value"

            # Int and Double are checked tgth bcs both internally call "IsNumber"
            if val_type == "int" or val_type == "double":
                if val_type == "int":
                    tmp[key].SetInt(10)
                    self.assertEqual(tmp[key].GetInt(),10)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetInt()

            if val_type == "double" or val_type == "int":
                if val_type == "double":
                    tmp[key].SetDouble(2.0)
                    self.assertEqual(tmp[key].GetDouble(),2.0)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetDouble()

            if val_type == "bool":
                tmp[key].SetBool(True)
                self.assertEqual(tmp[key].GetBool(),True)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetBool()

            if val_type == "string":
                tmp[key].SetString("hello")
                self.assertEqual(tmp[key].GetString(),"hello")
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetString()

            if val_type == "vector":
                vector = Vector(3)
                vector[0] = 5.2
                vector[1] = -3.1
                vector[2] = 4.33
                tmp[key].SetVector(vector)
                V = tmp[key].GetVector()
                self.assertEqual(V[0],5.2)
                self.assertEqual(V[1],-3.1)
                self.assertEqual(V[2],4.33)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetVector()

            if val_type == "matrix":
                matrix = Matrix(3,2)
                matrix[0,0] = 1.0
                matrix[0,1] = 2.0
                matrix[1,0] = 3.0
                matrix[1,1] = 4.0
                matrix[2,0] = 5.0
                matrix[2,1] = 6.0
                tmp[key].SetMatrix(matrix)
                A = tmp[key].GetMatrix()
                self.assertEqual(A[0,0],1.0)
                self.assertEqual(A[0,1],2.0)
                self.assertEqual(A[1,0],3.0)
                self.assertEqual(A[1,1],4.0)
                self.assertEqual(A[2,0],5.0)
                self.assertEqual(A[2,1],6.0)
            else:
                with self.assertRaises(RuntimeError):
                    tmp[key].GetMatrix()

    def test_add_methods(self):
        # This method checks all the "GetXXX" Methods if they throw an error
        tmp = Parameters("""{}""")

        key = "int"
        tmp.AddInt(key, 10)
        self.assertEqual(tmp[key].GetInt(),10)

        key = "double"
        tmp.AddDouble(key, 2.0)
        self.assertEqual(tmp[key].GetDouble(),2.0)

        key = "bool"
        tmp.AddBool(key, True)
        self.assertEqual(tmp[key].GetBool(),True)

        key = "string"
        tmp.AddString(key, "hello")
        self.assertEqual(tmp[key].GetString(),"hello")

        key = "vector"
        vector = Vector(3)
        vector[0] = 5.2
        vector[1] = -3.1
        vector[2] = 4.33
        tmp.AddVector(key, vector)
        V = tmp[key].GetVector()
        self.assertEqual(V[0],5.2)
        self.assertEqual(V[1],-3.1)
        self.assertEqual(V[2],4.33)

        key = "matrix"
        matrix = Matrix(3,2)
        matrix[0,0] = 1.0
        matrix[0,1] = 2.0
        matrix[1,0] = 3.0
        matrix[1,1] = 4.0
        matrix[2,0] = 5.0
        matrix[2,1] = 6.0
        tmp.AddMatrix(key, matrix)
        A = tmp[key].GetMatrix()
        self.assertEqual(A[0,0],1.0)
        self.assertEqual(A[0,1],2.0)
        self.assertEqual(A[1,0],3.0)
        self.assertEqual(A[1,1],4.0)
        self.assertEqual(A[2,0],5.0)
        self.assertEqual(A[2,1],6.0)


    def test_vector_interface(self):
        # Read and check Vectors from a Parameters-Object
        tmp = Parameters("""{
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
        }""")

        # Populate "valid_vectors" and test Parameters::Append
        tmp["valid_vectors"].Append(Vector([1, 2, 3]))
        tmp["valid_vectors"].Append([4, 5, 6])
        tmp["valid_vectors"].Append([0])
        tmp["valid_vectors"].Append([])

        # Check valid vectors
        for i in range(tmp["valid_vectors"].size()):
            valid_vector = tmp["valid_vectors"][i]
            self.assertTrue(valid_vector.IsVector())

        self.assertEqual(tmp["valid_vectors"][0].size(), 0)

        self.assertEqual(tmp["valid_vectors"][1].size(), 3)
        self.assertVectorAlmostEqual(tmp["valid_vectors"][1].GetVector(), [1, 2, 3])

        self.assertEqual(tmp["valid_vectors"][2].size(), 3)
        self.assertVectorAlmostEqual(tmp["valid_vectors"][2].GetVector(), [4, 5, 6])

        self.assertEqual(tmp["valid_vectors"][3].size(), 1)
        self.assertEqual(tmp["valid_vectors"][3].GetVector()[0], 0)

        # Check invalid vectors
        for i in range(tmp["false_vectors"].size()):
            false_vector = tmp["false_vectors"][i]
            self.assertFalse(false_vector.IsVector())

        # Check the GetVector Method also on the valid Matrices
        for i in range(tmp["valid_vectors"].size()):
            valid_vector = tmp["valid_vectors"][i]
            valid_vector.GetVector()

        # Check that the errors of the GetVector method are thrown correctly
        for i in range(tmp["false_vectors"].size()):
            false_vector = tmp["false_vectors"][i]
            with self.assertRaises(RuntimeError):
                false_vector.GetVector()

        # Manually assign and check a Vector
        vec = Vector(3)
        vec[0] = 1.32
        vec[1] = -2.22
        vec[2] = 5.5

        tmp.AddEmptyValue("vector_value")
        tmp["vector_value"].SetVector(vec)

        self.assertTrue(tmp["vector_value"].IsVector())

        V2 = tmp["vector_value"].GetVector()
        self.assertEqual(V2[0],1.32)
        self.assertEqual(V2[1],-2.22)
        self.assertEqual(V2[2],5.5)

    def test_matrix_interface(self):
        # Read and check Matrices from a Parameters-Object
        tmp = Parameters("""{
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
        }""")

        # Populate "valid_matrices" and test Parameters::Append
        tmp["valid_matrices"].Append(Matrix([[1, 2, 3], [4, 5, 6]]))
        tmp["valid_matrices"].Append(Matrix([[0]]))

        # Check valid matrices
        for i in range(tmp["valid_matrices"].size()):
            valid_matrix = tmp["valid_matrices"][i]
            self.assertTrue(valid_matrix.IsMatrix())

        # Check matrix values
        matrix = tmp["valid_matrices"][0].GetMatrix()
        self.assertEqual(matrix.Size1(), 1) # Number of rows in an empty matrix is non-zero! Is this intentional?
        self.assertEqual(matrix.Size2(), 0)

        matrix = tmp["valid_matrices"][3].GetMatrix()
        self.assertEqual(matrix.Size1(), 2)
        self.assertEqual(matrix.Size2(), 3)
        self.assertMatrixAlmostEqual(matrix, Matrix([[1, 2, 3], [4, 5, 6]]))

        matrix = tmp["valid_matrices"][4].GetMatrix()
        self.assertEqual(matrix.Size1(), 1)
        self.assertEqual(matrix.Size2(), 1)
        self.assertAlmostEqual(matrix[0, 0], 0)

        # Check invalid matrices
        for i in range(tmp["false_matrices"].size()):
            false_matrix = tmp["false_matrices"][i]
            self.assertFalse(false_matrix.IsMatrix())

        # Check the GetMatrix Method also on the valid Matrices
        for i in range(tmp["valid_matrices"].size()):
            valid_matrix = tmp["valid_matrices"][i]
            valid_matrix.GetMatrix()

        # Check that the errors of the GetMatrix method are thrown correctly
        for i in range(tmp["false_matrices"].size()):
            false_matrix = tmp["false_matrices"][i]
            with self.assertRaises(RuntimeError):
                false_matrix.GetMatrix()

        # Manually assign and check a Matrix
        mat = Matrix(3,2)
        mat[0,0] = 1.0
        mat[0,1] = 2.0
        mat[1,0] = 3.0
        mat[1,1] = 4.0
        mat[2,0] = 5.0
        mat[2,1] = 6.0

        tmp.AddEmptyValue("matrix_value")
        tmp["matrix_value"].SetMatrix(mat)

        self.assertTrue(tmp["matrix_value"].IsMatrix())

        A2 = tmp["matrix_value"].GetMatrix()
        self.assertEqual(A2[0,0],1.0)
        self.assertEqual(A2[0,1],2.0)
        self.assertEqual(A2[1,0],3.0)
        self.assertEqual(A2[1,1],4.0)
        self.assertEqual(A2[2,0],5.0)
        self.assertEqual(A2[2,1],6.0)

    def test_null_vs_null_validation(self):

        # supplied settings
        null_custom = Parameters("""{
        "parameter": null
        }""")

        # default settings
        null_default = Parameters("""{
        "parameter": null
        }""")

        #this should NOT raise, hence making the test to pass
        null_custom.ValidateAndAssignDefaults(null_default)

    def test_double_vs_null_validation(self):
        # supplied settings
        double_custom = Parameters("""{
        "parameter": 0.0
        }""")

        # default settings
        null_default = Parameters("""{
        "parameter": null
        }""")

        with self.assertRaises(RuntimeError):
            double_custom.ValidateAndAssignDefaults(null_default)

    def test_file_serialization(self):
        tmp = Parameters(defaults)
        check = tmp.WriteJsonString()

        file_name = "parameter_serialization"

        serializer = FileSerializer(file_name, SerializerTraceType.SERIALIZER_NO_TRACE)
        serializer.Save("ParametersSerialization",tmp)
        del(tmp)
        del(serializer)


        #unpickle data - note that here i override "serialized_data"
        serializer = FileSerializer(file_name,SerializerTraceType.SERIALIZER_NO_TRACE)

        loaded_parameters = Parameters()
        serializer.Load("ParametersSerialization",loaded_parameters)

        self.assertEqual(check, loaded_parameters.WriteJsonString())
        kratos_utils.DeleteFileIfExisting(file_name + ".rest")

    def test_get_string_array_valid(self):
        tmp = Parameters("""{
            "parameter": ["foo", "bar"]
        } """)
        v = tmp["parameter"].GetStringArray()
        self.assertEqual(len(v), 2)
        self.assertEqual(v[0], "foo")
        self.assertEqual(v[1], "bar")

    def test_get_string_array_invalid(self):
        tmp = Parameters("""{
            "parameter": ["foo", true]
        } """)
        with self.assertRaisesRegex(RuntimeError, r'Error: Argument must be a string'):
            tmp["parameter"].GetStringArray()

    def test_set_string_array_valid(self):
        initial = Parameters("""{
            "parameter": ["foo", "bar"]
        } """)
        string_array = initial["parameter"].GetStringArray()

        new_param = Parameters()
        new_param.AddEmptyValue("new_parameter")
        new_param["new_parameter"].SetStringArray(string_array)

        new_string_array = initial["parameter"].GetStringArray()

        self.assertListEqual(new_string_array, string_array)

    def test_add_string_array_valid(self):
        initial = Parameters("""{
            "parameter": ["foo", "bar"]
        } """)
        string_array = initial["parameter"].GetStringArray()

        new_param = Parameters()
        new_param.AddStringArray("new_parameter", string_array)

        new_string_array = initial["parameter"].GetStringArray()

        self.assertListEqual(new_string_array, string_array)

    def test_copy_values_from_existing_parameters(self):
        initial = Parameters("""{
            "parameter1": ["foo", "bar"],
            "parameter2": true,
            "parameter3": "Hello",
            "parameter4": 15
        } """)
        parameter_list = ["parameter2", "parameter3"]

        new_param = Parameters()
        new_param.CopyValuesFromExistingParameters(initial, parameter_list)

        self.assertFalse(new_param.Has("parameter1"))
        self.assertTrue(new_param.Has("parameter2"))
        self.assertEqual(new_param["parameter2"].GetBool(), True)
        self.assertTrue(new_param.Has("parameter3"))
        self.assertEqual(new_param["parameter3"].GetString(),"Hello")
        self.assertFalse(new_param.Has("parameter4"))

    @KratosUnittest.skipUnless(have_pickle_module, "Pickle module error: : " + pickle_message)
    def test_stream_serialization(self):
        tmp = Parameters(defaults)
        check = tmp.WriteJsonString()

        serializer = StreamSerializer(SerializerTraceType.SERIALIZER_NO_TRACE)
        serializer.Save("ParametersSerialization",tmp)
        del(tmp)

        #pickle dataserialized_data
        pickled_data = pickle.dumps(serializer, protocol=2) # Second argument is the protocol and is NECESSARY (according to pybind11 docs)
        del(serializer)

        #unpickle data - note that here i override "serialized_data"
        serializer = pickle.loads(pickled_data)

        loaded_parameters = Parameters()
        serializer.Load("ParametersSerialization",loaded_parameters)

        self.assertEqual(check, loaded_parameters.WriteJsonString())


if __name__ == '__main__':
    KratosUnittest.main()
