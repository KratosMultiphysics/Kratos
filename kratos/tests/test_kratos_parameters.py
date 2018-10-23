from __future__ import print_function, absolute_import, division
from KratosMultiphysics import Parameters
from KratosMultiphysics import Vector
from KratosMultiphysics import Matrix

import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys


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

class TestParameters(KratosUnittest.TestCase):

    def setUp(self):
        self.kp = Parameters(json_string)
        self.compact_expected_output = """{"bool_value":true,"double_value":2.0,"int_value":10,"level1":{"list_value":[3,"hi",false],"tmp":5.0},"string_value":"hello"}"""

        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

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

    def test_add_value(self):
        kp = Parameters("{}")
        kp.AddEmptyValue("new_double").SetDouble(1.0)

        self.assertTrue(kp.Has("new_double"))
        self.assertEqual(kp["new_double"].GetDouble(), 1.0)

    def test_iterators(self):
        kp = Parameters(json_string)

        #iteration by range
        nitems = 0
        for iterator in kp:
            nitems = nitems + 1
        self.assertEqual(nitems, 5)

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

    def test_is_methods(self):
        # This method checks all the "IsXXX" Methods
        tmp = Parameters("""{
            "int_value" : 10,
            "double_value": 2.0,
            "bool_value" : true,
            "string_value" : "hello",
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

        # Check the IsVector Method
        for i in range(tmp["valid_vectors"].size()):
            valid_vector = tmp["valid_vectors"][i]
            self.assertTrue(valid_vector.IsVector())

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

        # Check the IsMatrix Method
        for i in range(tmp["valid_matrices"].size()):
            valid_matrix = tmp["valid_matrices"][i]
            self.assertTrue(valid_matrix.IsMatrix())

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

if __name__ == '__main__':
    KratosUnittest.main()
