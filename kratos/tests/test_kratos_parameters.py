from __future__ import print_function, absolute_import, division
from KratosMultiphysics import Parameters
from KratosMultiphysics import Vector
from KratosMultiphysics import Matrix

import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys


# input string with ugly formatting
json_string = """
{
   "int_value" : 10,   "double_value": 2.0,   "bool_value" : true,   "string_value" : "hello",
   "level1":
   {
     "list_value":[ 3, "hi", false],
     "tmp" : 5.0
   }
}
"""

pretty_out = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": {
        "list_value": [
            3,
            "hi",
            false
        ],
        "tmp": 5.0
    }
}"""

pretty_out_after_change = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": {
        "list_value": [
            "changed",
            "hi",
            false
        ],
        "tmp": 5.0
    }
}"""

# here the level1 var is set to a double so that a validation error should be thrown
wrong_type = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": 0.0
}"""

# int value is badly spelt
wrong_spelling = """{
    "int_values": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": 0.0
}"""

# wrong on the first level
# error shall be only detective by recursive validation
wrong_lev2 = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": { "a":0.0 }
}"""

defaults = """
{
	"int_value": 10,
	"double_value": 2.0,
	"bool_value": false,
	"string_value": "hello",
	"level1": {
		"list_value": [
			3,
			"hi",
			false
		],
		"tmp": "here we expect a string"
	},
	"new_default_value": -123.0,
	"new_default_obj": {
		"aaa": "string",
		"bbb": false,
		"ccc": 22
	}
}
"""


expected_validation_output = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": {
        "list_value": [
            3,
            "hi",
            false
        ],
        "tmp": 5.0
    },
    "new_default_value": -123.0,
    "new_default_obj": {
        "aaa": "string",
        "bbb": false,
        "ccc": 22
    }
}"""

four_levels = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": { 
        "level2": {
            "level3": {
                "level4": {
                }
            }
        } 
    }
}"""

four_levels_variation = """{
    "int_value": 10,
    "double_value": 2.0,
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
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
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
    }
}"""

class TestParameters(KratosUnittest.TestCase):    

    def setUp(self):
        self.kp = Parameters(json_string)
        self.compact_expected_output = """{"int_value":10,"double_value":2.0,"bool_value":true,"string_value":"hello","level1":{"list_value":[3,"hi",false],"tmp":5.0}}"""
        
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
                          

        expected_values = ['10', '2.0', 'true', '"hello"', '{"list_value":[3,"hi",false],"tmp":5.0}']
        counter = 0
        
        for value in kp.values():
            self.assertEqual(value.WriteJsonString(), expected_values[counter])
            counter += 1

        #testing values
        expected_keys = ['int_value', 'double_value', 'bool_value', 'string_value', 'level1'] 
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
        
    def test_vector_interface(self):
        # Read and check a Vector from a Parameters-Object
        tmp = Parameters("""{
            "vector_value": [1.32,-2.22,5.5],
            "false_vector_value_1" : 5,
            "false_vector_value_2" : 7.2,
            "false_vector_value_3": [[2,1.5,3.3,6],[1,2,7.9,6]]
        }""")
        
        self.assertTrue(tmp["vector_value"].IsVector())
        
        V = tmp["vector_value"].GetVector()
        self.assertEqual(V[0],1.32)
        self.assertEqual(V[1],-2.22)
        self.assertEqual(V[2],5.5)
        
        # Manually assign and check a Vector
        vec = Vector(3)
        vec[0] = 1.32
        vec[1] = -2.22
        vec[2] = 5.5

        tmp.AddEmptyValue("vector_value2")
        tmp["vector_value2"].SetVector(vec)

        self.assertTrue(tmp["vector_value2"].IsVector())

        V2 = tmp["vector_value2"].GetVector()
        self.assertEqual(V2[0],1.32)
        self.assertEqual(V2[1],-2.22)
        self.assertEqual(V2[2],5.5)

        # Check the IsVector method
        self.assertFalse(tmp["false_vector_value_1"].IsVector()) # int
        self.assertFalse(tmp["false_vector_value_2"].IsVector()) # double
        self.assertFalse(tmp["false_vector_value_3"].IsVector()) # Matrix

        # check that the errors of the GetVector method are thrown correctly
        with self.assertRaises(RuntimeError):
            tmp["false_vector_value_1"].GetVector() # int
        with self.assertRaises(RuntimeError):
            tmp["false_vector_value_2"].GetVector() # double
        with self.assertRaises(RuntimeError):
            tmp["false_vector_value_3"].GetVector() # Matrix

    def test_matrix_interface(self):
        # Read and check a Matrix from a Parameters-Object
        tmp = Parameters("""{
            "matrix_value": [[1,2],[3,4],[5,6]],
            "false_matrix_value_1": 2,
            "false_matrix_value_2": 2.9,
            "false_matrix_value_3": [2, 1.5],
            "false_matrix_value_4": [[2, 1.5,3.3],[1,2]]
        }""")
        
        self.assertTrue(tmp["matrix_value"].IsMatrix())
        
        A = tmp["matrix_value"].GetMatrix()
        self.assertEqual(A[0,0],1.0)
        self.assertEqual(A[0,1],2.0)
        self.assertEqual(A[1,0],3.0)
        self.assertEqual(A[1,1],4.0)
        self.assertEqual(A[2,0],5.0)
        self.assertEqual(A[2,1],6.0)
        
        # Manually assign and check a Matrix
        mat = Matrix(3,2)
        mat[0,0] = 1.0
        mat[0,1] = 2.0
        mat[1,0] = 3.0
        mat[1,1] = 4.0
        mat[2,0] = 5.0
        mat[2,1] = 6.0

        tmp.AddEmptyValue("matrix_value2")
        tmp["matrix_value2"].SetMatrix(mat)
        
        self.assertTrue(tmp["matrix_value2"].IsMatrix())
        
        A2 = tmp["matrix_value2"].GetMatrix()
        self.assertEqual(A2[0,0],1.0)
        self.assertEqual(A2[0,1],2.0)
        self.assertEqual(A2[1,0],3.0)
        self.assertEqual(A2[1,1],4.0)
        self.assertEqual(A2[2,0],5.0)
        self.assertEqual(A2[2,1],6.0)

        # Check the IsMatrix Method
        self.assertFalse(tmp["false_matrix_value_1"].IsMatrix()) # int
        self.assertFalse(tmp["false_matrix_value_2"].IsMatrix()) # double
        self.assertFalse(tmp["false_matrix_value_3"].IsMatrix()) # Vector
        self.assertFalse(tmp["false_matrix_value_4"].IsMatrix()) # Mis-sized Matrix

        # check that the errors of the GetMatrix method are thrown correctly
        with self.assertRaises(RuntimeError):
            tmp["false_matrix_value_1"].GetMatrix() # int
        with self.assertRaises(RuntimeError):        
            tmp["false_matrix_value_2"].GetMatrix() # double
        with self.assertRaises(RuntimeError):        
            tmp["false_matrix_value_3"].GetMatrix() # Vector
        with self.assertRaises(RuntimeError):        
            tmp["false_matrix_value_4"].GetMatrix() # Mis-sized Matrix
        
        
if __name__ == '__main__':
    KratosUnittest.main()
