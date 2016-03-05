from __future__ import print_function, absolute_import, division

from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as KratosUnittest


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

#here the level1 var is set to a double so that a validation error should be thrown
wrong_type = """{
    "int_value": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": 0.0
}"""

#int value is badly spelt
wrong_spelling = """{
    "int_values": 10,
    "double_value": 2.0,
    "bool_value": true,
    "string_value": "hello",
    "level1": 0.0
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


class TestParameters(KratosUnittest.TestCase):

    def setUp(self):
        self.kp = Parameters(json_string)
        self.compact_expected_output = """{"int_value":10,"double_value":2.0,"bool_value":true,"string_value":"hello","level1":{"list_value":[3,"hi",false],"tmp":5.0}}"""


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
                
        #my_list = subparams["list_value"]
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
        kp1 =  Parameters(pretty_out_after_change)
        
        kp["bool_value"] = kp1["level1"]
        kp["bool_value"].PrettyPrintJsonString()
        self.assertEqual(kp["bool_value"].PrettyPrintJsonString() , kp1["level1"].PrettyPrintJsonString())
        
        
    def test_kratos_wrong_parameters(self):
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            self.kp["no_value"].GetInt()

    def test_validation_fails_due_to_wrong_type(self):
        kp = Parameters(wrong_type)
        defaults_params =  Parameters(defaults)
        
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            kp.ValidateAndAssignDefaults(defaults_params)

    def test_validation_fails_due_to_wrong_spelling(self):
        kp = Parameters(wrong_spelling)
        defaults_params =  Parameters(defaults)
        
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            kp.ValidateAndAssignDefaults(defaults_params)
            

     
    def test_validation_succeeds(self):
        kp = Parameters(json_string)
        defaults_params =  Parameters(defaults)
        defaults_params["level1"]["tmp"].SetDouble(2.0) ##this does not coincide with the value in kp, but is of the same type 
        
        kp.ValidateAndAssignDefaults(defaults_params)
        self.assertEqual(kp.PrettyPrintJsonString() , expected_validation_output)

        self.assertEqual(kp["level1"]["tmp"].GetDouble(), 5.0) ##not 2, since kp overwrites the defaults
     
            
if __name__ == '__main__':
    KratosUnittest.main()
