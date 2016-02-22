import unittest
from KratosMultiphysics import *

##input string with ugly formatting
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

pretty_out_after_change ="""{
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


class TestKratosParameters(unittest.TestCase):

    def test_kratos_parameters(self):
        
        kp = KratosParameters(json_string)
        
        compact_expected_output = """{"int_value":10,"double_value":2.0,"bool_value":true,"string_value":"hello","level1":{"list_value":[3,"hi",false],"tmp":5.0}}"""
        self.assertEqual(kp.WriteJsonString(), compact_expected_output)
            
        self.assertTrue(kp.Has("int_value"))
        self.assertFalse(kp.Has("unextisting_value"))
        
        self.assertEqual(kp.GetValue("int_value").GetInt(),10)
        self.assertEqual(kp.GetValue("double_value").GetDouble(),2.0)
        self.assertEqual(kp.GetValue("bool_value").GetBool(),True)
        self.assertEqual(kp.GetValue("string_value").GetString(),"hello")
        
        self.assertEqual(kp.PrettyPrintJsonString(), pretty_out)
        
        #now change one item in the sublist
        subparams = kp.GetValue("level1")
        my_list = subparams.GetValue("list_value")
        my_list.GetArrayItem(0).SetString("changed")
        self.assertEqual(kp.PrettyPrintJsonString(), pretty_out_after_change)
        
        #try to make a copy
        other_copy =  KratosParameters(kp)
        self.assertEqual(other_copy.PrettyPrintJsonString(), pretty_out_after_change)
        other_copy.GetValue("int_value").SetInt(-1)
        self.assertEqual(kp.GetValue("int_value").GetInt(),10)
        #self.assertEqual(other_copy.GetValue("int_value").GetString(),-1)
        
        #should check which errors are thrown!!

if __name__ == '__main__':
    unittest.main()