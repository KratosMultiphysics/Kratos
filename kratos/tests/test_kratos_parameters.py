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


class TestKratosParameters(KratosUnittest.TestCase):

    def setUp(self):
        self.kp = KratosParameters(json_string)
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
        other_copy = KratosParameters(self.kp)      

        self.assertEqual(
            other_copy.PrettyPrintJsonString(),
            original_out
        )

        other_copy["int_value"].SetInt(-1)
        self.assertEqual(self.kp["int_value"].GetInt(), 10)
        # self.assertEqual(other_copy["int_value").GetString(),-1)

    def test_set_value(self):
        kp = KratosParameters(json_string)
        kp1 =  KratosParameters(pretty_out_after_change)
        
        kp["bool_value"] = kp1["level1"]
        kp["bool_value"].PrettyPrintJsonString()
        self.assertEqual(kp["bool_value"].PrettyPrintJsonString() , kp1["level1"].PrettyPrintJsonString())
        
        
    def test_kratos_wrong_parameters(self):
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            self.kp["no_value"].GetInt()


if __name__ == '__main__':
    KratosUnittest.main()
