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
        self.compact_expected_output = """{
            "int_value":10,
            "double_value":2.0,
            "bool_value":true,
            "string_value":"hello",
            "level1":{
                "list_value":[3,"hi",false],
                "tmp":5.0
            }
        }"""

    def test_kratos_parameters(self):
        self.assertEqual(
            self.kp.WriteJsonString(),
            self.compact_expected_output
        )

        self.assertTrue(self.kp.Has("int_value"))
        self.assertFalse(self.kp.Has("unextisting_value"))

        self.assertEqual(self.kp.GetValue("int_value").GetInt(), 10)
        self.assertEqual(self.kp.GetValue("double_value").GetDouble(), 2.0)
        self.assertEqual(self.kp.GetValue("bool_value").GetBool(), True)
        self.assertEqual(self.kp.GetValue("string_value").GetString(), "hello")

        self.assertEqual(self.kp.PrettyPrintJsonString(), pretty_out)

    def test_kratos_change_parameters(self):
        # now change one item in the sublist
        subparams = self.kp.GetValue("level1")

        my_list = subparams.GetValue("list_value")
        my_list.GetArrayItem(0).SetString("changed")

        self.assertEqual(
            self.kp.PrettyPrintJsonString(),
            pretty_out_after_change
        )

    def test_kratos_copy_parameters(self):
        # try to make a copy
        other_copy = KratosParameters(self.kp)

        self.assertEqual(
            other_copy.PrettyPrintJsonString(),
            pretty_out_after_change
        )

        other_copy.GetValue("int_value").SetInt(-1)
        self.assertEqual(self.kp.GetValue("int_value").GetInt(), 10)
        # self.assertEqual(other_copy.GetValue("int_value").GetString(),-1)

    def test_kratos_wrong_parameters(self):
        # should check which errors are thrown!!
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            self.kp.GetValue("no_value").GetInt()

if __name__ == '__main__':
    KratosUnittest.main()
