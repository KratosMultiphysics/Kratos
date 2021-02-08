import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

import math
import os
import sys

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestPythonGenericFunctionUtility(KratosUnittest.TestCase):

    def test_PythonGenericFunctionUtility1(self):
        function1 = KM.PythonGenericFunctionUtility("x**2+y**2")

        self.assertTrue(function1.DependsOnSpace())
        self.assertFalse(function1.UseLocalSystem())
        self.assertEqual(function1.FunctionBody(), "x**2+y**2")
        self.assertEqual(function1.CallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 25)

        function2 = KM.PythonGenericFunctionUtility("3*t")

        self.assertFalse(function2.DependsOnSpace())
        self.assertFalse(function2.UseLocalSystem())
        self.assertEqual(function2.FunctionBody(), "3*t")
        self.assertEqual(function2.CallFunction(0.0,0.0,0.0,5.0,0.0,0.0,0.0), 15)

        function3 = KM.PythonGenericFunctionUtility("X**2+Y**2")

        self.assertTrue(function3.DependsOnSpace())
        self.assertFalse(function3.UseLocalSystem())
        self.assertEqual(function3.FunctionBody(), "X**2+Y**2")
        self.assertEqual(function3.CallFunction(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25)

        function4 = KM.PythonGenericFunctionUtility("(cos(x)+sin(y))*t")

        self.assertTrue(function4.DependsOnSpace())
        self.assertFalse(function4.UseLocalSystem())
        self.assertEqual(function4.FunctionBody(), "(cos(x)+sin(y))*t")
        self.assertEqual(function4.CallFunction(0.25,0.15,0.0,1.5,0.0,0.0,0.0), 1.5*(math.cos(0.25) + math.sin(0.15)))

    def test_PythonGenericFunctionUtility2(self):
        parameters = KM.Parameters ("""{
            "origin" : [0,0,0],
            "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

        }""")

        function = KM.PythonGenericFunctionUtility("x+2*y", parameters)

        self.assertTrue(function.DependsOnSpace())
        self.assertTrue(function.UseLocalSystem())
        self.assertEqual(function.FunctionBody(), "x+2*y")
        self.assertEqual(function.CallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 10)
        self.assertEqual(function.RotateAndCallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 11)

    def test_ApplyFunctionToNodesUtility(self):
        parameters = KM.Parameters ("""{
            "origin" : [0,0,0],
            "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

        }""")

        function = KM.PythonGenericFunctionUtility("x+2*y", parameters)

        self.assertTrue(function.DependsOnSpace())
        self.assertTrue(function.UseLocalSystem())
        self.assertEqual(function.FunctionBody(), "x+2*y")
        self.assertEqual(function.CallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 10)
        self.assertEqual(function.RotateAndCallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 11)

        this_model = KM.Model()
        model_part = this_model.CreateModelPart("Main", 2)
        current_process_info = model_part.ProcessInfo
        current_process_info[KM.DOMAIN_SIZE] = 2

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)

        model_part_io = KM.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        utility = KM.ApplyFunctionToNodesUtility(model_part.Nodes, function)
        utility.ApplyFunction(KM.VISCOSITY, 1.0)

        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KM.VISCOSITY) - (node.Y + 2.0 * node.X), 0.0)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
