import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

import math
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestGenericFunctionUtility(KratosUnittest.TestCase):

    def test_GenericFunctionUtility0(self):
        settings = KM.Parameters("""
        {
            "local_axes"      : {}
        }
        """
        )

        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")
        node = model_part.CreateNewNode(1, 1.00,0.00,0.00)
        node.X += 0.1
        current_time = model_part.ProcessInfo[KM.TIME]

        aux_function = KM.GenericFunctionUtility("x", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.1)

        aux_function = KM.GenericFunctionUtility("X", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.0)

        aux_function = KM.GenericFunctionUtility("X+Y", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.0)

        aux_function = KM.GenericFunctionUtility("x+X", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 2.1)

        aux_function = KM.GenericFunctionUtility("t", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 0.0)

    def test_GenericFunctionUtility1(self):
        function1 = KM.GenericFunctionUtility("x**2+y**2")

        self.assertTrue(function1.DependsOnSpace())
        self.assertFalse(function1.UseLocalSystem())
        self.assertEqual(function1.FunctionBody(), "x**2+y**2")
        self.assertEqual(function1.CallFunction(4.0,3.0,0.0,0.0,0.0,0.0,0.0), 25)

        function2 = KM.GenericFunctionUtility("3*t")

        self.assertFalse(function2.DependsOnSpace())
        self.assertFalse(function2.UseLocalSystem())
        self.assertEqual(function2.FunctionBody(), "3*t")
        self.assertEqual(function2.CallFunction(0.0,0.0,0.0,5.0,0.0,0.0,0.0), 15)

        function3 = KM.GenericFunctionUtility("X**2+Y**2")

        self.assertTrue(function3.DependsOnSpace())
        self.assertFalse(function3.UseLocalSystem())
        self.assertEqual(function3.FunctionBody(), "X**2+Y**2")
        self.assertEqual(function3.CallFunction(0.0,0.0,0.0,0.0,4.0,3.0,0.0), 25)

        function4 = KM.GenericFunctionUtility("(cos(x*pi)+sin(y*pi))*t")

        self.assertTrue(function4.DependsOnSpace())
        self.assertFalse(function4.UseLocalSystem())
        self.assertEqual(function4.FunctionBody(), "(cos(x*pi)+sin(y*pi))*t")
        self.assertEqual(function4.CallFunction(0.25,0.15,0.0,1.5,0.0,0.0,0.0), 1.5*(math.cos(0.25*math.pi) + math.sin(0.15*math.pi)))

    def test_GenericFunctionUtility2(self):
        parameters = KM.Parameters ("""{
            "origin" : [0,0,0],
            "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

        }""")

        function = KM.GenericFunctionUtility("x+2*y", parameters)

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

        function = KM.GenericFunctionUtility("x+2*y", parameters)

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
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)

        model_part_io = KM.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        utility = KM.ApplyFunctionToNodesUtility(model_part.Nodes, function)
        utility.ApplyFunction(KM.VISCOSITY, 1.0)

        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(KM.VISCOSITY) - (node.Y + 2.0 * node.X), 0.0)

    def test_ApplyFunctionToNodesUtilityTimeEvolutionTernary(self):
        parameters = KM.Parameters ("""{}""")
        function = KM.GenericFunctionUtility("1.5*(0.5*(1-cos(0.5*pi*t))*2.0)*(4.0/0.1681)*y*(0.41-y) if t<2.0 else 1.5*(2.0)*(4.0/0.1681)*y*(0.41-y)", parameters)

        this_model = KM.Model()
        model_part = this_model.CreateModelPart("Main", 2)
        current_process_info = model_part.ProcessInfo
        current_process_info[KM.DOMAIN_SIZE] = 2

        model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)
        model_part.AddNodalSolutionStepVariable(KM.VELOCITY)

        model_part_io = KM.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read"))
        model_part_io.ReadModelPart(model_part)

        current_process_info[KM.TIME] = 0.0
        time = current_process_info[KM.TIME]
        while time < 3.0:
            current_process_info[KM.TIME] = current_process_info[KM.TIME] + 1.0
            time = current_process_info[KM.TIME]

            utility = KM.ApplyFunctionToNodesUtility(model_part.Nodes, function)
            utility.ApplyFunction(KM.VISCOSITY, time)

            if time < 2.0:
                for node in model_part.Nodes:
                    self.assertEqual(node.GetSolutionStepValue(KM.VISCOSITY) - (1.5*(0.5*(1-math.cos(0.5*math.pi*time))*2.0)*(4.0/0.1681)*node.Y*(0.41-node.Y)), 0.0)
            else:
                for node in model_part.Nodes:
                    self.assertAlmostEqual(node.GetSolutionStepValue(KM.VISCOSITY), 1.5*(2.0)*(4.0/0.1681)*node.Y*(0.41-node.Y))

    def test_ApplyFunctionToNodesUtilityTimeEvolutionCTernaryFail(self):
        parameters = KM.Parameters ("""{
            "origin" : [0,0,0],
            "axes"   : [[0,1,0],[1,0,0],[0,0,1]]

        }""")

        with self.assertRaisesRegex(Exception, 'Parsing error in function: 1.5 if t<2.0 3.0 if defined, but not else'):
            KM.GenericFunctionUtility("1.5 if t<2.0 3.0", parameters)

    def test_GenericFunctionUtilityError(self):
        with self.assertRaisesRegex(Exception, 'Parsing error in function: \(0\)\*\(50\*\(expp\(t\)-1\)\)\nError occurred near here :             \^ \(char \[12\]\)\nCheck your locale \(e.g. if "." or "," is used as decimal point'):
            KM.GenericFunctionUtility("(0)*(50*(expp(t)-1))")


if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
