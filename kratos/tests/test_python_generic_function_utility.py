import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPythonGenericFunctionUtility(KratosUnittest.TestCase):

    def test_PythonGenericFunctionUtility(self):
        settings = KratosMultiphysics.Parameters("""
        {
            "local_axes"      : {}
        }
        """
        )

        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")
        node = model_part.CreateNewNode(1, 1.00,0.00,0.00)
        node.X += 0.1
        current_time = model_part.ProcessInfo[KratosMultiphysics.TIME]

        aux_function = KratosMultiphysics.PythonGenericFunctionUtility("x", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.1)

        aux_function = KratosMultiphysics.PythonGenericFunctionUtility("X", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.0)

        aux_function = KratosMultiphysics.PythonGenericFunctionUtility("X+Y", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 1.0)

        aux_function = KratosMultiphysics.PythonGenericFunctionUtility("x+X", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 2.1)

        aux_function = KratosMultiphysics.PythonGenericFunctionUtility("t", settings["local_axes"])
        value = aux_function.CallFunction(node.X , node.Y , node.Z, current_time, node.X0, node.Y0, node.Z0)

        self.assertEqual(value, 0.0)

if __name__ == '__main__':
    KratosUnittest.main()
