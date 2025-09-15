# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compute_pressure_coefficient_process


class ComputePressureCoefficientProcessTest(KratosUnittest.TestCase):
    @classmethod
    def _CreateModel(self):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("main")

        mpart.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        mpart.CreateNewNode(1, 0.0, 1.0, 0.0)
        mpart.CreateNewNode(2, 2.0, 0.0, 0.0)
        mpart.CreateNewCondition("LineCondition2D2N", 1, [1, 2], mpart.GetProperties()[0])

        return model

    @classmethod
    def _GetBlankParameters(cls):
        return KratosMultiphysics.Parameters("""
        {
            "python_module" : "compute_pressure_coefficient_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name" : "ComputePressureCoefficientProcess",
            "Parameters"    : {
                "model_part_name" : "main"
            }
        }
        """)

    def testPressureCoefficientOutputOnly(self):
        model = self._CreateModel()
        mpart = model["main"]
        mpart.GetNode(1).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e6)
        mpart.GetNode(2).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e4)

        params = self._GetBlankParameters()
        params["Parameters"].AddEmptyValue("freestream_pressure").SetDouble(1e5)
        params["Parameters"].AddEmptyValue("freestream_velocity").SetDouble(100)
        params["Parameters"].AddEmptyValue("freestream_density").SetDouble(1.2)
        params["Parameters"].AddEmptyValue("execution_step").SetString("ExecuteBeforeOutputStep")

        process = compute_pressure_coefficient_process.Factory(params, model)
        process.ExecuteInitialize()

        for node in mpart.Nodes:
            self.assertTrue(node.Has(KratosMultiphysics.PRESSURE_COEFFICIENT))
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 0.0)

        process.ExecuteFinalizeSolutionStep()

        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 0.0)
        self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 0.0)

        process.ExecuteBeforeOutputStep()

        self.assertAlmostEqual(mpart.GetNode(1).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 150.0)
        self.assertAlmostEqual(mpart.GetNode(2).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), -15.0)

    def testPressureCoefficientEveryStep(self):
        model = self._CreateModel()
        mpart = model["main"]
        mpart.GetNode(1).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e6)
        mpart.GetNode(2).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e4)

        params = self._GetBlankParameters()
        params["Parameters"].AddEmptyValue("freestream_pressure").SetDouble(1e5)
        params["Parameters"].AddEmptyValue("freestream_velocity").SetDouble(100)
        params["Parameters"].AddEmptyValue("freestream_density").SetDouble(1.2)
        params["Parameters"].AddEmptyValue("execution_step").SetString("ExecuteFinalizeSolutionStep")

        process = compute_pressure_coefficient_process.Factory(params, model)
        process.ExecuteInitialize()

        for node in mpart.Nodes:
            self.assertTrue(node.Has(KratosMultiphysics.PRESSURE_COEFFICIENT))
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 0.0)

        process.ExecuteFinalizeSolutionStep()

        self.assertAlmostEqual(mpart.GetNode(1).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 150.0)
        self.assertAlmostEqual(mpart.GetNode(2).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), -15.0)

        process.ExecuteBeforeOutputStep()

        self.assertAlmostEqual(mpart.GetNode(1).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 150.0)
        self.assertAlmostEqual(mpart.GetNode(2).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), -15.0)


if __name__ == '__main__':
    KratosUnittest.main()
