# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compute_aerodynamic_coefficients_process

from KratosMultiphysics import DENSITY, PRESSURE, PRESSURE_COEFFICIENT


class ComputeAerodynamicCoefficientsProcessTest(KratosUnittest.TestCase):
    @classmethod
    def _CreateModel(self, *variables):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("main")

        for var in variables:
            mpart.AddNodalSolutionStepVariable(var)

        mpart.CreateNewNode(1, 0.0, 0.0, 0.0)
        mpart.CreateNewNode(2, 1.0, 0.0, 0.0)
        mpart.CreateNewCondition("LineCondition2D2N", 1, [1, 2], mpart.GetProperties()[0])

        return model

    def _GetBlankParameters(cls):
        return KratosMultiphysics.Parameters("""
        {
            "python_module" : "compute_aerodynamic_coefficients_process",
            "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
            "process_name" : "ComputeAerodynamicCoefficientsProcess",
            "Parameters"    : {
                "model_part_name" : "main"
            }
        }
        """)

    def testNonsenseDensity(self):
        model = self._CreateModel(DENSITY, PRESSURE)
        mpart = model["main"]
        params = self._GetBlankParameters()
        params["Parameters"].AddEmptyValue("density_database").SetString("nonsensical_string_hello!")

        with self.assertRaises(RuntimeError) as context:
            _ = compute_aerodynamic_coefficients_process.Factory(params, model)

        self.assertIn("Error: Invalid database.", str(context.exception))
        self.assertIn("- historical", str(context.exception))
        self.assertIn("- non-historical", str(context.exception))
        self.assertIn("- properties", str(context.exception))

    def testHistorical(self):
        model = self._CreateModel(DENSITY, PRESSURE)
        mpart = model["main"]
        params = self._GetBlankParameters()
        params["Parameters"].AddEmptyValue("density_database").SetString("historical")

        process = compute_aerodynamic_coefficients_process.Factory(params, model)
        process.ExecuteInitialize()

        for node in mpart.Nodes:
            self.assertTrue(node.Has(PRESSURE_COEFFICIENT))
            self.assertAlmostEqual(node.GetValue(PRESSURE_COEFFICIENT), 0.0)

        # WIP: This test is not yet complete.




if __name__ == '__main__':
    KratosUnittest.main()



