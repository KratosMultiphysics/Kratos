# Import kratos core and applications
import KratosMultiphysics
from KratosMultiphysics import FluidDynamicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.FluidDynamicsApplication import compute_aerodynamic_coefficients_process


class ComputeAerodynamicCoefficientsProcessTest(KratosUnittest.TestCase):
    @classmethod
    def _CreateModel(self):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("main")

        mpart.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        mpart.CreateNewNode(1, 0.0, 1.0, 0.0)
        mpart.CreateNewNode(2, 2.0, 0.0, 0.0)
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

    def testPressureCoefficient(self):
        model = self._CreateModel()
        mpart = model["main"]
        mpart.GetNode(1).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e6)
        mpart.GetNode(2).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1e4)

        params = self._GetBlankParameters()
        params["Parameters"].AddEmptyValue("freestream_pressure").SetDouble(1e5)
        params["Parameters"].AddEmptyValue("freestream_velocity").SetDouble(100)
        params["Parameters"].AddEmptyValue("freestream_density").SetDouble(1.2)
        params["Parameters"].AddEmptyValue("compute_force_coefficient").SetBool(False)

        process = compute_aerodynamic_coefficients_process.Factory(params, model)
        process.ExecuteInitialize()

        for node in mpart.Nodes:
            self.assertTrue(node.Has(KratosMultiphysics.PRESSURE_COEFFICIENT))
            self.assertAlmostEqual(node.GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 0.0)

        process.ExecuteFinalizeSolutionStep()

        self.assertAlmostEqual(mpart.GetNode(1).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), 150.0)
        self.assertAlmostEqual(mpart.GetNode(2).GetValue(KratosMultiphysics.PRESSURE_COEFFICIENT), -15.0)

    def testForceCoefficients(self):
        model = self._CreateModel()
        mpart = model["main"]
        mpart.GetNode(1).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1.5e5)
        mpart.GetNode(2).SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1.6e5)

        params = self._GetBlankParameters()
        p_inf = 1e5
        rho_inf = 2.4
        v_inf = 100
        params["Parameters"].AddEmptyValue("freestream_pressure").SetDouble(p_inf)
        params["Parameters"].AddEmptyValue("freestream_velocity").SetDouble(v_inf)
        params["Parameters"].AddEmptyValue("freestream_density").SetDouble(rho_inf)
        params["Parameters"].AddEmptyValue("compute_force_coefficient").SetBool(True)


        process = compute_aerodynamic_coefficients_process.Factory(params, model)
        process.ExecuteInitialize()
        process.ExecuteFinalizeSolutionStep()

        mean_pressure = 1.55e5
        normal_x = 2.0
        normal_y = 1.0

        lift = normal_y * (mean_pressure - p_inf)
        drag = normal_x * (mean_pressure - p_inf)
        dynamic_pressure = 0.5 * rho_inf * v_inf**2
        cord = 1

        cl = lift / (dynamic_pressure * cord)
        cd = drag / (dynamic_pressure * cord)

        self.assertAlmostEqual(mpart.GetProperties()[0].GetValue(KratosMultiphysics.DRAG_COEFFICIENT), cd)
        self.assertAlmostEqual(mpart.GetProperties()[0].GetValue(FluidDynamicsApplication.LIFT_COEFFICIENT), cl)


if __name__ == '__main__':
    KratosUnittest.main()
