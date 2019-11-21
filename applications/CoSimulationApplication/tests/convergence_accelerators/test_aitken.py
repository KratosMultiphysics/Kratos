import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import numpy as np


class TestConvergenceAcceleratorAITKEN(KratosUnittest.TestCase):
    def test_convergence_accelerator_aitken(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        m = 10
        dz = 2.0
        r = 0.1
        x = 10.0
        xt0 = 10.5
        xt1 = 10.2
        xt2 = 10.1
        xt3 = 10.7
        xt4 = 9.9

        interface_settings = cs_data_structure.Parameters('{"wall": "AREA"}')

        # Create interface
        variable = vars(KM)["AREA"]
        model = cs_data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, x)
        interface = CoSimulationInterface(model, interface_settings)

        parameter_file_name = "convergence_accelerators/test_aitken.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())

        omega_max = settings["settings"]["omega_max"].GetDouble()

        convergence_accelerator = cs_tools.CreateInstance(settings)
        convergence_accelerator.Initialize()
        convergence_accelerator.InitializeSolutionStep()

        interface_r = interface.deepcopy()
        interface_r.SetNumpyArray(r * np.ones(m))
        interface_x = interface.deepcopy()
        interface_x.SetNumpyArray(x * np.ones(m))
        interface_xt0 = interface.deepcopy()
        interface_xt0.SetNumpyArray(xt0 * np.ones(m))
        interface_xt1 = interface.deepcopy()
        interface_xt1.SetNumpyArray(xt1 * np.ones(m))
        interface_xt2 = interface.deepcopy()
        interface_xt2.SetNumpyArray(xt2 * np.ones(m))
        interface_xt3 = interface.deepcopy()
        interface_xt3.SetNumpyArray(xt3 * np.ones(m))
        interface_xt4 = interface.deepcopy()
        interface_xt4.SetNumpyArray(xt4 * np.ones(m))

        # Test value of self.added
        is_ready = convergence_accelerator.IsReady()
        self.assertFalse(is_ready)

        # Test Update()
        convergence_accelerator.Update(interface_x, interface_xt0)
        is_ready = convergence_accelerator.IsReady()
        self.assertTrue(is_ready)
        omega = convergence_accelerator.omega
        self.assertEqual(omega, omega_max)
        convergence_accelerator.Update(interface_x, interface_xt1)
        is_ready = convergence_accelerator.IsReady()
        self.assertTrue(is_ready)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, 5 / 3 * omega_max, 10)
        convergence_accelerator.Update(interface_x, interface_xt2)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, 10 / 3 * omega_max, 10)

        # Test Predict()
        interface_dx = convergence_accelerator.Predict(interface_r)
        omega = interface_dx.GetNumpyArray()[0] / r
        self.assertEqual(omega, convergence_accelerator.omega)
        convergence_accelerator.Predict(interface_r)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, 10 / 3 * omega_max, 10)

        # New solution step
        convergence_accelerator.FinalizeSolutionStep()
        convergence_accelerator.InitializeSolutionStep()

        # Test value of self.added
        is_ready = convergence_accelerator.IsReady()
        self.assertFalse(is_ready)

        # Test Update()
        convergence_accelerator.Update(interface_x, interface_xt0)
        is_ready = convergence_accelerator.IsReady()
        self.assertTrue(is_ready)
        omega = convergence_accelerator.omega
        self.assertEqual(omega, omega_max)
        convergence_accelerator.Update(interface_x, interface_xt3)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, -5 / 2 * omega_max, 10)

        # New solution step
        convergence_accelerator.FinalizeSolutionStep()
        convergence_accelerator.InitializeSolutionStep()

        # Test Update()
        convergence_accelerator.Update(interface_x, interface_xt0)
        omega = convergence_accelerator.omega
        self.assertEqual(omega, -omega_max)
        convergence_accelerator.Update(interface_x, interface_xt4)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, -5 / 6 * omega_max, 10)

        # New solution step
        convergence_accelerator.FinalizeSolutionStep()
        convergence_accelerator.InitializeSolutionStep()

        # Test Update()
        convergence_accelerator.Update(interface_x, interface_xt0)
        omega = convergence_accelerator.omega
        self.assertAlmostEqual(omega, -5 / 6 * omega_max, 10)


if __name__ == '__main__':
    KratosUnittest.main()
