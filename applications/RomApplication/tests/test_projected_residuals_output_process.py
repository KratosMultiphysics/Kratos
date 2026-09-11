import os
import shutil
import types
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities

if kratos_utilities.CheckIfApplicationsAvailable("FluidDynamicsApplication", "ConvectionDiffusionApplication"):
    import KratosMultiphysics.FluidDynamicsApplication
    import KratosMultiphysics.ConvectionDiffusionApplication

@KratosUnittest.skipIfApplicationsNotAvailable("FluidDynamicsApplication", "ConvectionDiffusionApplication")
class TestProjectedResidualsOutputProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.work_folder = "coupled_fluid_thermal_test_files/ROM/"
        self.parameters_filename = "../ProjectParameters.json"

        self.output_path_fluid = Path("rom_data/ResidualsFluid")
        self.output_path_thermal = Path("rom_data/ResidualsThermal")

    def testResidualsOutputProcess(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.parameters_filename, 'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            if not parameters.Has("output_processes"):
                parameters.AddEmptyValue("output_processes")
            if not parameters["output_processes"].Has("rom_output"):
                parameters["output_processes"].AddEmptyArray("rom_output")

            residuals_from_fluid_solver = """{
                "python_module": "projected_residuals_output_process",
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "ProjectedResidualsOutputProcess",
                "Parameters": {
                    "model_part_name": "FluidModelPart",
                    "output_control_type": "step",
                    "output_interval": 1,
                    "output_path": "rom_data/ResidualsFluid",
                    "sub_solver_name": "fluid_solver",
                    "range_of_entity_ids_to_fetch_residuals_projected": ["1", "50"]
                }
            }"""
            residuals_from_thermal_solver = """{
                "python_module": "projected_residuals_output_process",
                "kratos_module": "KratosMultiphysics.RomApplication",
                "process_name": "ProjectedResidualsOutputProcess",
                "Parameters": {
                    "model_part_name": "ThermalModelPart",
                    "output_control_type": "step",
                    "output_interval": 1,
                    "output_path": "rom_data/ResidualsThermal",
                    "sub_solver_name": "thermal_solver",
                    "range_of_entity_ids_to_fetch_residuals_projected": ["1", "End"]
                }
            }"""
            parameters["output_processes"]["rom_output"].Append(KratosMultiphysics.Parameters(residuals_from_fluid_solver))
            parameters["output_processes"]["rom_output"].Append(KratosMultiphysics.Parameters(residuals_from_thermal_solver))

            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)


            solution_to_impose = np.load("ExpectedOutputCoupledROM.npy")
            def SolveSolutionStep(cls):
                computing_model_part = cls._solver.GetComputingModelPart()
                step_index = computing_model_part.ProcessInfo[KratosMultiphysics.STEP] - 1

                if step_index < solution_to_impose.shape[1]:
                    current_step_solution = solution_to_impose[:, step_index]
                    variables_array = [
                        KratosMultiphysics.PRESSURE,
                        KratosMultiphysics.TEMPERATURE,
                        KratosMultiphysics.VELOCITY_X,
                        KratosMultiphysics.VELOCITY_Y
                    ]

                    counter = 0
                    for node in computing_model_part.Nodes:
                        for var in variables_array:
                            if node.SolutionStepsDataHas(var):
                                node.SetSolutionStepValue(var, current_step_solution[counter])
                                counter += 1

                return True

            self.simulation.SolveSolutionStep = types.MethodType(SolveSolutionStep, self.simulation)
            self.simulation.Run()
            steps_to_check = [1, 2, 3]

            for step in steps_to_check:
                # Fluid Residuals, checks for from element 1 to 50
                expected_file = Path(f"rom_data/ExpectedResidualsFluid/Residual_{step}.npy")
                obtained_file = self.output_path_fluid / f"Residual_{step}.npy"

                self.assertTrue(obtained_file.exists(), msg=f"File {obtained_file} was not generated.")
                expected = np.load(expected_file)
                obtained = np.load(obtained_file)
                for i in range(obtained.shape[0]):
                    for j in range(obtained.shape[1]):
                        self.assertAlmostEqual(obtained[i,j], expected[i,j])

                # Thermal Residuals, checks for from element 1 to End
                expected_file = Path(f"rom_data/ExpectedResidualsThermal/Residual_{step}.npy")
                obtained_file = self.output_path_thermal / f"Residual_{step}.npy"

                self.assertTrue(obtained_file.exists(), msg=f"File {obtained_file} was not generated.")
                expected = np.load(expected_file)
                obtained = np.load(obtained_file)
                for i in range(obtained.shape[0]):
                    for j in range(obtained.shape[1]):
                        self.assertAlmostEqual(obtained[i,j], expected[i,j])


    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Clean root time files
            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

            # Recursively delete the generated residuals folders so the workspace stays clean
            if self.output_path_fluid.exists():
                shutil.rmtree(self.output_path_fluid)
            if self.output_path_thermal.exists():
                shutil.rmtree(self.output_path_thermal)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()