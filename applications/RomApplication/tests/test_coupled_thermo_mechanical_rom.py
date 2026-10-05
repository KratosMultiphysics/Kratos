import os
import types
import numpy as np
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities


if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "ConvectionDiffusionApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
    import KratosMultiphysics.ConvectionDiffusionApplication

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "ConvectionDiffusionApplication")
class TestCoupledThermoMechanicalRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testCoupledThermoMechanicalRom2D(self):
        self.work_folder = "coupled_thermo_mechanical_test_files/ROM/"
        self._RunAndCheck("ExpectedOutputCoupledROM.npy")

    def testCoupledThermoMechanicalHRom2D(self):
        self.work_folder = "coupled_thermo_mechanical_test_files/HROM/"
        self._RunAndCheck("ExpectedOutputCoupledHROM.npy")

    def testCoupledThermoMechanicalAnnEnhancedRom2D(self):
        # ANN-enhanced ROM run without the RomManager: each coupled solver has its own basis and network
        self.work_folder = "coupled_thermo_mechanical_test_files/ROM_ANN/"
        self._RunAndCheck("ExpectedOutputCoupledROM_ANN.npy")

    def _GetProjectParameters(self):
        # The problem (parameters, mdpa and materials) is the thermo-mechanical test of the ConvectionDiffusionApplication.
        # Its paths are relative to the tests folder of that application, so they are made absolute
        convection_diffusion_tests_folder = Path(__file__).resolve().parents[2] / "ConvectionDiffusionApplication" / "tests"
        with open(convection_diffusion_tests_folder / "thermo_mechanical_tests/thermo_mechanical/coupled_problem_test_parameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        parameters.AddString("analysis_stage", "KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis")
        for sub_solver_settings_name in ["structural_solver_settings", "thermal_solver_settings"]:
            sub_solver_settings = parameters["solver_settings"][sub_solver_settings_name]
            for file_name in [sub_solver_settings["model_import_settings"]["input_filename"], sub_solver_settings["material_import_settings"]["materials_filename"]]:
                file_name.SetString(str(convection_diffusion_tests_folder / file_name.GetString()))
        parameters["processes"].RemoveValue("json_check_process") # It checks against the FOM results
        return parameters

    def _RunSimulation(self):
        # Set up simulation
        model = KratosMultiphysics.Model()
        self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, self._GetProjectParameters())

        # Patch the RomAnalysis class to save the selected time steps results
        self.simulation.selected_time_step_solution_container = []

        def FinalizeSolutionStep(cls):
            super(type(self.simulation), cls).FinalizeSolutionStep()

            variables_array = [KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.TEMPERATURE]
            array_of_results = rom_testing_utilities.GetNodalResults(cls._solver.GetComputingModelPart(), variables_array)
            cls.selected_time_step_solution_container.append(array_of_results)

        self.simulation.FinalizeSolutionStep  = types.MethodType(FinalizeSolutionStep, self.simulation)

        # Run test case
        self.simulation.Run()

        return np.array(self.simulation.selected_time_step_solution_container).T

    def _RunAndCheck(self, expected_output_filename):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            obtained_snapshot_matrix = self._RunSimulation()

            # Check results
            expected_output = np.load(expected_output_filename)
            self.assertEqual(obtained_snapshot_matrix.shape, expected_output.shape)
            for i in range(expected_output.shape[1]):
                up = sum((expected_output[:,i] - obtained_snapshot_matrix[:,i])**2)
                down = sum((expected_output[:,i])**2)
                l2 = np.sqrt(up/down)
                self.assertLess(l2, self.relative_tolerance)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Cleaning
            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
