import os
import types
import numpy as np

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.RomApplication.rom_manager import RomManager
from pathlib import Path

def CustomizeSimulation(cls, global_model, parameters):

    class CustomSimulation(cls):

        def GetFinalData(self):
            super().GetFinalData()
            #taken from the GeoMechanics app tests folder
            return {f"displacement": rom_testing_utilities.GetVectorNodalResults(self._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT), f"nodal_area": rom_testing_utilities.GetNodalAreaVector(self._GetSolver().GetComputingModelPart()) }

    return CustomSimulation(global_model, parameters)


general_rom_manager_parameters = KratosMultiphysics.Parameters("""{
            "projection_strategy": "galerkin"
        }""")



@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestStructuralRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    def testStructuralStaticRom2D(self):
        self.work_folder = "structural_static_test_files/ROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name=parameters_filename, CustomizeSimulation = CustomizeSimulation, general_rom_manager_parameters=general_rom_manager_parameters)
            rom_manager.RunROM()
            # Check results
            expected_output = np.load(expected_output_filename)
            obtained_output = rom_manager.QoI_Run_ROM[0]["displacement"]
            aux_nodal_area = rom_manager.QoI_Run_ROM[0]["nodal_area"]
            nodal_area = [aux_nodal_area[int((i-1)/2) if i%2 != 0 else int(i/2)] for i in range(2*len(aux_nodal_area))]

            numerator = 0.0
            denominator = 0.0
            for i in range(len(obtained_output)):
                if expected_output[i] > 1.0e-12:
                    numerator += nodal_area[i]*((1 - (obtained_output[i]/expected_output[i]))**2)
                    denominator += nodal_area[i]
            l2 = np.sqrt(numerator/denominator)*100
            self.assertLess(l2, self.relative_tolerance)

    def testStructuralStaticHRom2D(self):
        self.work_folder = "structural_static_test_files/HROM/"
        parameters_filename = "ProjectParametersHROM.json"
        expected_output_filename = "ExpectedOutputHROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            rom_manager = RomManager(project_parameters_name=parameters_filename, CustomizeSimulation = CustomizeSimulation)
            rom_manager.RunHROM(use_full_model_part = True)

            # Check results
            expected_output = np.load(expected_output_filename)
            obtained_output = rom_manager.QoI_Run_HROM[0]["displacement"]
            aux_nodal_area = rom_manager.QoI_Run_HROM[0]["nodal_area"]
            nodal_area = [aux_nodal_area[int((i-1)/2) if i%2 != 0 else int(i/2)] for i in range(2*len(aux_nodal_area))]

            numerator = 0.0
            denominator = 0.0
            for i in range(len(obtained_output)):
                if expected_output[i] > 1.0e-12:
                    numerator += nodal_area[i]*((1 - (obtained_output[i]/expected_output[i]))**2)
                    denominator += nodal_area[i]
            l2 = np.sqrt(numerator/denominator)*100
            self.assertLess(l2, self.relative_tolerance)

    def testStructuralDynamicRom2D(self):
        self.work_folder = "structural_dynamic_test_files/ROM/"
        parameters_filename = "../ProjectParameters.json"
        expected_output_filename = "ExpectedOutputROM.npy"

        time_snapshots = [2,4,6,8,10]

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up and run simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Patch the RomAnalysis class to save the selected time steps results
            def Initialize(cls):
                super(type(self.simulation), cls).Initialize()
                cls.selected_time_step_solution_container = []

            def FinalizeSolutionStep(cls):
                super(type(self.simulation), cls).FinalizeSolutionStep()

                time = cls._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
                if np.any(np.isclose(time, time_snapshots)):
                    array_of_displacements = rom_testing_utilities.GetVectorNodalResults(self.simulation._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)
                    cls.selected_time_step_solution_container.append(array_of_displacements)

            self.simulation.Initialize  = types.MethodType(Initialize, self.simulation)
            self.simulation.FinalizeSolutionStep  = types.MethodType(FinalizeSolutionStep, self.simulation)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            n_values = len(self.simulation.selected_time_step_solution_container[0])
            n_snapshots = len(self.simulation.selected_time_step_solution_container)
            obtained_snapshot_matrix = np.zeros((n_values, n_snapshots))
            for i in range(n_snapshots):
                snapshot_i= np.array(self.simulation.selected_time_step_solution_container[i])
                obtained_snapshot_matrix[:,i] = snapshot_i.transpose()

            tolerance = 1.0e-10
            for i in range (n_snapshots):
                up = sum((expected_output[:,i] - obtained_snapshot_matrix[:,i])**2)
                down = sum((expected_output[:,i])**2)
                l2 = np.sqrt(up/down)
                self.assertLess(l2, tolerance)

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Cleaning
            for file_name in os.listdir():
                if file_name.endswith(".time"):
                    kratos_utilities.DeleteFileIfExisting(file_name)
        workfolders = ["structural_static_test_files/HROM/","structural_static_test_files/ROM/"]
        for workfolder in workfolders:
            with KratosUnittest.WorkFolderScope(workfolder, __file__):
                kratos_utilities.DeleteDirectoryIfExisting(Path('./rom_data/rom_database'))

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
