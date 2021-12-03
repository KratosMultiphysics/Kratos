import os
import json
import types

try:
    import numpy as np
    numpy_available = True
except:
    numpy_available = False

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
import KratosMultiphysics.RomApplication as KratosROM
from KratosMultiphysics.RomApplication.rom_analysis import CreateRomAnalysisInstance
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities
if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication

@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
class TestStructuralRom(KratosUnittest.TestCase):

    def setUp(self):
        self.relative_tolerance = 1.0e-12

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def testStructuralStaticRom2D(self):
        self.work_folder = "structural_static_test_files"
        parameters_filename = "ProjectParametersROM.json"
        expected_output_filename = "ExpectedOutputROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Run test case
            self.simulation.Run()

            # Check results
            expected_output = np.load(expected_output_filename)
            obtained_output = rom_testing_utilities.GetVectorNodalResults(self.simulation._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)
            aux_nodal_area = rom_testing_utilities.GetNodalAreaVector(self.simulation._GetSolver().GetComputingModelPart())
            nodal_area = [aux_nodal_area[int((i-1)/2) if i%2 != 0 else int(i/2)] for i in range(2*len(aux_nodal_area))]

            numerator = 0.0
            denominator = 0.0
            for i in range(len(obtained_output)):
                if expected_output[i] > 1.0e-12:
                    numerator += nodal_area[i]*((1 - (obtained_output[i]/expected_output[i]))**2)
                    denominator += nodal_area[i]
            l2 = np.sqrt(numerator/denominator)*100
            self.assertLess(l2, self.relative_tolerance)

    @KratosUnittest.skipUnless(numpy_available, "numpy is required for RomApplication")
    def testStructuralStaticHRom2D(self):
        self.work_folder = "structural_static_test_files"
        parameters_filename = "ProjectParametersHROM.json"
        expected_output_filename = "ExpectedOutputHROM.npy"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up and run simulation
            with open(parameters_filename,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            self.simulation = rom_testing_utilities.SetUpSimulationInstance(model, parameters)

            # Patch the RomAnalysis class to set the HROM weights
            def ModifyAfterSolverInitialize(cls):
                super(type(self.simulation), cls).ModifyAfterSolverInitialize()
                computing_model_part = cls._GetSolver().GetComputingModelPart()
                print(computing_model_part)

                with open('ElementsAndWeights.json') as f:
                    HR_data = json.load(f)
                    print(HR_data)
                    for key in HR_data["Elements"].keys():
                        computing_model_part.GetElement(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, HR_data["Elements"][key])
                    for key in HR_data["Conditions"].keys():
                        computing_model_part.GetCondition(int(key)+1).SetValue(KratosROM.HROM_WEIGHT, HR_data["Conditions"][key])
                print("Done")

            self.simulation.ModifyAfterSolverInitialize  = types.MethodType(ModifyAfterSolverInitialize, self.simulation)

            print("here!")

            # Run test case
            self.simulation.Run()

            print("here2")

            # Check results
            expected_output = np.load(expected_output_filename)
            obtained_output = rom_testing_utilities.GetVectorNodalResults(self.simulation._GetSolver().GetComputingModelPart(), KratosMultiphysics.DISPLACEMENT)
            aux_nodal_area = rom_testing_utilities.GetNodalAreaVector(self.simulation._GetSolver().GetComputingModelPart())
            nodal_area = [aux_nodal_area[int((i-1)/2) if i%2 != 0 else int(i/2)] for i in range(2*len(aux_nodal_area))]

            numerator = 0.0
            denominator = 0.0
            for i in range(len(obtained_output)):
                if expected_output[i] > 1.0e-12:
                    numerator += nodal_area[i]*((1 - (obtained_output[i]/expected_output[i]))**2)
                    denominator += nodal_area[i]
            l2 = np.sqrt(numerator/denominator)*100
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
