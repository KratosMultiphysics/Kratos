import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_analysis
from KratosMultiphysics.json_output_process import JsonOutputProcess
from KratosMultiphysics.from_json_check_result_process import FromJsonCheckResultProcess

import os


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)


class DiscontinuityPropagationConvectionDiffusionAnalysis(convection_diffusion_analysis.ConvectionDiffusionAnalysis):
    """Derived convection-diffusion analysis stage."""

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        ## Set initial and boundary conditions according to the test configuration
        model_part = self.model["ThermalModelPart"]
        t = model_part.ProcessInfo[KratosMultiphysics.TIME]

        for node in self.model["ThermalModelPart.Top_Wall"].Nodes:
            if node.X <= 4.0 * t:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1.0)
                node.Fix(KratosMultiphysics.TEMPERATURE)

        for node in self.model["ThermalModelPart.Left_Wall"].Nodes:
            if node.Y >= 1.0 - 2.5 * t:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1.0)
                node.Fix(KratosMultiphysics.TEMPERATURE)


class DiscontinuityPropagationUnitSquareTest(UnitTest.TestCase):
    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            ## Set up the custom source term input
            self.model = KratosMultiphysics.Model()
            with open(self.parameters,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            ## Run test
            self.conv_diffusion_test = DiscontinuityPropagationConvectionDiffusionAnalysis(self.model, parameters)
            self.conv_diffusion_test.Run()

    def checkResults(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            if self.print_reference_values:
                json_output_settings = KratosMultiphysics.Parameters(r'''{
                        "output_variables": ["TEMPERATURE"],
                        "output_file_name": "",
                        "model_part_name": "ThermalModelPart",
                        "time_frequency": 0.0
                }''')
                json_output_settings["output_file_name"].SetString(GetFilePath(self.work_folder + "/" + self.reference_file + "_results.json"))
                json_output_process = JsonOutputProcess(self.model, json_output_settings)
                json_output_process.ExecuteInitialize()
                json_output_process.ExecuteBeforeSolutionLoop()
                json_output_process.ExecuteFinalizeSolutionStep()
            else:
                json_check_parameters = KratosMultiphysics.Parameters(r'''{
                    "check_variables"      : ["TEMPERATURE"],
                    "input_file_name"      : "",
                    "model_part_name"      : "ThermalModelPart",
                    "time_frequency"       : 0.0
                }''')
                json_check_parameters["input_file_name"].SetString(GetFilePath(self.work_folder + "/" + self.reference_file + "_results.json"))
                json_check_process = FromJsonCheckResultProcess(self.model, json_check_parameters)
                json_check_process.ExecuteInitialize()
                json_check_process.ExecuteBeforeSolutionLoop()
                json_check_process.ExecuteFinalizeSolutionStep()

    def testDiscontinuityPropagation(self):
        self.print_reference_values = False
        self.work_folder = "discontinuity_propagation_unit_square_test"
        self.parameters = "conv_diffusion_test_parameters.json"
        # self.reference_file = "discontinuity_propagation_unit_square_wo_crosswind_test"
        self.reference_file = "discontinuity_propagation_unit_square_w_crosswind_test"
        self.runTest()
        self.checkResults()

if __name__ == '__main__':
    test = DiscontinuityPropagationUnitSquareTest()
    test.testDiscontinuityPropagation()

