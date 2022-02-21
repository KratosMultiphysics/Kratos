import KratosMultiphysics
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface
from KratosMultiphysics.FreeSurfaceApplication.free_surface_analysis import FreeSurfaceAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
import pathlib
import shutil


class KratosFreeSurfaceGeneralTests(KratosUnittest.TestCase):

    def ExecuteReservoirTests(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.setUp()
            self.setUpProblem()
            self.runTest()
            self.tearDown()
            self.checkResults(self.reference_file)

    def test_Reservoir2D(self):
        self.__RunFreeSurfaceAnalysis("Reservoir2D", "reference_reservoir_2D.csv")

    def test_Reservoir3D(self):
        self.__RunFreeSurfaceAnalysis("Reservoir3D", "reference_reservoir_3D.csv")

    def __RunFreeSurfaceAnalysis(self, test_name: str, reference_name: str):
        with KratosUnittest.WorkFolderScope(test_name, __file__):
            parameters_path = pathlib.Path("ProjectParameters.json")
            with open(parameters_path, 'r') as file:
                parameters = KratosMultiphysics.Parameters(file.read())

            model = KratosMultiphysics.Model()

            analysis = FreeSurfaceAnalysis(model, parameters)
            analysis.Run()

            self.main_model_part = analysis._GetSolver().GetComputingModelPart()
            self.checkResults(reference_name)

    def setUp(self):
        self.check_tolerance = 1e-5
        self.print_reference_values = False

    def checkResults(self, reference_file_name: str):
        """Write results to a temporary file and compare its contents to the reference file."""
        output_file_name = "free_surface_general_tests_tmp.csv"
        with open(output_file_name, 'w') as output_file:
            output_file.write("#ID, DISTANCE\n")
            for node in self.main_model_part.Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,0)
                output_file.write("{0}, {1}\n".format(node.Id, dist))

        # Overwrite reference values if requested
        if self.print_reference_values:
            shutil.copy(output_file_name, reference_file_name)

        # Compare to reference file, and delete temporary output
        comparison_parameters = KratosMultiphysics.Parameters("""{
            "comparison_type"       : "csv_file",
            "remove_output_file"    : true
        }""")
        comparison_parameters.AddString("reference_file_name", reference_file_name)
        comparison_parameters.AddString("output_file_name", output_file_name)
        comparison_parameters.AddDouble("tolerance", self.check_tolerance)
        CompareTwoFilesCheckProcess(comparison_parameters).Execute()




if __name__ == "__main__":
    KratosUnittest.main()