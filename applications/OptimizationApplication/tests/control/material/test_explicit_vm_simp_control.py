import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

class TestExplicitVertexMorphingSimpControl(kratos_unittest.TestCase):
    def test_SimpControl(self):
        with kratos_unittest.WorkFolderScope("../../mdpas", __file__):
            with open("../control/material/optimization_parameters.json",'r') as parameter_file:
                parameters = KM.Parameters(parameter_file.read())
            model = KM.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

            CompareTwoFilesCheckProcess(KM.Parameters("""
            {
                "reference_file_name"   : "../control/material/summary_orig.csv",
                "output_file_name"      : "summary.csv",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")).Execute()

if __name__ == "__main__":
    KM.Tester.SetVerbosity(KM.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()