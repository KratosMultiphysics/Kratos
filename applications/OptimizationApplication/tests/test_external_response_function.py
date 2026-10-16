import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting, DeleteDirectoryIfExisting

import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestExternalResponseFunction(kratos_unittest.TestCase):
    restart_files_path = "external_response_function_restart"

    @staticmethod
    def _ReadParameters(max_iter: int, summary_file_name: str = "", restart_settings: str = "") -> Kratos.Parameters:
        with open("auxiliary_files/external_response_function_optimization_parameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)
        if summary_file_name:
            parameters["processes"]["optimization_data_processes"].AddValue("output_processes", Kratos.Parameters("""[{
                "type"    : "optimization_problem_ascii_output_process",
                "module"  : "KratosMultiphysics.OptimizationApplication.processes",
                "settings": {
                    "output_file_name"    : \"""" + summary_file_name + """\",
                    "write_kratos_version": false,
                    "write_time_stamp"    : false
                }
            }]"""))
        if restart_settings:
            parameters.AddValue("restart_settings", Kratos.Parameters(restart_settings))
        return parameters

    @staticmethod
    def _GetDesign(model: Kratos.Model) -> 'list[float]':
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model["CustomMP"].Nodes, KratosOA.CUSTOM_DESIGN_VARIABLE)
        ta.CollectData()
        return list(ta.data)

    @staticmethod
    def _ReadSummaryRows(output_file_name: str) -> 'list[str]':
        with open(output_file_name, "r") as file_input:
            return [line for line in file_input if not line.startswith("#")]

    def test_steepest_descent_analysis(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            model = Kratos.Model()
            OptimizationAnalysis(model, self._ReadParameters(10)).Run()
            self.assertVectorAlmostEqual(self._GetDesign(model), [2,3,4,5])

    def test_steepest_descent_restart(self):
        """Checkpoints with the per-point buffered values of ExternalResponseFunctionControl, and a
        resumed run reproduces the uninterrupted summary rows and design."""
        with kratos_unittest.WorkFolderScope(".", __file__):
            reference_model = Kratos.Model()
            OptimizationAnalysis(reference_model, self._ReadParameters(10, "external_reference.csv")).Run()

            restart_settings = """{
                "save_restart"     : true,
                "load_restart"     : true,
                "restart_file_name": \"""" + self.restart_files_path + """/restart_<step>.pkl"
            }"""
            # no checkpoint yet: iterations 0..4 from scratch, then resume from step 4 up to 10
            OptimizationAnalysis(Kratos.Model(), self._ReadParameters(4, "external_restart.csv", restart_settings)).Run()
            resume_model = Kratos.Model()
            OptimizationAnalysis(resume_model, self._ReadParameters(10, "external_restart.csv", restart_settings)).Run()

            self.assertVectorAlmostEqual(self._GetDesign(resume_model), self._GetDesign(reference_model))
            self.assertEqual(self._ReadSummaryRows("external_restart.csv"), self._ReadSummaryRows("external_reference.csv"))

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("external_reference.csv")
            DeleteFileIfExisting("external_restart.csv")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
