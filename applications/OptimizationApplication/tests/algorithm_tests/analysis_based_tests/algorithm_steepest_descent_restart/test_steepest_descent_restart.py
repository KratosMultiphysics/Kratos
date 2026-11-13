import json
import subprocess
import sys
from unittest import mock

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSteepestDescentRestart(kratos_unittest.TestCase):
    restart_files_path = "Optimization_Restart_test"
    checkpoint_only_restart_files_path = "Optimization_Restart_checkpoint_only_test"
    fresh_interpreter_restart_files_path = "Optimization_Restart_fresh_interpreter_test"
    fresh_interpreter_result = "fresh_interpreter_result.json"
    reference_summary = "summary_reference.csv"
    restart_summary = "summary_restart.csv"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

    def _AddSummaryOutput(self, parameters: Kratos.Parameters, output_file_name: str) -> None:
        parameters["processes"]["optimization_data_processes"].AddValue("output_processes", Kratos.Parameters("""[{
            "type"    : "optimization_problem_ascii_output_process",
            "module"  : "KratosMultiphysics.OptimizationApplication.processes",
            "settings": {
                "output_file_name"    : \"""" + output_file_name + """\",
                "write_kratos_version": false,
                "write_time_stamp"    : false
            }
        }]"""))

    @staticmethod
    def _ReadSummaryRows(output_file_name: str) -> 'list[list[float]]':
        with open(output_file_name, "r") as file_input:
            return [[float(v) for v in line.split(",")] for line in file_input if not line.startswith("#")]

    def _RunToConvergence(self, parameters: Kratos.Parameters) -> OptimizationAnalysis:
        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()
        return analysis

    def test_steepest_descent_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..6 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 6)
            self._AddSummaryOutput(reference_parameters, self.reference_summary)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_control_field = list(reference_algorithm.GetCurrentControlField().data)

            # checkpoint run: solve iterations 0..5 only, writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 5)
            self._AddSummaryOutput(checkpoint_parameters, self.restart_summary)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"           : true,
                "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_save_frequency" : 1
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 4 and run iterations 5..6; the summary
            # keeps rows 0..4 of the checkpoint run, drops its row 5 and appends the resumed rows.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 6)
            self._AddSummaryOutput(resume_parameters, self.restart_summary)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"       : true,
                "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_load_step"  : 4
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_control_field = list(resume_algorithm.GetCurrentControlField().data)

            self.assertAlmostEqual(reference_obj_value, resume_obj_value, places=9)
            # YOUNG_MODULUS-scale field (~1e10): compare with a delta relative to its magnitude
            max_abs_reference_value = max(abs(value) for value in reference_control_field)
            self.assertVectorAlmostEqual(reference_control_field, resume_control_field, places=None, delta=max_abs_reference_value * 1e-9)

            reference_rows = self._ReadSummaryRows(self.reference_summary)
            restart_rows = self._ReadSummaryRows(self.restart_summary)
            self.assertEqual([row[0] for row in restart_rows], list(range(7)))
            for reference_row, restart_row in zip(reference_rows, restart_rows):
                self.assertVectorAlmostEqual(reference_row, restart_row, places=None, delta=1e-6 * max(abs(v) for v in reference_row))

    def test_resume_in_fresh_interpreter(self):
        """Resumes in a new Python process that imported no application besides OptimizationApplication:
        the checkpoint must import the applications whose elements and laws the Model holds."""
        with kratos_unittest.WorkFolderScope(".", __file__):
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 6)
            reference_algorithm = self._RunToConvergence(reference_parameters).GetAlgorithm()

            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 4)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"     : true,
                "restart_file_name": \"""" + self.fresh_interpreter_restart_files_path + """/restart_<step>.pkl"
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            resume_script = f"""
import json, sys
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
assert "KratosMultiphysics.StructuralMechanicsApplication" not in sys.modules
with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(6)
parameters.AddValue("restart_settings", Kratos.Parameters('{{"load_restart": true, "restart_file_name": "{self.fresh_interpreter_restart_files_path}/restart_<step>.pkl"}}'))
analysis = OptimizationAnalysis(Kratos.Model(), parameters)
analysis.Run()
with open("{self.fresh_interpreter_result}", "w") as file_output:
    json.dump({{"obj": analysis.GetAlgorithm().GetOptimizedObjectiveValue(), "control_field": list(analysis.GetAlgorithm().GetCurrentControlField().data)}}, file_output)
"""
            process = subprocess.run([sys.executable, "-c", resume_script], capture_output=True, text=True)
            self.assertEqual(process.returncode, 0, msg=process.stderr[-3000:])
            with open(self.fresh_interpreter_result, "r") as file_input:
                result = json.load(file_input)

            self.assertAlmostEqual(reference_algorithm.GetOptimizedObjectiveValue(), result["obj"], places=9)
            reference_control_field = list(reference_algorithm.GetCurrentControlField().data)
            max_abs_reference_value = max(abs(value) for value in reference_control_field)
            self.assertVectorAlmostEqual(reference_control_field, result["control_field"], places=None, delta=max_abs_reference_value * 1e-9)

    def test_resume_with_checkpoint_only_component(self):
        """Data of a component absent from the resume run is restored verbatim: nodes (also empty ones)
        are recreated from node_metadata, tensors relink to the restored Model, and the run proceeds."""
        with kratos_unittest.WorkFolderScope(".", __file__):
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 6)
            reference_algorithm = self._RunToConvergence(reference_parameters).GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()

            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 5)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"     : true,
                "restart_file_name": \"""" + self.checkpoint_only_restart_files_path + """/restart_<step>.pkl"
            }"""))
            checkpoint_model = Kratos.Model()
            checkpoint_analysis = OptimizationAnalysis(checkpoint_model, checkpoint_parameters)
            checkpoint_analysis.Initialize()
            checkpoint_problem = checkpoint_analysis.optimization_problem
            orphan_view = ComponentDataView("checkpoint_only_control", checkpoint_problem)
            orphan_field = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(checkpoint_model["Structure"].Elements, Kratos.YOUNG_MODULUS)
            orphan_field.CollectData()
            orphan_view.GetUnBufferedData().SetValue("field", orphan_field)
            orphan_view.GetUnBufferedData().SetValue("nested/scalar", 3.5)
            checkpoint_problem.GetProblemDataContainer()[f"{orphan_view.GetDataPath()}/buffered"] = BufferedDict(3, False)
            orphan_view.GetBufferedData().SetValue("history", 7)
            # a component whose 3-slot history buffer holds no values yet
            ComponentDataView("checkpoint_only_empty", checkpoint_problem).SetDataBuffer(3)
            checkpoint_analysis.Check()
            checkpoint_analysis.GetAlgorithm().Solve()
            checkpoint_analysis.Finalize()

            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 6)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"     : true,
                "restart_file_name": \"""" + self.checkpoint_only_restart_files_path + """/restart_<step>.pkl",
                "restart_load_step": 5
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            self.assertAlmostEqual(reference_obj_value, resume_analysis.GetAlgorithm().GetOptimizedObjectiveValue(), places=9)

            resume_problem = resume_analysis.optimization_problem
            orphan_view = ComponentDataView("checkpoint_only_control", resume_problem)
            restored_field = orphan_view.GetUnBufferedData()["field"]
            self.assertIsInstance(restored_field, KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor)
            self.assertEqual(restored_field.GetContainer(), resume_analysis.model["Structure"].Elements)
            self.assertVectorAlmostEqual(restored_field.data, orphan_field.data)
            self.assertEqual(orphan_view.GetUnBufferedData()["nested/scalar"], 3.5)
            restored_buffered = orphan_view.GetBufferedData()
            self.assertEqual(restored_buffered.GetBufferSize(), 3)
            self.assertFalse(restored_buffered.GetClearBufferWhenAdvancing())
            self.assertTrue(any(restored_buffered.HasValue("history", i) for i in range(3)))
            empty_view = ComponentDataView("checkpoint_only_empty", resume_problem)
            self.assertTrue(empty_view.HasDataBuffer())
            self.assertEqual(empty_view.GetBufferedData().GetBufferSize(), 3)
            self.assertEqual(empty_view.GetBufferedData().GetValueItems(), {})

    def test_restart_with_unsupported_algorithm_raises(self):
        """Restart settings on an algorithm without SUPPORTS_RESTART raise in the constructor."""
        with kratos_unittest.WorkFolderScope(".", __file__):
            parameters = self._ReadParameters()
            parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"     : true,
                "restart_file_name": \"""" + self.restart_files_path + """/restart_<step>.pkl"
            }"""))
            with mock.patch.object(AlgorithmSteepestDescent, "SUPPORTS_RESTART", False):
                with self.assertRaisesRegex(RuntimeError, "does not call SaveRestartFile"):
                    OptimizationAnalysis(Kratos.Model(), parameters)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting(cls.reference_summary)
            DeleteFileIfExisting(cls.restart_summary)
            DeleteDirectoryIfExisting(cls.restart_files_path)
            DeleteDirectoryIfExisting(cls.checkpoint_only_restart_files_path)
            DeleteDirectoryIfExisting(cls.fresh_interpreter_restart_files_path)
            DeleteFileIfExisting(cls.fresh_interpreter_result)

if __name__ == "__main__":
    kratos_unittest.main()
