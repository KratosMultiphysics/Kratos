import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSteepestDescentRestart(kratos_unittest.TestCase):
    restart_files_path = "Optimization_Restart_test"
    checkpoint_only_restart_files_path = "Optimization_Restart_checkpoint_only_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

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
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_control_field = list(reference_algorithm.GetCurrentControlField().data)

            # checkpoint run: solve iterations 0..5 only, writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 5)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"           : true,
                "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_save_frequency" : 1
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 5 and run one more iteration to 6.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 6)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"       : true,
                "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_load_step"  : 5
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_control_field = list(resume_algorithm.GetCurrentControlField().data)

            self.assertAlmostEqual(reference_obj_value, resume_obj_value, places=9)
            # YOUNG_MODULUS-scale field (~1e10): compare with a delta relative to its magnitude
            max_abs_reference_value = max(abs(value) for value in reference_control_field)
            self.assertVectorAlmostEqual(reference_control_field, resume_control_field, places=None, delta=max_abs_reference_value * 1e-9)

    def test_resume_with_checkpoint_only_component(self):
        """Data of a component absent from the resume run is restored verbatim: nodes are recreated
        from node_metadata, tensors relink to the restored Model, and the run proceeds."""
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

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)
            DeleteDirectoryIfExisting(cls.checkpoint_only_restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
