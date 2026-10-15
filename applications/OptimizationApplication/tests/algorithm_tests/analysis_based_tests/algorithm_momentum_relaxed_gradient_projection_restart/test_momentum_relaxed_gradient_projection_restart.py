import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestMomentumRelaxedGradientProjectionRestart(kratos_unittest.TestCase):
    """Restart test for momentum RGP (momentum lives in the buffered data). The explicit filter is
    deliberate: the implicit filter's HelmholtzElement currently crashes on Serializer load."""
    restart_files_path = "Optimization_Restart_test"

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

    def _GetPhysicalField(self, analysis: OptimizationAnalysis) -> 'list[float]':
        """Thickness read off the mesh."""
        control, = analysis.GetAlgorithm().master_control.GetListOfControls()
        return list(control.GetPhysicalField().data)

    def test_momentum_relaxed_gradient_projection_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            # reference run: solve iterations 0..3 in one uninterrupted process.
            reference_parameters = self._ReadParameters()
            self._SetMaxIter(reference_parameters, 4)
            reference_analysis = self._RunToConvergence(reference_parameters)
            reference_algorithm = reference_analysis.GetAlgorithm()
            reference_obj_value = reference_algorithm.GetOptimizedObjectiveValue()
            reference_physical_field = self._GetPhysicalField(reference_analysis)

            # checkpoint run: solve iterations 0..2 only (momentum is already non-trivial by
            # then), writing a restart checkpoint every step.
            checkpoint_parameters = self._ReadParameters()
            self._SetMaxIter(checkpoint_parameters, 3)
            checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "save_restart"           : true,
                "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_save_frequency" : 1
            }"""))
            self._RunToConvergence(checkpoint_parameters)

            # resume run: restore the checkpoint from iteration 3 and run one more iteration to 4.
            resume_parameters = self._ReadParameters()
            self._SetMaxIter(resume_parameters, 4)
            resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
                "load_restart"       : true,
                "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
                "restart_load_step"  : 3
            }"""))
            resume_analysis = self._RunToConvergence(resume_parameters)
            resume_algorithm = resume_analysis.GetAlgorithm()
            resume_obj_value = resume_algorithm.GetOptimizedObjectiveValue()
            resume_physical_field = self._GetPhysicalField(resume_analysis)

            # filter + sigmoidal projection + constraint solve give a ~1e-9 relative noise floor
            self.assertAlmostEqual(reference_obj_value, resume_obj_value, delta=abs(reference_obj_value) * 1e-7)
            max_abs_reference_value = max(abs(value) for value in reference_physical_field)
            self.assertVectorAlmostEqual(reference_physical_field, resume_physical_field, places=None, delta=max_abs_reference_value * 1e-6)

            # the control resumes from its checkpointed phi fields, not from ProjectBackward() of the mesh
            reference_control, = reference_algorithm.master_control.GetListOfControls()
            resume_control, = resume_algorithm.master_control.GetListOfControls()
            self.assertVectorAlmostEqual(reference_control.physical_phi_field.data, resume_control.physical_phi_field.data, places=6)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
