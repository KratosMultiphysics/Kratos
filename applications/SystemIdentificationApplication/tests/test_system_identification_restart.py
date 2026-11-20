import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.kratos_utilities import DeleteDirectoryIfExisting
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestSystemIdentificationRestart(kratos_unittest.TestCase):
    """Restart with sensors and a filtered SysId control: sensors are re-derived (not checkpointed),
    and the primal/adjoint model parts come back together with the whole Model."""
    restart_files_path = "Optimization_Restart_test"

    def _ReadParameters(self) -> Kratos.Parameters:
        with open("auxiliary_files/system_identification_restart/optimization_parameters.json", "r") as file_input:
            return Kratos.Parameters(file_input.read())

    def _SetMaxIter(self, parameters: Kratos.Parameters, max_iter: int) -> None:
        parameters["algorithm_settings"]["settings"]["conv_settings"]["max_iter"].SetInt(max_iter)

    def _RunToConvergence(self, parameters: Kratos.Parameters) -> OptimizationAnalysis:
        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)
        analysis.Run()
        return analysis

    def _GetPhysicalField(self, analysis: OptimizationAnalysis) -> 'list[float]':
        """Material field read off the mesh."""
        control, = analysis.GetAlgorithm().master_control.GetListOfControls()
        return list(control.GetPhysicalField().data)

    def _RunReferenceAndResume(self, update_parameters) -> 'tuple[OptimizationAnalysis, OptimizationAnalysis]':
        """Runs iterations 0..2 uninterrupted (reference), and 0..1 with checkpoints followed by a resume
        from checkpoint 2 up to 3. update_parameters(parameters) is applied to all three runs."""
        DeleteDirectoryIfExisting(self.restart_files_path)

        reference_parameters = self._ReadParameters()
        update_parameters(reference_parameters)
        self._SetMaxIter(reference_parameters, 3)
        reference_analysis = self._RunToConvergence(reference_parameters)

        checkpoint_parameters = self._ReadParameters()
        update_parameters(checkpoint_parameters)
        self._SetMaxIter(checkpoint_parameters, 2)
        checkpoint_parameters.AddValue("restart_settings", Kratos.Parameters("""{
            "save_restart"           : true,
            "restart_file_name"      : \"""" + self.restart_files_path + """/restart_<step>.pkl",
            "restart_save_frequency" : 1
        }"""))
        self._RunToConvergence(checkpoint_parameters)

        resume_parameters = self._ReadParameters()
        update_parameters(resume_parameters)
        self._SetMaxIter(resume_parameters, 3)
        resume_parameters.AddValue("restart_settings", Kratos.Parameters("""{
            "load_restart"       : true,
            "restart_file_name"  : \"""" + self.restart_files_path + """/restart_<step>.pkl",
            "restart_load_step"  : 2
        }"""))
        resume_analysis = self._RunToConvergence(resume_parameters)
        return reference_analysis, resume_analysis

    def _CheckResumeMatchesReference(self, reference_analysis: OptimizationAnalysis, resume_analysis: OptimizationAnalysis) -> None:
        self.assertAlmostEqual(reference_analysis.GetAlgorithm().GetOptimizedObjectiveValue(), resume_analysis.GetAlgorithm().GetOptimizedObjectiveValue(), places=9)
        # the adjoint solve uses an iterative (amgcl) solver: compare relative to the field magnitude
        reference_physical_field = self._GetPhysicalField(reference_analysis)
        max_abs_reference_value = max(abs(value) for value in reference_physical_field)
        self.assertVectorAlmostEqual(reference_physical_field, self._GetPhysicalField(resume_analysis), places=None, delta=max_abs_reference_value * 1e-6)

    def test_system_identification_restart(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            self._CheckResumeMatchesReference(*self._RunReferenceAndResume(lambda _: None))

    def test_system_identification_restart_with_active_clamping(self):
        """Upper bound just above the initial YOUNG_MODULUS: phi overshoots 1 and the clamping
        projection maps it to the bound, so the resume must continue from the checkpointed phi."""
        def update_parameters(parameters: Kratos.Parameters) -> None:
            parameters["controls"][0]["settings"]["control_variable_bounds"].SetVector([0.0, 1.02e10])
            parameters["algorithm_settings"]["settings"]["line_search"]["init_step"].SetDouble(0.1)

        with kratos_unittest.WorkFolderScope(".", __file__):
            reference_analysis, resume_analysis = self._RunReferenceAndResume(update_parameters)
            self._CheckResumeMatchesReference(reference_analysis, resume_analysis)

            reference_control, = reference_analysis.GetAlgorithm().master_control.GetListOfControls()
            resume_control, = resume_analysis.GetAlgorithm().master_control.GetListOfControls()
            self.assertGreater(max(reference_control.physical_phi_field.data), 1.0)
            self.assertVectorAlmostEqual(reference_control.physical_phi_field.data, resume_control.physical_phi_field.data, places=6)

    @classmethod
    def tearDownClass(cls) -> None:
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Structure.time")
            DeleteFileIfExisting("AdjointStructure.time")
            DeleteDirectoryIfExisting(cls.restart_files_path)

if __name__ == "__main__":
    kratos_unittest.main()
