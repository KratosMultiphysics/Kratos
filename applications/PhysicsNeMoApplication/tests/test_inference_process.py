from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestInferenceProcess(KratosUnittest.TestCase):
    """Uses a scripted AFFINE model 2*x + 1 so outputs are exactly
    predictable: input VELOCITY (N,3) -> output ACCELERATION (N,3).

    The offset is deliberate. For a homogeneous-linear stand-in, f(0) = 0,
    a transform that is silently dropped between the model and the Kratos
    write is indistinguishable from one that is applied - which is the
    structural reason a run of de-normalization bugs survived this suite.
    """

    def setUp(self):
        self.checkpoint = Path("test_inference_model.pt")

        class Affine(torch.nn.Module):
            def forward(self, x):
                return 2.0 * x + 1.0

        torch.jit.script(Affine()).save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        for i in range(4):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.VELOCITY, [i + 1.0, 0.0, -1.0])

    def tearDown(self):
        KratosUtilities.DeleteFileIfExisting(str(self.checkpoint))

    def _CreateProcess(self, execution_point="finalize_solution_step", output_interval=1):
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file" : "test_inference_model.pt",
                    "device"          : "cpu"
                },
                "input_fields"    : [ { "variable_name" : "VELOCITY", "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "ACCELERATION", "data_location" : "node_historical" } ],
                "execution_point" : "finalize_solution_step",
                "output_interval" : 1
            }
        }""")
        settings["Parameters"]["execution_point"].SetString(execution_point)
        settings["Parameters"]["output_interval"].SetInt(output_interval)
        return inference_process.Factory(settings, self.model)

    def _AssertPrediction(self):
        for node in self.model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertVectorAlmostEqual(acceleration, [2.0 * v + 1.0 for v in velocity])

    def test_FinalizeSolutionStepExecution(self):
        process = self._CreateProcess("finalize_solution_step")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteInitializeSolutionStep()  # wrong point: must do nothing
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), [0.0, 0.0, 0.0])
        process.ExecuteFinalizeSolutionStep()
        self._AssertPrediction()

    def test_InitializeSolutionStepExecution(self):
        process = self._CreateProcess("initialize_solution_step")
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteInitializeSolutionStep()
        self._AssertPrediction()

    def test_OutputIntervalGating(self):
        process = self._CreateProcess("finalize_solution_step", output_interval=2)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()  # step 1: not due
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), [0.0, 0.0, 0.0])
        self.model_part.ProcessInfo[Kratos.STEP] = 2
        process.ExecuteFinalizeSolutionStep()  # step 2: due
        self._AssertPrediction()

    def test_InvalidExecutionPointRaises(self):
        with self.assertRaisesRegex(ValueError, "Unsupported execution point"):
            self._CreateProcess("before_output")


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestOutputNormalizationThroughProcess(KratosUnittest.TestCase):
    """The card de-normalization at its CALL SITE, not just its maths.

    model_registry's unit tests cover the transform exactly, but nothing
    drove a real process with an "output_normalization" card - so removing
    the call from _WriteOutputs, or flipping scale_only in
    _WriteUncertainty, passed the entire suite. These fail on both.

    The stand-in is AFFINE (2x + 1) on purpose: for a homogeneous-linear
    model a missing affine de-normalization is indistinguishable from a
    present one, which is exactly why this class of bug kept surviving.
    """

    _MEAN, _STD = 100.0, 3.0

    def setUp(self):
        self.checkpoints = []
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.NODAL_PAUX)
        for i in range(4):
            node = self.model_part.CreateNewNode(i + 1, float(i), 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, float(i + 1))

    def tearDown(self):
        for path in self.checkpoints:
            KratosUtilities.DeleteFileIfExisting(str(path))
            KratosUtilities.DeleteFileIfExisting(str(path) + ".card.json")

    def _SaveAffine(self, path, scale, offset, card=True):
        class Affine(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.scale = scale
                self.offset = offset

            def forward(self, x):
                return self.scale * x + self.offset

        torch.jit.script(Affine()).save(str(path))
        if card:
            model_registry.SaveModelCard(path, {
                "output_normalization": {"type": "mean_std",
                                         "mean": [self._MEAN], "std": [self._STD]}})
        self.checkpoints.append(path)
        return path

    def _Pressures(self):
        return numpy.array([node.GetSolutionStepValue(Kratos.PRESSURE)
                            for node in self.model_part.Nodes])

    def _Written(self, variable):
        return numpy.array([node.GetSolutionStepValue(variable)
                            for node in self.model_part.Nodes])

    def test_SingleModelOutputIsDeNormalized(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        checkpoint = self._SaveAffine(Path("test_norm_single.pt"), 2.0, 1.0)
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % checkpoint)
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        # raw prediction 2p + 1, de-normalized as x*std + mean
        expected = (2.0 * self._Pressures() + 1.0) * self._STD + self._MEAN
        numpy.testing.assert_allclose(self._Written(Kratos.TEMPERATURE), expected, rtol=1e-9)

    def test_EnsembleMeanIsShiftedAndSpreadIsNot(self):
        """Two mutations at once.

        The mean must be scaled AND shifted; the standard deviation must be
        scaled and NOT shifted. This also covers the ensemble path, which
        reaches the model through _GetEnsembleModels rather than _GetModel.
        """
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        a = self._SaveAffine(Path("test_norm_ens_a.pt"), 2.0, 1.0)
        self._SaveAffine(Path("test_norm_ens_b.pt"), 4.0, 1.0, card=False)
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : {
                    "checkpoint_file"  : "%s",
                    "checkpoint_files" : ["test_norm_ens_a.pt", "test_norm_ens_b.pt"],
                    "device"           : "cpu"
                },
                "uncertainty"     : {
                    "method"             : "ensemble",
                    "uncertainty_fields" : [ { "variable_name" : "NODAL_PAUX", "data_location" : "node_historical" } ]
                },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % a)
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()

        pressures = self._Pressures()
        members = numpy.stack([2.0 * pressures + 1.0, 4.0 * pressures + 1.0])
        raw_mean = members.mean(axis=0)
        # ddof=1: the ensemble reduction uses the SAMPLE standard deviation,
        # as test_uncertainty_utils::test_EnsembleMeanAndStd also assumes
        raw_std = members.std(axis=0, ddof=1)

        numpy.testing.assert_allclose(
            self._Written(Kratos.TEMPERATURE), raw_mean * self._STD + self._MEAN, rtol=1e-6)
        # scaled, NOT shifted - shifting a spread by the training mean is
        # meaningless, and is the mutation the source comment warns about
        numpy.testing.assert_allclose(
            self._Written(Kratos.NODAL_PAUX), raw_std * self._STD, rtol=1e-6)
        self.assertLess(numpy.abs(self._Written(Kratos.NODAL_PAUX)).max(), self._MEAN)

    def test_InputsAreNormalizedFromTheCard(self):
        """The symmetric half: a model trained on standardized features
        expects standardized features. Both keys on one card, so the
        expected value is the full chain: (2 * (p - m_in) / s_in + 1) *
        s_out + m_out - an affine stand-in again, for the same reason."""
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        checkpoint = self._SaveAffine(Path("test_norm_input.pt"), 2.0, 1.0, card=False)
        in_mean, in_std = 2.5, 0.5
        model_registry.SaveModelCard(checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [in_mean], "std": [in_std]},
            "output_normalization": {"type": "mean_std", "mean": [self._MEAN], "std": [self._STD]}})
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % checkpoint)
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        expected = (2.0 * (self._Pressures() - in_mean) / in_std + 1.0) * self._STD + self._MEAN
        numpy.testing.assert_allclose(self._Written(Kratos.TEMPERATURE), expected, rtol=1e-9)

    def test_WithoutACardNothingIsTransformed(self):
        from KratosMultiphysics.PhysicsNeMoApplication.processes.inference import inference_process
        checkpoint = self._SaveAffine(Path("test_norm_none.pt"), 2.0, 1.0, card=False)
        settings = Kratos.Parameters("""{
            "Parameters": {
                "model_part_name" : "Main",
                "model_settings"  : { "checkpoint_file" : "%s", "device" : "cpu" },
                "input_fields"    : [ { "variable_name" : "PRESSURE",    "data_location" : "node_historical" } ],
                "output_fields"   : [ { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ]
            }
        }""" % checkpoint)
        process = inference_process.Factory(settings, self.model)
        self.model_part.ProcessInfo[Kratos.STEP] = 1
        process.ExecuteFinalizeSolutionStep()
        numpy.testing.assert_allclose(
            self._Written(Kratos.TEMPERATURE), 2.0 * self._Pressures() + 1.0, rtol=1e-9)


if __name__ == '__main__':
    KratosUnittest.main()
