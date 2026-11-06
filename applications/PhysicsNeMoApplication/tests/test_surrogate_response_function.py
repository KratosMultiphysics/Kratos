"""A trained surrogate standing in for a Kratos response function.

The mirror of the co-simulation solver wrapper, and tested the same way: the
value and the gradient are checked against closed forms rather than against
"it ran". The model is an affine map, so J and dJ/dX are known exactly and a
sign, a missing chain rule or a transposed layout cannot hide.

The interface itself is also part of the contract - a driver reaches this
class only through ``ResponseFunctionInterface``, so ``RunCalculation`` is
exercised as a whole, and the refusals a Kratos response makes
(``GetNodalGradient`` for the wrong variable, ``GetElementalGradient`` at
all) are pinned with the message that tells the caller what to do instead.
"""

import os
import shutil
import tempfile
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge

try:
    import torch
    have_torch = True
except ImportError:
    have_torch = False

have_diffusion = kratos_utils.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")

_CASE_DIR = Path(__file__).parent / "adjoint_cases"
# the affine model's coordinate weights, i.e. the exact dJ/dX it implies
_COORDINATE_WEIGHTS = (2.0, -3.0, 5.0)
_FEATURE_WEIGHT = 7.0


@KratosUnittest.skipUnless(have_torch, "Missing required python module: torch.")
class TestSurrogateResponseFunction(KratosUnittest.TestCase):
    """The surrogate gradient mode, against a closed form."""

    def setUp(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import (
            surrogate_response_function)
        self.module = surrogate_response_function

        self.checkpoint = Path("test_surrogate_response_model.pt")

        class Affine(torch.nn.Module):
            """(1, N, 4) [x, y, z, T] -> (1, N, 1), a known linear map."""

            def __init__(self):
                super().__init__()
                self.weights = torch.nn.Parameter(
                    torch.tensor(list(_COORDINATE_WEIGHTS) + [_FEATURE_WEIGHT],
                                 dtype=torch.float64))

            def forward(self, x):
                return (x * self.weights).sum(dim=-1, keepdim=True)

        torch.jit.script(Affine()).save(str(self.checkpoint))

        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        self.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        for i in range(5):
            node = self.model_part.CreateNewNode(i + 1, float(i), 2.0 * i, -float(i))
            node.SetSolutionStepValue(Kratos.TEMPERATURE, 10.0 + i)

    def tearDown(self):
        kratos_utils.DeleteFileIfExisting(str(self.checkpoint))

    def _Settings(self, **overrides):
        settings = Kratos.Parameters("""{
            "model_part_name"       : "Main",
            "model_settings"        : {
                "checkpoint_file" : "test_surrogate_response_model.pt",
                "device"          : "cpu"
            },
            "input_fields"          : [
                { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" }
            ],
            "output_fields"         : [
                { "variable_name" : "PRESSURE", "data_location" : "node_historical" }
            ],
            "objective"             : {
                "type" : "weighted_sum", "variable_name" : "PRESSURE"
            },
            "gradient_mode"         : "surrogate",
            "model_interface"       : "generic",
            "normalize_coordinates" : false
        }""")
        for key, value in overrides.items():
            if isinstance(value, bool):
                settings[key].SetBool(value)
            else:
                settings[key].SetString(value)
        return settings

    def _Create(self, **overrides):
        return self.module.CreateResponseFunction(
            "surrogate", self._Settings(**overrides), self.model)

    def _ExactValue(self):
        total = 0.0
        for node in self.model_part.Nodes:
            total += (_COORDINATE_WEIGHTS[0] * node.X + _COORDINATE_WEIGHTS[1] * node.Y
                      + _COORDINATE_WEIGHTS[2] * node.Z
                      + _FEATURE_WEIGHT * node.GetSolutionStepValue(Kratos.TEMPERATURE))
        return total

    def test_ValueMatchesTheClosedForm(self):
        response = self._Create()
        response.CalculateValue()
        self.assertAlmostEqual(response.GetValue(), self._ExactValue(), places=9)

    def test_PredictionIsWrittenOntoTheModelPart(self):
        # the gradient is taken around the state the surrogate wrote, so that
        # state has to actually be on the model part
        response = self._Create()
        response.CalculateValue()
        for node in self.model_part.Nodes:
            expected = (_COORDINATE_WEIGHTS[0] * node.X + _COORDINATE_WEIGHTS[1] * node.Y
                        + _COORDINATE_WEIGHTS[2] * node.Z
                        + _FEATURE_WEIGHT * node.GetSolutionStepValue(Kratos.TEMPERATURE))
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE), expected,
                                   places=9)

    def test_GradientMatchesTheClosedForm(self):
        response = self._Create()
        response.RunCalculation(calculate_gradient=True)

        gradient = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)
        self.assertEqual(sorted(gradient), [node.Id for node in self.model_part.Nodes])
        for node_id, value in gradient.items():
            with self.subTest(node=node_id):
                for axis in range(3):
                    self.assertAlmostEqual(value[axis], _COORDINATE_WEIGHTS[axis], places=9)

    def test_GradientIsWrittenIntoShapeSensitivity(self):
        # what makes it readable by tooling that knows nothing about this app
        response = self._Create()
        response.RunCalculation(calculate_gradient=True)
        written = adjoint_bridge.ReadSensitivityField(
            self.model_part, Kratos.SHAPE_SENSITIVITY, "node_non_historical")
        for row in range(written.shape[0]):
            for axis in range(3):
                self.assertAlmostEqual(written[row, axis], _COORDINATE_WEIGHTS[axis], places=9)

    def test_NormalizedCoordinatesGetTheChainRule(self):
        # autograd returns dJ/dx_norm; the physical gradient is that over the
        # bounding-box extent, and forgetting the division is a silent error
        # of exactly that factor
        response = self._Create(normalize_coordinates=True)
        response.RunCalculation(calculate_gradient=True)

        coordinates = numpy.array([[node.X, node.Y, node.Z] for node in self.model_part.Nodes])
        extent = coordinates.max(axis=0) - coordinates.min(axis=0)
        gradient = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)
        for node_id, value in gradient.items():
            for axis in range(3):
                with self.subTest(node=node_id, axis=axis):
                    self.assertAlmostEqual(value[axis],
                                           _COORDINATE_WEIGHTS[axis] / extent[axis], places=9)

    def test_OutputNormalizationScalesValueAndGradient(self):
        """A card's "output_normalization" reaches the written field AND
        the autograd gradient: J is a function of the physical prediction,
        so dJ/dX carries the training scale. A de-normalization applied
        only on the write path would leave the gradient wrong by exactly
        std, and the value reported after CalculateGradient inconsistent
        with the field CalculateValue wrote."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = 100.0, 3.0
        model_registry.SaveModelCard(self.checkpoint, {
            "output_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        response = self._Create()
        response.RunCalculation(calculate_gradient=True)

        n_nodes = self.model_part.NumberOfNodes()
        self.assertAlmostEqual(response.GetValue(), std * self._ExactValue() + mean * n_nodes,
                               places=8)
        for node in self.model_part.Nodes:
            raw = (_COORDINATE_WEIGHTS[0] * node.X + _COORDINATE_WEIGHTS[1] * node.Y
                   + _COORDINATE_WEIGHTS[2] * node.Z
                   + _FEATURE_WEIGHT * node.GetSolutionStepValue(Kratos.TEMPERATURE))
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE),
                                   std * raw + mean, places=8)
        gradient = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)
        for node_id, value in gradient.items():
            for axis in range(3):
                with self.subTest(node=node_id, axis=axis):
                    self.assertAlmostEqual(value[axis], std * _COORDINATE_WEIGHTS[axis], places=9)

    def test_InputNormalizationFromTheCardReachesValueAndGradient(self):
        """The field feature is standardized, the coordinates are not: J
        changes by the feature term, dJ/dX (the coordinate weights) does
        not - and the gradient path reports the same value."""
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
        mean, std = 10.0, 2.0
        model_registry.SaveModelCard(self.checkpoint, {
            "input_normalization": {"type": "mean_std", "mean": [mean], "std": [std]}})
        self.addCleanup(kratos_utils.DeleteFileIfExisting, str(self.checkpoint) + ".card.json")

        def Expected(node):
            return (_COORDINATE_WEIGHTS[0] * node.X + _COORDINATE_WEIGHTS[1] * node.Y
                    + _COORDINATE_WEIGHTS[2] * node.Z
                    + _FEATURE_WEIGHT * (node.GetSolutionStepValue(Kratos.TEMPERATURE) - mean) / std)

        expected = sum(Expected(node) for node in self.model_part.Nodes)
        response = self._Create()
        # the VALUE path first, on its own: CalculateGradient recomputes the
        # value through the autograd objective and would mask a value path
        # that forgot the key (a mutation that survived the first draft)
        response.CalculateValue()
        self.assertAlmostEqual(response.GetValue(), expected, places=9)
        for node in self.model_part.Nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(Kratos.PRESSURE), Expected(node),
                                   places=9)
        response.CalculateGradient()
        self.assertAlmostEqual(response.GetValue(), expected, places=9)
        for value in response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY).values():
            for axis in range(3):
                self.assertAlmostEqual(value[axis], _COORDINATE_WEIGHTS[axis], places=9)

    def test_FlatInterfaceIsRefusedForTheSurrogateGradient(self):
        with self.assertRaisesRegex(ValueError, "COORDINATES"):
            self._Create(model_interface="flat")

    def test_ExactModeNeedsDofFields(self):
        with self.assertRaisesRegex(ValueError, "dof_fields"):
            self._Create(gradient_mode="exact", model_interface="flat")

    def test_GetValueBeforeCalculateValueIsRefused(self):
        with self.assertRaisesRegex(RuntimeError, "CalculateValue"):
            self._Create().GetValue()

    def test_GetNodalGradientBeforeCalculateGradientIsRefused(self):
        response = self._Create()
        response.CalculateValue()
        with self.assertRaisesRegex(RuntimeError, "CalculateGradient"):
            response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)

    def test_WrongNodalVariableIsRefused(self):
        response = self._Create()
        response.RunCalculation(calculate_gradient=True)
        with self.assertRaisesRegex(RuntimeError, "sensitivity_variable"):
            response.GetNodalGradient(Kratos.VELOCITY)

    def test_ElementalGradientPointsAtTheRealAdjointStack(self):
        response = self._Create()
        with self.assertRaisesRegex(RuntimeError, "adjoint_bridge"):
            response.GetElementalGradient(Kratos.KratosGlobals.GetVariable("DENSITY"))

    def test_ObjectiveOutsideTheOutputFieldsIsNamed(self):
        settings = self._Settings()
        settings["objective"]["variable_name"].SetString("TEMPERATURE")
        response = self.module.CreateResponseFunction("surrogate", settings, self.model)
        with self.assertRaisesRegex(ValueError, "output_fields"):
            response.CalculateValue()

    def test_UnknownGradientModeIsNamed(self):
        with self.assertRaisesRegex(ValueError, "gradient_mode"):
            self._Create(gradient_mode="telepathy")

    def test_ReachableThroughTheAdjointBridgeModulePath(self):
        # the point of matching Kratos's factory signature: a driver that
        # resolves responses by module path needs no special case
        settings = Kratos.Parameters("""{
            "response_id"     : "surrogate",
            "response_module" :
                "KratosMultiphysics.PhysicsNeMoApplication.deployment.surrogate_response_function"
        }""")
        settings.AddValue("response_settings", self._Settings())
        response = adjoint_bridge.CreateResponseFunction(settings, self.model)
        response.RunCalculation(calculate_gradient=True)
        fields = adjoint_bridge.EvaluateResponse(
            response, self.model_part, run_lifecycle=False)
        self.assertAlmostEqual(fields.value, self._ExactValue(), places=9)
        for row in range(fields.nodal["SHAPE_SENSITIVITY"].shape[0]):
            self.assertAlmostEqual(fields.nodal["SHAPE_SENSITIVITY"][row, 0],
                                   _COORDINATE_WEIGHTS[0], places=9)


@KratosUnittest.skipUnless(have_torch and have_diffusion,
                           "Missing torch or ConvectionDiffusionApplication.")
class TestSurrogateResponseFunctionExactMode(KratosUnittest.TestCase):
    """The FEM-adjoint gradient mode, on a real solved case."""

    def setUp(self):
        from KratosMultiphysics.PhysicsNeMoApplication.deployment import (
            surrogate_response_function)
        self.module = surrogate_response_function

        self.work_dir = Path(tempfile.mkdtemp(prefix="physicsnemo_surrogate_response_"))
        for path in _CASE_DIR.iterdir():
            shutil.copy(path, self.work_dir / path.name)
        self.previous_dir = Path.cwd()
        os.chdir(self.work_dir)

        # identity on TEMPERATURE: the surrogate leaves the solved state
        # alone, so the exact gradient must reproduce a direct adjoint call
        class Identity(torch.nn.Module):
            def forward(self, x):
                return x

        self.checkpoint = Path("identity.pt")
        torch.jit.script(Identity()).save(str(self.checkpoint))

        import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401
        from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import (
            ConvectionDiffusionAnalysis)
        with open("diffusion_primal.json") as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        self.model = Kratos.Model()
        analysis = ConvectionDiffusionAnalysis(self.model, parameters)
        analysis.Initialize()
        analysis.RunSolutionLoop()
        self.model_part = self.model["ThermalModelPart"]
        for node in self.model_part.GetSubModelPart("ImposedTemperature2D_left").Nodes:
            node.Fix(Kratos.TEMPERATURE)

    def tearDown(self):
        os.chdir(self.previous_dir)
        kratos_utils.DeleteDirectoryIfExisting(str(self.work_dir))

    def test_ExactModeReproducesADirectAdjointCall(self):
        from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
        from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils

        settings = Kratos.Parameters("""{
            "model_part_name" : "ThermalModelPart",
            "model_settings"  : { "checkpoint_file" : "identity.pt", "device" : "cpu" },
            "input_fields"    : [
                { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
            "output_fields"   : [
                { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
            "dof_fields"      : [
                { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" } ],
            "objective"       : {
                "type" : "weighted_sum", "variable_name" : "TEMPERATURE",
                "model_part_name" : "HeatFlux2D_right" },
            "gradient_mode"   : "exact",
            "model_interface" : "flat"
        }""")
        response = self.module.CreateResponseFunction("thermal", settings, self.model)
        response.RunCalculation(calculate_gradient=True)

        assembler = differentiable_residual.TangentAssembler(self.model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("TEMPERATURE", "node_historical")])
        weights = adjoint_bridge.MakeObjectiveWeights(Kratos.Parameters("""{
            "type" : "weighted_sum", "variable_name" : "TEMPERATURE",
            "model_part_name" : "HeatFlux2D_right" }"""), self.model_part)
        reference = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dof_map.FieldsToDofVector(weights), fd_step=1e-6)

        self.assertAlmostEqual(response.GetValue(), 870.0, places=8)
        self.assertGreater(numpy.abs(reference).max(), 1.0)
        gradient = response.GetNodalGradient(Kratos.SHAPE_SENSITIVITY)
        for row, node in enumerate(self.model_part.Nodes):
            for axis in range(3):
                with self.subTest(node=node.Id, axis=axis):
                    self.assertAlmostEqual(gradient[node.Id][axis], reference[row, axis],
                                           places=10)


if __name__ == '__main__':
    KratosUnittest.main()
