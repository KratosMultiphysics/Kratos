"""The adjoint bridge driven against Kratos's real adjoint stacks.

test_adjoint_bridge.py pins the conversion against stubs; this pins it
against the thing itself, on two applications that implement the same core
``ResponseFunctionInterface`` with entirely different machinery:

- **StructuralMechanics** - ``AdjointFiniteDifferencing*`` elements, a
  separate ``Kratos.Model`` for the adjoint part, ``SensitivityBuilder``
  producing SHAPE_SENSITIVITY.
- **ConvectionDiffusion** - ``AdjointDiffusionElement``, primal and adjoint
  in the same model, a ``LocalTemperatureAverageResponseFunction``.

In both, the bridge's array is compared node by node against
``sensitivity_utils.ComputeShapeSensitivityField`` - this application's own,
independently implemented adjoint. Two implementations agreeing across two
physics is a much stronger statement than either against finite differences,
and it is what makes the bridge trustworthy as a source of *training data*.

The cases live in ``adjoint_cases/`` and the tests run in a throwaway copy of
it: a Kratos response function reads its primal settings from a file by name
and re-reads the mdpa for its own adjoint part, so it only works from the
case's own working directory.
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
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils
from KratosMultiphysics.PhysicsNeMoApplication.processes import adjoint_sensitivity_process

have_structural = kratos_utils.CheckIfApplicationsAvailable(
    "StructuralMechanicsApplication", "LinearSolversApplication",
    "ConstitutiveLawsApplication")
have_diffusion = kratos_utils.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")

_CASE_DIR = Path(__file__).parent / "adjoint_cases"
_TIP_NODE = 9


class _CaseDirectory:
    """Mixin running each test in a throwaway copy of adjoint_cases/.

    Deliberately NOT a TestCase subclass: the registration guard in
    test_suite_registration.py flags every TestCase class that no suite
    adds, and an abstract base would be a false positive that teaches
    people to ignore it.
    """

    def setUp(self):
        self.work_dir = Path(tempfile.mkdtemp(prefix="physicsnemo_adjoint_bridge_"))
        for path in _CASE_DIR.iterdir():
            shutil.copy(path, self.work_dir / path.name)
        self.previous_dir = Path.cwd()
        os.chdir(self.work_dir)

    def tearDown(self):
        os.chdir(self.previous_dir)
        kratos_utils.DeleteDirectoryIfExisting(str(self.work_dir))

    @staticmethod
    def _NodeRows(model_part):
        return {node.Id: row for row, node in enumerate(model_part.Nodes)}


@KratosUnittest.skipUnless(
    have_structural, "Missing StructuralMechanics/ConstitutiveLaws/LinearSolvers applications.")
class TestAdjointBridgeStructuralMechanics(_CaseDirectory, KratosUnittest.TestCase):
    """The cantilever, through StructuralMechanics' own adjoint stack."""

    _RESPONSE = """{
        "response_id"          : "tip_displacement",
        "response_application" : "structural_mechanics",
        "response_settings"    : {
            "response_type"                    : "adjoint_nodal_displacement",
            "gradient_mode"                    : "semi_analytic",
            "step_size"                        : 1e-7,
            "primal_settings"                  : "cantilever_primal.json",
            "adjoint_settings"                 : "auto",
            "primal_data_transfer_with_python" : true,
            "response_part_name"               : "Tip",
            "direction"                        : [0.0, 0.0, 1.0],
            "traced_dof"                       : "DISPLACEMENT",
            "sensitivity_settings"             : {
                "sensitivity_model_part_name"               : "Design",
                "nodal_solution_step_sensitivity_variables" : ["SHAPE_SENSITIVITY"],
                "build_mode"                                : "static"
            }
        }
    }"""

    def _EvaluateThroughTheBridge(self):
        model = Kratos.Model()
        response = adjoint_bridge.CreateResponseFunction(
            Kratos.Parameters(self._RESPONSE), model)
        # the primal part exists as soon as the response is constructed, and
        # it is the one whose Nodes order defines the rows
        model_part = model["Structure"]
        return adjoint_bridge.EvaluateResponse(response, model_part), model_part

    @staticmethod
    def _ShippedField():
        """dJ/dX from this application's own adjoint, for the same objective."""
        from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
            StructuralMechanicsAnalysis)

        with open("cantilever_primal.json") as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        model = Kratos.Model()
        analysis = StructuralMechanicsAnalysis(model, parameters)
        analysis.Initialize()
        analysis.RunSolutionLoop()
        model_part = model["Structure"]
        # the BC process releases its DOFs in ExecuteFinalizeSolutionStep, and
        # the adjoint needs the CONSTRAINED operator - re-fix them
        for node in model_part.GetSubModelPart("Support").Nodes:
            for component in (Kratos.DISPLACEMENT_X, Kratos.DISPLACEMENT_Y,
                              Kratos.DISPLACEMENT_Z):
                node.Fix(component)

        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("DISPLACEMENT", "node_historical")])
        weights = adjoint_bridge.MakeObjectiveWeights(Kratos.Parameters("""{
            "type" : "traced_node", "variable_name" : "DISPLACEMENT",
            "node_id" : %d, "direction" : [0.0, 0.0, 1.0] }""" % _TIP_NODE), model_part)
        value = adjoint_bridge.EvaluateObjective(model_part, weights, Kratos.DISPLACEMENT)
        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dof_map.FieldsToDofVector(weights), fd_step=1e-6)
        return value, field, model_part

    def test_BridgeGradientMatchesTheShippedAdjoint(self):
        fields, bridge_part = self._EvaluateThroughTheBridge()
        value, shipped, shipped_part = self._ShippedField()

        self.assertEqual(list(fields.node_ids), [node.Id for node in shipped_part.Nodes])
        self.assertAlmostEqual(fields.value, value, places=12)

        bridge = fields.nodal["SHAPE_SENSITIVITY"]
        self.assertEqual(bridge.shape, shipped.shape)
        self.assertGreater(numpy.abs(shipped).max(), 1e-9)  # a real gradient
        for row in range(bridge.shape[0]):
            for axis in range(3):
                with self.subTest(node=int(fields.node_ids[row]), axis=axis):
                    # Kratos differentiates FORWARD, the shipped field
                    # centrally, so Kratos carries the larger step error and
                    # sets the tolerance
                    self.assertAlmostEqual(
                        bridge[row, axis], shipped[row, axis],
                        delta=1e-4 * abs(shipped[row, axis]) + 1e-12)

    def test_ObjectiveWeightsReproduceTheResponseValue(self):
        # the two sides must be computing the SAME J; a normalization
        # difference is exactly what makes two correct adjoints disagree
        fields, _ = self._EvaluateThroughTheBridge()
        value, _, model_part = self._ShippedField()
        self.assertAlmostEqual(fields.value, value, places=12)
        self.assertAlmostEqual(
            value, model_part.GetNode(_TIP_NODE).GetSolutionStepValue(Kratos.DISPLACEMENT_Z),
            places=14)

    def test_ElementalGradientIsRefusedWithTheAcceptedList(self):
        model = Kratos.Model()
        response = adjoint_bridge.CreateResponseFunction(
            Kratos.Parameters(self._RESPONSE), model)
        adjoint_bridge.EvaluateResponse(response, model["Structure"])
        with self.assertRaisesRegex(RuntimeError, "sensitivity_settings"):
            adjoint_bridge.EvaluateResponse(
                response, model["Structure"], nodal_variables=(), run_lifecycle=False,
                elemental_variables=(Kratos.KratosGlobals.GetVariable("DENSITY"),))


@KratosUnittest.skipUnless(have_diffusion, "Missing ConvectionDiffusionApplication.")
class TestAdjointBridgeConvectionDiffusion(_CaseDirectory, KratosUnittest.TestCase):
    """A second physics, with an entirely different adjoint implementation."""

    _RESPONSE_PART = "HeatFlux2D_right"

    def _EvaluateThroughTheBridge(self):
        import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401 (registers vars)

        with open("diffusion_response.json") as parameter_file:
            response_settings = Kratos.Parameters(parameter_file.read())["response_settings"]
        settings = Kratos.Parameters("""{
            "response_id"          : "point_temperature",
            "response_application" : "convection_diffusion"
        }""")
        settings.AddValue("response_settings", response_settings)

        model = Kratos.Model()
        response = adjoint_bridge.CreateResponseFunction(settings, model)
        model_part = model["ThermalModelPart"]
        return adjoint_bridge.EvaluateResponse(response, model_part), model_part

    @staticmethod
    def _ShippedField():
        import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401
        from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import (
            ConvectionDiffusionAnalysis)

        with open("diffusion_primal.json") as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        model = Kratos.Model()
        analysis = ConvectionDiffusionAnalysis(model, parameters)
        analysis.Initialize()
        analysis.RunSolutionLoop()
        model_part = model["ThermalModelPart"]
        for node in model_part.GetSubModelPart("ImposedTemperature2D_left").Nodes:
            node.Fix(Kratos.TEMPERATURE)

        assembler = differentiable_residual.TangentAssembler(model_part)
        dof_map = differentiable_residual.DofFieldMap(
            assembler, [("TEMPERATURE", "node_historical")])
        weights = adjoint_bridge.MakeObjectiveWeights(Kratos.Parameters("""{
            "type" : "weighted_sum", "variable_name" : "TEMPERATURE",
            "model_part_name" : "HeatFlux2D_right" }"""), model_part)
        value = adjoint_bridge.EvaluateObjective(model_part, weights, Kratos.TEMPERATURE)
        field = sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dof_map.FieldsToDofVector(weights), fd_step=1e-6)
        return value, field, model_part

    def test_BridgeGradientMatchesTheShippedAdjointUpToTheObjectiveNormalization(self):
        fields, _ = self._EvaluateThroughTheBridge()
        value, shipped, model_part = self._ShippedField()

        # Kratos's "point_temperature" response AVERAGES the temperature over
        # the traced part; MakeObjectiveWeights' "weighted_sum" sums it. The
        # ratio is the node count, and pinning it is the point: two adjoints
        # of two different objectives would "disagree" for a reason that has
        # nothing to do with either being wrong.
        n_traced = model_part.GetSubModelPart(self._RESPONSE_PART).NumberOfNodes()
        self.assertEqual(n_traced, 3)
        self.assertAlmostEqual(fields.value, value / n_traced, places=10)

        bridge = fields.nodal["SHAPE_SENSITIVITY"]
        scaled = shipped / n_traced
        self.assertGreater(numpy.abs(scaled).max(), 1e-6)
        for row in range(bridge.shape[0]):
            for axis in range(3):
                with self.subTest(node=int(fields.node_ids[row]), axis=axis):
                    self.assertAlmostEqual(
                        bridge[row, axis], scaled[row, axis],
                        delta=1e-6 * abs(scaled[row, axis]) + 1e-7)

    def test_OutOfPlaneGradientIsExactlyZero(self):
        # a 2-D case: the Z row must not merely be small
        fields, _ = self._EvaluateThroughTheBridge()
        self.assertEqual(float(numpy.abs(fields.nodal["SHAPE_SENSITIVITY"][:, 2]).max()), 0.0)


@KratosUnittest.skipUnless(have_diffusion, "Missing ConvectionDiffusionApplication.")
class TestAdjointSensitivityProcess(_CaseDirectory, KratosUnittest.TestCase):
    """The process that puts the field on the model part."""

    def _SolvedThermalCase(self):
        import KratosMultiphysics.ConvectionDiffusionApplication  # noqa: F401
        from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import (
            ConvectionDiffusionAnalysis)

        with open("diffusion_primal.json") as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())
        model = Kratos.Model()
        analysis = ConvectionDiffusionAnalysis(model, parameters)
        analysis.Initialize()
        analysis.RunSolutionLoop()
        model_part = model["ThermalModelPart"]
        for node in model_part.GetSubModelPart("ImposedTemperature2D_left").Nodes:
            node.Fix(Kratos.TEMPERATURE)
        return model, model_part

    @staticmethod
    def _ProcessSettings(**overrides):
        settings = Kratos.Parameters("""{
            "Parameters" : {
                "model_part_name"    : "ThermalModelPart",
                "sensitivity_source" : "shipped",
                "dof_fields"         : [
                    { "variable_name" : "TEMPERATURE", "data_location" : "node_historical" }
                ],
                "objective"          : {
                    "type"            : "weighted_sum",
                    "variable_name"   : "TEMPERATURE",
                    "model_part_name" : "HeatFlux2D_right"
                }
            }
        }""")
        for key, value in overrides.items():
            settings["Parameters"].AddString(key, value)
        return settings

    def test_ProcessWritesTheFieldIntoTheVariable(self):
        model, model_part = self._SolvedThermalCase()
        process = adjoint_sensitivity_process.Factory(self._ProcessSettings(), model)
        process.Execute()

        self.assertAlmostEqual(process.last_value, 870.0, places=8)
        written = adjoint_bridge.ReadSensitivityField(
            model_part, Kratos.SHAPE_SENSITIVITY, "node_non_historical")
        self.assertTrue(numpy.array_equal(written, process.last_field))
        self.assertGreater(numpy.abs(written).max(), 1.0)

    def test_ExportedSampleCarriesTheGradient(self):
        # the whole reason the process writes a VARIABLE rather than a file:
        # the ordinary exporter then carries the gradient with no new code,
        # under the key the training path reads
        from KratosMultiphysics.PhysicsNeMoApplication.training.streaming_dataset import (
            GatherSampleArrays)

        model, model_part = self._SolvedThermalCase()
        adjoint_sensitivity_process.Factory(self._ProcessSettings(), model).Execute()

        arrays = GatherSampleArrays(
            model_part, [("TEMPERATURE", "node_historical"),
                         ("SHAPE_SENSITIVITY", "node_non_historical")])
        self.assertIn("SHAPE_SENSITIVITY__node_non_historical", arrays)
        self.assertEqual(arrays["SHAPE_SENSITIVITY__node_non_historical"].shape,
                         (model_part.NumberOfNodes(), 3))

    def test_DesignRestrictionZeroesTheRest(self):
        model, model_part = self._SolvedThermalCase()
        full = adjoint_sensitivity_process.Factory(self._ProcessSettings(), model)
        full.Execute()

        restricted = adjoint_sensitivity_process.Factory(
            self._ProcessSettings(design_sub_model_part_name="HeatFlux2D_right"), model)
        restricted.Execute()

        design_ids = {node.Id
                      for node in model_part.GetSubModelPart("HeatFlux2D_right").Nodes}
        for row, node in enumerate(model_part.Nodes):
            if node.Id in design_ids:
                for axis in range(3):
                    self.assertAlmostEqual(restricted.last_field[row, axis],
                                           full.last_field[row, axis],
                                           delta=1e-9 * abs(full.last_field[row, axis]) + 1e-18)
            else:
                self.assertEqual(float(numpy.abs(restricted.last_field[row]).max()), 0.0)

    def test_ObjectiveOutsideTheDofFieldsIsNamed(self):
        model, _ = self._SolvedThermalCase()
        settings = self._ProcessSettings()
        settings["Parameters"]["objective"]["variable_name"].SetString("VELOCITY")
        process = adjoint_sensitivity_process.Factory(settings, model)
        with self.assertRaisesRegex(ValueError, "dof_fields"):
            process.Execute()

    def test_UnknownSourceIsNamed(self):
        model, _ = self._SolvedThermalCase()
        with self.assertRaisesRegex(ValueError, "response_function"):
            adjoint_sensitivity_process.Factory(
                self._ProcessSettings(sensitivity_source="clairvoyance"), model)

    def test_ResponseFunctionSourceWritesTheBridgeField(self):
        model, model_part = self._SolvedThermalCase()

        with open("diffusion_response.json") as parameter_file:
            response_settings = Kratos.Parameters(parameter_file.read())["response_settings"]
        settings = Kratos.Parameters("""{
            "Parameters" : {
                "model_part_name"    : "ThermalModelPart",
                "sensitivity_source" : "response_function",
                "response_settings"  : {
                    "response_id"          : "point_temperature",
                    "response_application" : "convection_diffusion"
                }
            }
        }""")
        settings["Parameters"]["response_settings"].AddValue(
            "response_settings", response_settings)

        process = adjoint_sensitivity_process.Factory(settings, model)
        process.Execute()

        written = adjoint_bridge.ReadSensitivityField(
            model_part, Kratos.SHAPE_SENSITIVITY, "node_non_historical")
        self.assertTrue(numpy.array_equal(written, process.last_field))
        self.assertAlmostEqual(process.last_value, 290.0, places=8)


if __name__ == '__main__':
    KratosUnittest.main()
