"""Unit tests for the adjoint bridge's contracts.

Everything here is pure Kratos + numpy: no torch, no compiled optimization
application, and the response functions are stubs. What is being pinned is
the *conversion* - Kratos hands out ``{entity_id: value}`` dicts and this
application works in row-ordered arrays, and the two must line up with every
other gather in the tree. The real adjoint stacks are exercised in
test_adjoint_integration.py.
"""

import numpy

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge


class _StubResponse:
    """The smallest thing satisfying ResponseFunctionInterface's shape."""

    def __init__(self, value=1.0, nodal=None, elemental=None):
        self.value = value
        self.nodal = nodal or {}
        self.elemental = elemental or {}
        self.calls = []

    def Initialize(self):
        self.calls.append("Initialize")

    def InitializeSolutionStep(self):
        self.calls.append("InitializeSolutionStep")

    def CalculateValue(self):
        self.calls.append("CalculateValue")

    def CalculateGradient(self):
        self.calls.append("CalculateGradient")

    def GetValue(self):
        return self.value

    def GetNodalGradient(self, variable):
        if variable.Name() not in self.nodal:
            raise RuntimeError(f"No gradient for {variable.Name()}!")
        return self.nodal[variable.Name()]

    def GetElementalGradient(self, variable):
        if variable.Name() not in self.elemental:
            raise RuntimeError(f"No gradient for {variable.Name()}!")
        return self.elemental[variable.Name()]


def CreateResponseFunction(response_id, response_settings, model):
    """Module-path factory used by the escape-hatch test (this module)."""
    return _StubResponse(value=response_settings["value"].GetDouble())


def _MakeModelPart(model, name="Test"):
    model_part = model.CreateModelPart(name, 1)
    model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
    # ids deliberately out of order and not contiguous: the conversion must
    # key off ids, and the ROWS must follow the Nodes iteration order
    for node_id, x in ((10, 0.0), (3, 1.0), (7, 2.0), (5, 3.0)):
        model_part.CreateNewNode(node_id, x, 0.0, 0.0)
    properties = model_part.CreateNewProperties(1)
    model_part.CreateNewElement("Element2D2N", 100, [10, 3], properties)
    model_part.CreateNewElement("Element2D2N", 42, [7, 5], properties)
    return model_part


class TestAdjointBridgeConversion(KratosUnittest.TestCase):
    """{id: value} dicts -> row-ordered arrays."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = _MakeModelPart(self.model)
        self.node_ids = [node.Id for node in self.model_part.Nodes]

    def test_NodalGradientFollowsNodesOrderNotIdOrder(self):
        gradient = {node_id: [float(node_id), 0.0, -float(node_id)]
                    for node_id in self.node_ids}
        response = _StubResponse(value=2.5, nodal={"SHAPE_SENSITIVITY": gradient})

        fields = adjoint_bridge.EvaluateResponse(response, self.model_part)

        self.assertAlmostEqual(fields.value, 2.5)
        self.assertEqual(list(fields.node_ids), self.node_ids)
        array = fields.nodal["SHAPE_SENSITIVITY"]
        self.assertEqual(array.shape, (len(self.node_ids), 3))
        for row, node_id in enumerate(self.node_ids):
            self.assertAlmostEqual(array[row, 0], float(node_id))
            self.assertAlmostEqual(array[row, 2], -float(node_id))

    def test_DictIterationOrderDoesNotMatter(self):
        # the whole point: a dict carries no order, and the adjoint model part
        # a response iterates need not be the primal one
        forward = {node_id: [float(node_id), 0.0, 0.0] for node_id in self.node_ids}
        reversed_dict = {node_id: forward[node_id] for node_id in reversed(self.node_ids)}

        a = adjoint_bridge.EvaluateResponse(
            _StubResponse(nodal={"SHAPE_SENSITIVITY": forward}), self.model_part)
        b = adjoint_bridge.EvaluateResponse(
            _StubResponse(nodal={"SHAPE_SENSITIVITY": reversed_dict}), self.model_part)
        self.assertTrue(numpy.array_equal(a.nodal["SHAPE_SENSITIVITY"],
                                          b.nodal["SHAPE_SENSITIVITY"]))

    def test_PartialGradientLeavesUnnamedNodesAtZero(self):
        # restricting the sensitivity model part to a design surface is the
        # normal case, not an error
        gradient = {self.node_ids[1]: [1.0, 2.0, 3.0]}
        fields = adjoint_bridge.EvaluateResponse(
            _StubResponse(nodal={"SHAPE_SENSITIVITY": gradient}), self.model_part)
        array = fields.nodal["SHAPE_SENSITIVITY"]
        self.assertAlmostEqual(array[1, 1], 2.0)
        for row in (0, 2, 3):
            self.assertAlmostEqual(float(numpy.abs(array[row]).max()), 0.0)

    def test_ElementalGradientFollowsElementsOrder(self):
        element_ids = [element.Id for element in self.model_part.Elements]
        gradient = {element_id: 0.5 * element_id for element_id in element_ids}
        fields = adjoint_bridge.EvaluateResponse(
            _StubResponse(elemental={"YOUNG_MODULUS": gradient}),
            self.model_part, nodal_variables=(),
            elemental_variables=(Kratos.KratosGlobals.GetVariable("YOUNG_MODULUS"),))
        array = fields.nodal, fields.elemental["YOUNG_MODULUS"]
        self.assertEqual(array[1].shape, (len(element_ids), 1))
        for row, element_id in enumerate(element_ids):
            self.assertAlmostEqual(array[1][row, 0], 0.5 * element_id)

    def test_ForeignIdIsNamed(self):
        gradient = {12345: [1.0, 0.0, 0.0]}
        with self.assertRaisesRegex(RuntimeError, "different mesh"):
            adjoint_bridge.EvaluateResponse(
                _StubResponse(nodal={"SHAPE_SENSITIVITY": gradient}), self.model_part)

    def test_UnsupportedVariableSaysWhatIsAvailable(self):
        response = _StubResponse(nodal={"SHAPE_SENSITIVITY": {}})
        with self.assertRaisesRegex(RuntimeError, "sensitivity_settings"):
            adjoint_bridge.EvaluateResponse(
                response, self.model_part, nodal_variables=(Kratos.VELOCITY,))

    def test_LifecycleIsDrivenOnceAndCanBeSkipped(self):
        response = _StubResponse(nodal={"SHAPE_SENSITIVITY": {}})
        adjoint_bridge.EvaluateResponse(response, self.model_part)
        self.assertEqual(response.calls, ["Initialize", "InitializeSolutionStep",
                                          "CalculateValue", "CalculateGradient"])

        # a second evaluation must not re-Initialize: that restarts both the
        # primal and the adjoint analysis a real response owns
        adjoint_bridge.EvaluateResponse(response, self.model_part, initialize=False)
        self.assertEqual(response.calls.count("Initialize"), 1)

        already = _StubResponse(nodal={"SHAPE_SENSITIVITY": {}})
        adjoint_bridge.EvaluateResponse(already, self.model_part, run_lifecycle=False)
        self.assertEqual(already.calls, [])


class TestAdjointBridgeFactory(KratosUnittest.TestCase):
    """Dispatch to the applications' factories, and the module-path escape."""

    def test_ModulePathFactoryIsUsed(self):
        settings = Kratos.Parameters("""{
            "response_module"   : "test_adjoint_bridge",
            "response_settings" : { "value" : 7.5 }
        }""")
        response = adjoint_bridge.CreateResponseFunction(settings, Kratos.Model())
        self.assertAlmostEqual(response.GetValue(), 7.5)

    def test_UnknownApplicationListsTheKnownOnes(self):
        settings = Kratos.Parameters('{ "response_application" : "sorcery" }')
        with self.assertRaisesRegex(ValueError, "structural_mechanics"):
            adjoint_bridge.CreateResponseFunction(settings, Kratos.Model())

    def test_BothSelectorsIsRefused(self):
        settings = Kratos.Parameters("""{
            "response_application" : "structural_mechanics",
            "response_module"      : "test_adjoint_bridge"
        }""")
        with self.assertRaisesRegex(ValueError, "not both"):
            adjoint_bridge.CreateResponseFunction(settings, Kratos.Model())

    def test_MissingModuleIsNamed(self):
        settings = Kratos.Parameters('{ "response_module" : "no_such_response_module" }')
        with self.assertRaisesRegex(ImportError, "no_such_response_module"):
            adjoint_bridge.CreateResponseFunction(settings, Kratos.Model())

    def test_NonResponseObjectIsRejected(self):
        # numpy exposes CreateResponseFunction? no - but a module that does
        # and returns the wrong thing must be caught rather than failing later
        settings = Kratos.Parameters("""{
            "response_module"   : "test_adjoint_bridge_bad_factory",
            "response_settings" : { }
        }""")
        import sys
        import types
        module = types.ModuleType("test_adjoint_bridge_bad_factory")
        module.CreateResponseFunction = lambda *args: object()
        sys.modules["test_adjoint_bridge_bad_factory"] = module
        try:
            with self.assertRaisesRegex(TypeError, "ResponseFunctionInterface"):
                adjoint_bridge.CreateResponseFunction(settings, Kratos.Model())
        finally:
            del sys.modules["test_adjoint_bridge_bad_factory"]


class TestAdjointBridgeFields(KratosUnittest.TestCase):
    """Reading and writing sensitivity variables directly."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = _MakeModelPart(self.model)

    def test_WriteThenReadRoundTrips(self):
        field = numpy.arange(3.0 * self.model_part.NumberOfNodes()).reshape(-1, 3)
        adjoint_bridge.WriteSensitivityField(self.model_part, field, Kratos.SHAPE_SENSITIVITY)
        read = adjoint_bridge.ReadSensitivityField(
            self.model_part, Kratos.SHAPE_SENSITIVITY, "node_non_historical")
        self.assertTrue(numpy.array_equal(read, field))

        # and the rows really are the Nodes order, checked entity by entity
        for row, node in enumerate(self.model_part.Nodes):
            written = node.GetValue(Kratos.SHAPE_SENSITIVITY)
            for axis in range(3):
                self.assertAlmostEqual(written[axis], field[row, axis])

    def test_ScalarFieldKeepsATrailingAxis(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, float(node.Id))
        read = adjoint_bridge.ReadSensitivityField(
            self.model_part, Kratos.TEMPERATURE, "node_historical")
        self.assertEqual(read.shape, (self.model_part.NumberOfNodes(), 1))


class TestAdjointObjectiveWeights(KratosUnittest.TestCase):
    """The shared objective definition."""

    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = _MakeModelPart(self.model)
        self.sub = self.model_part.CreateSubModelPart("Design")
        self.sub.AddNodes([3, 7])
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [0.0, 0.0, float(node.Id)])
            node.SetSolutionStepValue(Kratos.TEMPERATURE, float(node.Id))
            node.SetValue(Kratos.NODAL_AREA, 2.0)

    def _Weights(self, json_text):
        return adjoint_bridge.MakeObjectiveWeights(
            Kratos.Parameters(json_text), self.model_part)

    def test_TracedNodeProjectsOnDirection(self):
        weights = self._Weights("""{
            "type" : "traced_node", "variable_name" : "DISPLACEMENT",
            "node_id" : 7, "direction" : [0.0, 0.0, 1.0] }""")
        value = adjoint_bridge.EvaluateObjective(
            self.model_part, weights, Kratos.DISPLACEMENT)
        self.assertAlmostEqual(value, 7.0)

    def test_TracedScalarNeedsNoDirection(self):
        weights = self._Weights("""{
            "type" : "traced_node", "variable_name" : "TEMPERATURE", "node_id" : 5 }""")
        self.assertEqual(weights.shape, (self.model_part.NumberOfNodes(), 1))
        self.assertAlmostEqual(
            adjoint_bridge.EvaluateObjective(self.model_part, weights, Kratos.TEMPERATURE), 5.0)

    def test_WeightedSumOverASubModelPart(self):
        weights = self._Weights("""{
            "type" : "weighted_sum", "variable_name" : "TEMPERATURE",
            "model_part_name" : "Design" }""")
        self.assertAlmostEqual(
            adjoint_bridge.EvaluateObjective(self.model_part, weights, Kratos.TEMPERATURE),
            3.0 + 7.0)

    def test_WeightVariableTurnsTheSumIntoAnIntegral(self):
        weights = self._Weights("""{
            "type" : "weighted_sum", "variable_name" : "TEMPERATURE",
            "weight_variable_name" : "NODAL_AREA" }""")
        self.assertAlmostEqual(
            adjoint_bridge.EvaluateObjective(self.model_part, weights, Kratos.TEMPERATURE),
            2.0 * (10 + 3 + 7 + 5))

    def test_UnknownTypeIsNamed(self):
        with self.assertRaisesRegex(ValueError, "weighted_sum"):
            self._Weights('{ "type" : "guesswork", "variable_name" : "TEMPERATURE" }')

    def test_VectorObjectiveVariableIsRefused(self):
        with self.assertRaisesRegex(ValueError, "neither a double"):
            self._Weights("""{
                "type" : "traced_node", "variable_name" : "CAUCHY_STRESS_VECTOR",
                "node_id" : 5 }""")


if __name__ == '__main__':
    KratosUnittest.main()
