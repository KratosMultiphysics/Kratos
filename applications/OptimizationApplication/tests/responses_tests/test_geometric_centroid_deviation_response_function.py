import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.responses.geometric_centroid_deviation_response_function import GeometricCentroidDeviationResponseFunction

class TestGeometricCentroidDeviationResponseFunction(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("../model_part_utils_test/quads", cls.model_part)

        cls.response_function = GeometricCentroidDeviationResponseFunction("geo_centroid", cls.model, Kratos.Parameters("""{"evaluated_model_part_names": ["test"]}"""))
        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def _CheckSensitivity(self, update_method, expression_sensitivity_retrieval_method, delta, precision):
        ref_value = self.response_function.CalculateValue()
        for node in self.model_part.Nodes:
            update_method(node, delta)
            value = self.response_function.CalculateValue()
            sensitivity = (value - ref_value)/delta
            update_method(node, -delta)
            self.assertAlmostEqual(sensitivity, expression_sensitivity_retrieval_method(node), precision)

    def _UpdateNodalPositions(self, direction, entity, delta):
        if direction == 0:
            entity.X += delta
        if direction == 1:
            entity.Y += delta
        if direction == 2:
            entity.Z += delta

    def test_CalculateValue(self):
        self.assertAlmostEqual(self.ref_value, 0.0, 12)

    def test_CalculateShapeSensitivity(self):
        sensitivity = KratosOA.CollectiveExpression([Kratos.Expression.NodalExpression(self.model_part)])
        self.response_function.CalculateGradient({KratosOA.SHAPE: sensitivity})
        Kratos.Expression.VariableExpressionIO.Write(sensitivity.GetContainerExpressions()[0], KratosOA.SHAPE, False)

        # calculate nodal shape sensitivities
        self._CheckSensitivity(
            lambda x, y: self._UpdateNodalPositions(0, x, y),
            lambda x: x.GetValue(KratosOA.SHAPE_X),
            1e-6,
            4)

        self._CheckSensitivity(
            lambda x, y: self._UpdateNodalPositions(1, x, y),
            lambda x: x.GetValue(KratosOA.SHAPE_Y),
            1e-6,
            4)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()