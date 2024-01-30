import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitFilter(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

    def test_NodalExplicitFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("design")
        vm_filter = KratosOA.NodalExplicitFilter(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitFilterElementConsistency(self):
        model_part = self.model_part.GetSubModelPart("structure")
        vm_filter = KratosOA.NodalExplicitFilter(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitFilterElementConsistencyWithDamping(self):
        model_part = self.model_part.GetSubModelPart("structure")
        fixed_model_part = self.model_part.GetSubModelPart("fixed")
        vm_filter = KratosOA.NodalExplicitFilter(model_part, fixed_model_part, "linear", "sigmoidal", 1000)
        filter_radius =  Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, 2.0)
        vm_filter.SetFilterRadius(filter_radius)

        constant_field_value = Kratos.Array3([2.0, 2.0, 2.0])
        unfiltered_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, constant_field_value)
        filtered_data = vm_filter.FilterField(unfiltered_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_data), 1e-9, 8)

        integration_weights = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(integration_weights, constant_field_value)
        vm_filter.GetIntegrationWeights(integration_weights)
        integrated_constant_field = integration_weights * unfiltered_field
        damped_filtered_integrated_data = vm_filter.FilterIntegratedField(integrated_constant_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(damped_filtered_integrated_data), 1e-9, 8)

    def test_ConditionExplicitFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part)

    def test_ElementExplicitFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Elements, self.model_part)

    def __RunConsistencyTest(self, vm_filter, entities, model_part):
        if isinstance(entities, Kratos.NodesArray):
            container_expression_type = Kratos.Expression.NodalExpression
        elif isinstance(entities, Kratos.ConditionsArray):
            container_expression_type = Kratos.Expression.ConditionExpression
        elif isinstance(entities, Kratos.ElementsArray):
            container_expression_type = Kratos.Expression.ElementExpression

        filter_radius = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, 2.0)
        vm_filter.SetFilterRadius(filter_radius)

        constant_field_value = Kratos.Array3([2.0, 2.0, 2.0])

        unfiltered_field = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, constant_field_value)
        filtered_data = vm_filter.FilterField(unfiltered_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_data), math.sqrt(12 * len(entities)), 12)

        integration_weights = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(integration_weights, constant_field_value)
        vm_filter.GetIntegrationWeights(integration_weights)
        integrated_constant_field = integration_weights * unfiltered_field
        filtered_integrated_data = vm_filter.FilterIntegratedField(integrated_constant_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_integrated_data), math.sqrt(12 * len(entities)), 12)

if __name__ == "__main__":
    kratos_unittest.main()