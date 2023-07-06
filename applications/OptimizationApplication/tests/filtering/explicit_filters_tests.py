import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitVertexMorphingFilter(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

    def test_NodalExplicitVertexMorphingFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("design")
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.RADIUS, model_part.Nodes)

        for condition in model_part.Conditions:
            domain_size_contribution = condition.GetGeometry().DomainSize() / 3
            for node in condition.GetGeometry():
                node.SetValue(Kratos.RADIUS, node.GetValue(Kratos.RADIUS) + domain_size_contribution)

        vm_filter = KratosOA.NodalExplicitVertexMorphingFilter(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part, lambda x: x.GetValue(Kratos.RADIUS))

    def test_NodalExplicitVertexMorphingFilterElementConsistency(self):
        model_part = self.model_part.GetSubModelPart("structure")
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.RADIUS, model_part.Nodes)

        for element in model_part.Elements:
            domain_size_contribution = element.GetGeometry().DomainSize() / 4
            for node in element.GetGeometry():
                node.SetValue(Kratos.RADIUS, node.GetValue(Kratos.RADIUS) + domain_size_contribution)

        vm_filter = KratosOA.NodalExplicitVertexMorphingFilter(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part, lambda x: x.GetValue(Kratos.RADIUS))

    def test_ConditionExplicitVertexMorphingFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitVertexMorphingFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part, lambda x: x.GetGeometry().DomainSize())

    def test_ElementExplicitVertexMorphingFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitVertexMorphingFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Elements, self.model_part, lambda x: x.GetGeometry().DomainSize())

    def __RunConsistencyTest(self, vm_filter, entities, model_part, domain_size_lambda):
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

        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormL2(filtered_data), math.sqrt(12 * len(entities)), 12)

        for entity in entities:
            domain_size_contribution = domain_size_lambda(entity)
            entity.SetValue(Kratos.VELOCITY, constant_field_value * domain_size_contribution)

        unfiltered_integrated_field = container_expression_type(model_part)
        if isinstance(vm_filter, KratosOA.NodalExplicitVertexMorphingFilter):
            Kratos.Expression.VariableExpressionIO.Read(unfiltered_integrated_field, Kratos.VELOCITY, False)
        else:
            Kratos.Expression.VariableExpressionIO.Read(unfiltered_integrated_field, Kratos.VELOCITY)

        self.assertAlmostEqual(KratosOA.ExpressionUtils.NormInf(filtered_data - vm_filter.FilterIntegratedField(unfiltered_integrated_field)), 0.0, 12)

if __name__ == "__main__":
    kratos_unittest.main()