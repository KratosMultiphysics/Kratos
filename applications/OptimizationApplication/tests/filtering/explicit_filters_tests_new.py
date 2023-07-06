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
            ReadModelPart("shell_tri_non_uni", cls.model_part)

    def test_NodalExplicitVertexMorphingFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("shell")
        fixed_model_part = self.model_part.GetSubModelPart("top_edge")
        Kratos.VariableUtils().SetNonHistoricalVariableToZero(Kratos.RADIUS, model_part.Nodes)

        # for condition in model_part.Conditions:
        #     domain_size_contribution = condition.GetGeometry().DomainSize() / 3
        #     for node in condition.GetGeometry():
        #         node.SetValue(Kratos.RADIUS, node.GetValue(Kratos.RADIUS) + domain_size_contribution)

        vtu_output = Kratos.VtuOutput(model_part, True, Kratos.VtuOutput.ASCII, 9)

        vm_filter = KratosOA.NodalExplicitVertexMorphingFilter(model_part, fixed_model_part, "cosine", 1000)

        filter_radius = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(filter_radius, 1.5)
        vm_filter.SetFilterRadius(filter_radius)

        integration_weight_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(integration_weight_field, Kratos.Array3([0.0, 0.0, 0.0]))
        vm_filter.GetIntegrationWeights(integration_weight_field)
        vtu_output.AddContainerExpression("integration_weight", integration_weight_field)


        unfiltered_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, Kratos.Array3([0.0, 0.0, 1.0]))
        unfiltered_field *= integration_weight_field
        vtu_output.AddContainerExpression("unfiltered_field", unfiltered_field)

        filtered_data = vm_filter.FilterIntegratedField(unfiltered_field)

        vtu_output.AddContainerExpression("filtered_field", filtered_data)

        vtu_output.PrintOutput("temp")

        # print(integration_weight_field.Evaluate())
        kl

        constant_field_value = Kratos.Array3([2.0, 2.0, 2.0])

        unfiltered_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, constant_field_value)
        filtered_data = vm_filter.FilterField(unfiltered_field)

        print(KratosOA.ExpressionUtils.NormL2(filtered_data))


if __name__ == "__main__":
    kratos_unittest.main()