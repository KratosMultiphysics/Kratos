import math
import typing
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitFilterConsistency(kratos_unittest.TestCase):
    FilterUtilsType = typing.Union[KratosOA.NodalExplicitFilterUtils, KratosOA.ConditionExplicitFilterUtils, KratosOA.ElementExplicitFilterUtils]
    EntityContainerType = typing.Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

    def test_NodalExplicitFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("design")
        vm_filter = KratosOA.NodalExplicitFilterUtils(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitFilterElementConsistency(self):
        model_part = self.model_part.GetSubModelPart("structure")
        vm_filter = KratosOA.NodalExplicitFilterUtils(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_ConditionExplicitFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitFilterUtils(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part)

    def test_ElementExplicitFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitFilterUtils(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Elements, self.model_part)

    def __RunConsistencyTest(self, vm_filter: FilterUtilsType, entities: EntityContainerType, model_part: Kratos.ModelPart):
        if isinstance(entities, Kratos.NodesArray):
            container_expression_type = Kratos.Expression.NodalExpression
        elif isinstance(entities, Kratos.ConditionsArray):
            container_expression_type = Kratos.Expression.ConditionExpression
        elif isinstance(entities, Kratos.ElementsArray):
            container_expression_type = Kratos.Expression.ElementExpression

        temp_exp = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(temp_exp, 2.0)
        vm_filter.SetFilterRadius(temp_exp.Clone())
        Kratos.Expression.LiteralExpressionIO.SetData(temp_exp, Kratos.Array3([1.0, 1.0, 1.0]))
        vm_filter.SetDampingCoefficients(temp_exp.Clone())
        vm_filter.Update()

        constant_field_value = Kratos.Array3([2.0, 2.0, 2.0])

        unfiltered_field = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, constant_field_value)
        filtered_data = vm_filter.ForwardFilterField(unfiltered_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_data), math.sqrt(12 * len(entities)), 12)

        physical_sensitivity_numpy_array = []
        damping_coeffs_numpy_array = []
        for entity in entities:
            physical_sensitivity_numpy_array.append([entity.Id, entity.Id + 1, entity.Id + 2])
            damping_coeffs_numpy_array.append([(entity.Id % 10) / 10, ((entity.Id + 1) % 10) / 10, ((entity.Id + 2) % 10) / 10])
        physical_sensitivity_numpy_array = numpy.array(physical_sensitivity_numpy_array, dtype=numpy.float64)
        damping_coeffs_numpy_array = numpy.array(damping_coeffs_numpy_array, dtype=numpy.float64)

        physical_sensitivity_field = container_expression_type(model_part)
        damping_coeffs = container_expression_type(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(physical_sensitivity_field, physical_sensitivity_numpy_array)
        Kratos.Expression.CArrayExpressionIO.Read(damping_coeffs, damping_coeffs_numpy_array)
        control_update = physical_sensitivity_field * -1.0
        vm_filter.SetDampingCoefficients(damping_coeffs)
        vm_filter.Update()

        control_sensitivity_field = vm_filter.BackwardFilterField(physical_sensitivity_field)
        physical_update = vm_filter.ForwardFilterField(control_update)

        self.assertAlmostEqual(
            Kratos.Expression.Utils.InnerProduct(physical_sensitivity_field, physical_update),
            Kratos.Expression.Utils.InnerProduct(control_sensitivity_field, control_update), 9)

        integration_weights = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(integration_weights, constant_field_value)
        vm_filter.GetIntegrationWeights(integration_weights)
        integrated_physical_sensitivity_field = physical_sensitivity_field * integration_weights

        temp = vm_filter.BackwardFilterIntegratedField(integrated_physical_sensitivity_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(temp - control_sensitivity_field), 0.0, 9)

if __name__ == "__main__":
    kratos_unittest.main()