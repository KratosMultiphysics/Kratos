import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitSmoothDampingFilter(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

    def test_NodalExplicitSmoothDampingFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("design")
        vm_filter = KratosOA.NodalExplicitSmoothDampingFilterUtils(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitSmoothDampingFilterElementConsistency(self):
        model_part = self.model_part.GetSubModelPart("structure")
        vm_filter = KratosOA.NodalExplicitSmoothDampingFilterUtils(model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitSmoothDampingFilterElementConsistencyWithDamping(self):
        model_part = self.model_part.GetSubModelPart("structure")
        fixed_model_part = self.model_part.GetSubModelPart("fixed")
        vm_filter = KratosOA.NodalExplicitSmoothDampingFilterUtils(model_part, fixed_model_part, "linear", "sigmoidal", 1000)
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

    def test_ConditionExplicitSmoothDampingFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitSmoothDampingFilterUtils(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part)

    def test_ElementExplicitSmoothDampingFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitSmoothDampingFilterUtils(self.model_part, "linear", 1000)
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

    def test_Factory(self) -> None:
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("../mdpas/shell", model_part)

        optimization_problem = OptimizationProblem(0)
        filter_data = ComponentDataView("test", optimization_problem)
        filter_data.SetDataBuffer(1)
        settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_smooth_damping_filter",
            "filter_function_type"      : "linear",
            "damping_function_type"     : "cosine",
            "max_nodes_in_filter_radius": 100000,
            "echo_level"                : 0,
            "filter_radius_settings":{
                "filter_radius_type": "constant",
                "filter_radius"     : 0.5
            },
            "filtering_boundary_conditions": {
                "test.DISPLACEMENT_left_edge"  : [ true, false, false]
            }
        }""")
        vm_filter = FilterFactory(model, "test", KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, settings)
        vm_filter.SetComponentDataView(ComponentDataView("test", optimization_problem))
        vm_filter.Initialize()

        nodal_neighbours = Kratos.Expression.NodalExpression(model_part)
        KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(nodal_neighbours)

        Kratos.NormalCalculationUtils().CalculateNormalsInElements(model_part, Kratos.NORMAL)
        element_exp = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.VariableExpressionIO.Read(element_exp, Kratos.NORMAL)
        domain_size_exp = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)
        physical_element_gradient = Kratos.Expression.Utils.Scale(element_exp, domain_size_exp)

        physical_space_gradient = Kratos.Expression.NodalExpression(model_part)
        KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(physical_space_gradient, physical_element_gradient, nodal_neighbours)

        control_space_gradient = vm_filter.BackwardFilterField(physical_space_gradient)
        control_update = control_space_gradient * (1e-2 / Kratos.Expression.Utils.NormInf(control_space_gradient))
        physical_update = vm_filter.ForwardFilterField(control_update)

        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(control_space_gradient), 0.02268505190702854)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(physical_update), 0.010249567403956243)

if __name__ == "__main__":
    kratos_unittest.main()