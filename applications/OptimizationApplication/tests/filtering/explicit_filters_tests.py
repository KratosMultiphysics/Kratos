import typing
import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.ShapeOptimizationApplication as KratosSOA
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitFilterSolid(kratos_unittest.TestCase):
    ExplicitFilterTypes = typing.Union[KratosOA.NodalExplicitFilter, KratosOA.ConditionExplicitFilter, KratosOA.ElementExplicitFilter]
    EntityContainerTypes = typing.Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]

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

    def test_ConditionExplicitFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part)

    def test_ElementExplicitFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitFilter(self.model_part, "linear", 1000)
        self.__RunConsistencyTest(vm_filter, self.model_part.Elements, self.model_part)

    def __RunConsistencyTest(self, vm_filter: ExplicitFilterTypes, entities: EntityContainerTypes, model_part: Kratos.ModelPart):
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
        filtered_data = vm_filter.FilterField(unfiltered_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_data), math.sqrt(12 * len(entities)), 12)

        integration_weights = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(integration_weights, constant_field_value)
        vm_filter.GetIntegrationWeights(integration_weights)
        integrated_constant_field = integration_weights * unfiltered_field
        filtered_integrated_data = vm_filter.FilterIntegratedField(integrated_constant_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_integrated_data), math.sqrt(12 * len(entities)), 12)

class TestExplicitFilterShell(kratos_unittest.TestCase):
    ExplicitFilterTypes = typing.Union[KratosOA.NodalExplicitFilter, KratosOA.ConditionExplicitFilter, KratosOA.ElementExplicitFilter]
    EntityContainerTypes = typing.Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosSOA.DF1DX_MAPPED)
        cls.model_part.AddNodalSolutionStepVariable(KratosSOA.SHAPE_UPDATE)
        cls.model_part.AddNodalSolutionStepVariable(KratosSOA.CONTROL_POINT_UPDATE)
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("../mdpas/shell", cls.model_part)

    def test_ConsistenceFiltering(self):
        pass

    def test_Reference(self):
        N = 10

        neighbour_nodes_exp = Kratos.Expression.NodalExpression(self.model_part)
        KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(neighbour_nodes_exp)

        filter_settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_filter",
            "filter_function_type"      : "linear",
            "max_nodes_in_filter_radius": 100000,
            "filter_radius_settings":{
                "filter_radius_type": "constant",
                "filter_radius"     : 0.5
            },
            "filtering_boundary_conditions": {
                "damping_type"              : "nearest_entity",
                "damping_function_type"     : "cosine",
                "damped_model_part_settings": {
                    "test.DISPLACEMENT_left_edge" : [true, true, true],
                    "test.DISPLACEMENT_right_edge": [true, true, true]
                }
            }
        }""")
        opt_problem = OptimizationProblem()
        comp_data_view = ComponentDataView("test", opt_problem)

        explicit_vm_filter = FilterFactory(self.model, "test", KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, filter_settings)
        explicit_vm_filter.SetComponentDataView(comp_data_view)
        explicit_vm_filter.Initialize()
        explicit_vm_filter.Check()

        vtu_output = Kratos.VtuOutput(self.model_part)

        step_size = 5e-2
        soa_params = Kratos.Parameters("""{
            "filter_function_type": "linear",
            "filter_radius": 0.5,
            "max_nodes_in_filter_radius": 100
        }""")
        soa_filter = KratosSOA.MapperVertexMorphingMatrixFree(self.model_part, self.model_part, soa_params)
        soa_filter.Initialize()
        for i in range(N):
            explicit_vm_filter.Update()
            soa_filter.Update()

            # compute entity normals
            Kratos.NormalCalculationUtils().CalculateNormalsInElements(self.model_part, Kratos.NORMAL)
            normal_exp = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(normal_exp, Kratos.NORMAL)

            # compute entity sizes
            domain_size_exp = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)

            # compute entity domain size normals
            normal_exp = Kratos.Expression.Utils.Scale(normal_exp, domain_size_exp)

            # compute nodal normals with magnitude
            nodal_normal_exp = Kratos.Expression.NodalExpression(self.model_part)
            KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(nodal_normal_exp, normal_exp, neighbour_nodes_exp)

            # assuming the physical space gradient
            gradient = nodal_normal_exp
            Kratos.Expression.VariableExpressionIO.Write(gradient, Kratos.SHAPE_SENSITIVITY, True)
            soa_filter.InverseMap(Kratos.SHAPE_SENSITIVITY, KratosSOA.DF1DX_MAPPED)


            control_space_gradient = explicit_vm_filter.UnfilterField(gradient)
            update = control_space_gradient * step_size / Kratos.Expression.Utils.NormInf(control_space_gradient)
            Kratos.Expression.VariableExpressionIO.Write(update, KratosSOA.CONTROL_POINT_UPDATE, True)
            soa_filter.Map(KratosSOA.CONTROL_POINT_UPDATE, KratosSOA.SHAPE_UPDATE)

            physical_update = explicit_vm_filter.FilterField(update)

            vtu_output.AddContainerExpression("normal", normal_exp)
            vtu_output.AddContainerExpression("normal", nodal_normal_exp)
            vtu_output.AddContainerExpression("control_space_update", update)
            vtu_output.AddContainerExpression("update", physical_update)
            vtu_output.AddContainerExpression("control_space_gradient", control_space_gradient)
            vtu_output.AddHistoricalVariable(Kratos.SHAPE_SENSITIVITY)
            vtu_output.AddHistoricalVariable(KratosSOA.DF1DX_MAPPED)
            vtu_output.AddHistoricalVariable(KratosSOA.SHAPE_UPDATE)
            vtu_output.AddHistoricalVariable(KratosSOA.CONTROL_POINT_UPDATE)
            vtu_output.PrintOutput(f"output_{i}")

            # apply the shape change
            nodal_position_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.NodalPositionExpressionIO.Read(nodal_position_exp, Kratos.Configuration.Current)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_position_exp + physical_update, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_position_exp + physical_update, Kratos.Configuration.Current)

        explicit_vm_filter.Finalize()

if __name__ == "__main__":
    kratos_unittest.main()