import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

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

        tmp = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(tmp, 2.0)
        vm_filter.SetFilterRadius(tmp.Clone())
        Kratos.Expression.LiteralExpressionIO.SetData(tmp, 1.0)
        vm_filter.SetDampingCoefficients(tmp.Clone())
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

class TestExplicitFilterFactory(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("shell", cls.model_part)

        cls.initial_nodal_pos = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(cls.initial_nodal_pos, Kratos.Configuration.Initial)

    def setUp(self) -> None:
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_pos, Kratos.Configuration.Initial)
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_pos, Kratos.Configuration.Current)

    def __RunTestCase(self, filter_function_type: str, damping_function_type: str, ref_file: str) -> None:
        settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_vertex_morphing",
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
                    "test.DISPLACEMENT_left_edge"  : [ true, false, false],
                    "test.DISPLACEMENT_top_edge"   : [false,  true, false],
                    "test.DISPLACEMENT_right_edge" : [false, false,  true],
                    "test.DISPLACEMENT_bottom_edge": [false,  true,  true]
                }
            }
        }""")
        settings["filter_function_type"].SetString(filter_function_type)
        settings["filtering_boundary_conditions"]["damping_function_type"].SetString(damping_function_type)
        vm_filter = Factory(self.model, "test", KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, settings)
        vm_filter.Initialize()

        vtu_output = Kratos.VtuOutput(self.model_part, binary_output=Kratos.VtuOutput.ASCII, precision=3)
        step_size = 5e-2
        for _ in range(10):
            Kratos.NormalCalculationUtils().CalculateNormals(self.model_part)
            normal_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(normal_exp, Kratos.NORMAL, True)

            update = normal_exp * (step_size / Kratos.Expression.Utils.NormInf(normal_exp))
            filtered_update = vm_filter.FilterField(update)

            # update the mesh
            nodal_coords = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.NodalPositionExpressionIO.Read(nodal_coords, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_coords + filtered_update, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_coords + filtered_update, Kratos.Configuration.Current)

        vtu_output.PrintOutput(f"output_{ref_file}")
        params = Kratos.Parameters("""{
            "reference_file_name"   : "explicit_filter_reference_1.vtu.orig",
            "output_file_name"      : "explicit_filter_reference.vtu",
            "remove_output_file"    : true,
            "comparison_type"       : "deterministic"
        }""")
        params["reference_file_name"].SetString(ref_file)
        params["output_file_name"].SetString(f"output_{ref_file}.vtu")
        CompareTwoFilesCheckProcess(params).Execute()

    def test_FilterCosine(self):
        self.__RunTestCase("cosine", "cosine", "explicit_filter_reference_cosine.vtu.orig")

    def test_FilterLinear(self):
        self.__RunTestCase("linear", "cosine", "explicit_filter_reference_linear.vtu.orig")

    def test_FilterGaussian(self):
        self.__RunTestCase("gaussian", "cosine", "explicit_filter_reference_gaussian.vtu.orig")

    def test_FilterQuartic(self):
        self.__RunTestCase("quartic", "cosine", "explicit_filter_reference_quartic.vtu.orig")

if __name__ == "__main__":
    kratos_unittest.main()