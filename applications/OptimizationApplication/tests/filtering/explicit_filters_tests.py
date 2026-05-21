import math
import typing
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestExplicitFilterConsistency(kratos_unittest.TestCase):
    FilterUtilsType = typing.Union[KratosOA.NodeExplicitFilterUtils, KratosOA.ConditionExplicitFilterUtils, KratosOA.ElementExplicitFilterUtils]
    DampingType = typing.Union[KratosOA.NodeExplicitDamping, KratosOA.ConditionExplicitDamping, KratosOA.ElementExplicitDamping]
    EntityContainerType = typing.Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]

    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("solid", cls.model_part)

    def test_NodalExplicitFilterConditionConsistency(self):
        model_part = self.model_part.GetSubModelPart("design")
        vm_filter = KratosOA.NodeExplicitFilterUtils(model_part, "linear", 1000, 0)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_NodalExplicitFilterElementConsistency(self):
        model_part = self.model_part.GetSubModelPart("structure")
        vm_filter = KratosOA.NodeExplicitFilterUtils(model_part, "linear", 1000, 0)
        self.__RunConsistencyTest(vm_filter, model_part.Nodes, model_part)

    def test_ConditionExplicitFilterConsistency(self):
        vm_filter = KratosOA.ConditionExplicitFilterUtils(self.model_part, "linear", 1000, 0)
        self.__RunConsistencyTest(vm_filter, self.model_part.Conditions, self.model_part)

    def test_ElementExplicitFilterConsistency(self):
        vm_filter = KratosOA.ElementExplicitFilterUtils(self.model_part, "linear", 1000, 0)
        self.__RunConsistencyTest(vm_filter, self.model_part.Elements, self.model_part)

    def test_NearestEntityDamping(self):
        damping = KratosOA.NearestNodeExplicitDamping(self.model, Kratos.Parameters("""{"damping_function_type": "cosine", "damped_model_part_settings": { "test.fixed": [true] }}"""), 1)
        self.__RunMatrixTest(damping)

    def test_IntegratedNearestEntityDamping(self):
        damping = KratosOA.IntegratedNearestNodeExplicitDamping(self.model, Kratos.Parameters("""{"damping_function_type": "cosine", "damped_model_part_settings": { "test.fixed": [true] }}"""), 1)
        self.__RunMatrixTest(damping)

    def __RunConsistencyTest(self, vm_filter: FilterUtilsType, entities: EntityContainerType, model_part: Kratos.ModelPart):
        if isinstance(entities, Kratos.NodesArray):
            container_expression_type = Kratos.Expression.NodalExpression
            damping_type = KratosOA.NearestNodeExplicitDamping
        elif isinstance(entities, Kratos.ConditionsArray):
            container_expression_type = Kratos.Expression.ConditionExpression
            damping_type = KratosOA.NearestConditionExplicitDamping
        elif isinstance(entities, Kratos.ElementsArray):
            container_expression_type = Kratos.Expression.ElementExpression
            damping_type = KratosOA.NearestElementExplicitDamping

        radius_exp = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(radius_exp, 2.0)
        vm_filter.SetRadius(radius_exp.Clone())
        damping = damping_type(self.model, Kratos.Parameters("""{"damping_function_type": "cosine"}"""), 3)
        damping.SetRadius(radius_exp.Clone())
        vm_filter.SetDamping(damping)

        damping.Update()
        vm_filter.Update()

        constant_field_value = Kratos.Array3([2.0, 2.0, 2.0])

        unfiltered_field = container_expression_type(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(unfiltered_field, constant_field_value)
        filtered_data = vm_filter.ForwardFilterField(unfiltered_field)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(filtered_data), math.sqrt(12 * len(entities)), 12)

        # check with damping
        damping = damping_type(self.model, Kratos.Parameters("""{"damping_function_type": "cosine", "damped_model_part_settings": { "test.fixed": [true, true, false] }}"""), 3)
        damping.SetRadius(radius_exp.Clone())
        vm_filter.SetDamping(damping)
        damping.Update()

        physical_sensitivity_numpy_array = []
        for entity in entities:
            physical_sensitivity_numpy_array.append([entity.Id, entity.Id + 1, entity.Id + 2])
        physical_sensitivity_numpy_array = numpy.array(physical_sensitivity_numpy_array, dtype=numpy.float64)

        physical_sensitivity_field = container_expression_type(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(physical_sensitivity_field, physical_sensitivity_numpy_array)
        control_update = physical_sensitivity_field * -1.0
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

    def __RunMatrixTest(self, damping: DampingType) -> None:
        model_part = self.model_part
        vm_filter = KratosOA.NodeExplicitFilterUtils(model_part, "linear", 1000, 0)

        radius_exp = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(radius_exp, 0.7)
        vm_filter.SetRadius(radius_exp.Clone())


        damping.SetRadius(radius_exp.Clone())
        vm_filter.SetDamping(damping)
        damping.Update()
        vm_filter.Update()

        physical_sensitivity_numpy_array = []
        for entity in model_part.Nodes:
            physical_sensitivity_numpy_array.append(entity.Id)
        physical_sensitivity_numpy_array = numpy.array(physical_sensitivity_numpy_array, dtype=numpy.float64)

        physical_sensitivity_field = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(physical_sensitivity_field, physical_sensitivity_numpy_array)
        control_update = physical_sensitivity_field * -1.0

        # compute nodal area
        Kratos.CalculateNonHistoricalNodalAreaProcess(model_part).Execute()
        nodal_area_exp = Kratos.Expression.NodalExpression(model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_area_exp, Kratos.NODAL_AREA, False)
        nodal_area = nodal_area_exp.Evaluate()

        integrated_physical_sensitivity_field = Kratos.Expression.Utils.Scale(physical_sensitivity_field, nodal_area_exp)

        A = Kratos.Matrix()
        damping.CalculateMatrix(A, 0)
        D = numpy.array(A)
        vm_filter.CalculateMatrix(A)
        A = numpy.array(A)

        p =  control_update.Evaluate()
        df_dx = physical_sensitivity_field.Evaluate()

        # test forward filtering
        x = D @ A @ p
        x_vm = vm_filter.ForwardFilterField(control_update).Evaluate()
        self.assertVectorAlmostEqual(x, x_vm)

        # now test the backward filtering
        df_dp = A.T @ D @ df_dx
        df_dp_vm = vm_filter.BackwardFilterField(physical_sensitivity_field).Evaluate()
        self.assertVectorAlmostEqual(df_dp, df_dp_vm)

        self.assertAlmostEqual(df_dx.dot(x), df_dp.dot(p))

        # now test the integrated backward filtering
        int_df_dx = integrated_physical_sensitivity_field.Evaluate()
        df_dp_1 = A.T @ D @ (int_df_dx / nodal_area)
        df_dp_1_vm = vm_filter.BackwardFilterIntegratedField(integrated_physical_sensitivity_field).Evaluate()
        self.assertVectorAlmostEqual(df_dp_1, df_dp_vm)
        self.assertVectorAlmostEqual(df_dp_1, df_dp_1_vm)

class TestExplicitFilterReference(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        with kratos_unittest.WorkFolderScope(".", __file__):
            ReadModelPart("../mdpas/shell", cls.model_part)

        cls.optimization_problem = OptimizationProblem(0)

        cls.initial_nodal_pos = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(cls.initial_nodal_pos, Kratos.Configuration.Initial)
        cls.filter_data = ComponentDataView("test", cls.optimization_problem)
        cls.filter_data.SetDataBuffer(1)

        cls.vtu_output = Kratos.VtuOutput(cls.model_part, binary_output=Kratos.VtuOutput.ASCII, precision=6)

    def setUp(self) -> None:
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_pos, Kratos.Configuration.Initial)
        Kratos.Expression.NodalPositionExpressionIO.Write(self.initial_nodal_pos, Kratos.Configuration.Current)

    def test_FilterCosine(self):
        self.__RunTestCase("cosine", "cosine", "explicit_filter_reference_cosine.vtu")

    def test_FilterConstant(self):
        self.__RunTestCase("constant", "cosine", "explicit_filter_reference_constant.vtu")

    def test_FilterLinear(self):
        self.__RunTestCase("linear", "cosine", "explicit_filter_reference_linear.vtu")

    def test_FilterGaussian(self):
        self.__RunTestCase("gaussian", "cosine", "explicit_filter_reference_gaussian.vtu")

    def test_FilterQuartic(self):
        self.__RunTestCase("quartic", "cosine", "explicit_filter_reference_quartic.vtu")

    def test_FilterSigmoidal(self):
        self.__RunTestCase("sigmoidal", "cosine", "explicit_filter_reference_sigmoidal.vtu")

    def __RunTestCase(self, filter_function_type: str, damping_function_type: str, ref_file: str) -> None:
        settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_filter",
            "filter_function_type"      : "linear",
            "max_nodes_in_filter_radius": 100000,
            "echo_level"                : 0,
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
        vm_filter = FilterFactory(self.model, "test", KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, settings)
        vm_filter.SetComponentDataView(ComponentDataView("test", self.optimization_problem))
        vm_filter.Initialize()

        nodal_neighbours = Kratos.Expression.NodalExpression(self.model_part)
        KratosOA.ExpressionUtils.ComputeNumberOfNeighbourElements(nodal_neighbours)

        step_size = 5e-2
        for i in range(10):
            Kratos.NormalCalculationUtils().CalculateNormalsInElements(self.model_part, Kratos.NORMAL)
            element_exp = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.VariableExpressionIO.Read(element_exp, Kratos.NORMAL)
            domain_size_exp = Kratos.Expression.ElementExpression(self.model_part)
            Kratos.Expression.DomainSizeExpressionIO.Read(domain_size_exp)
            physical_element_gradient = Kratos.Expression.Utils.Scale(element_exp, domain_size_exp)

            physical_space_gradient = Kratos.Expression.NodalExpression(self.model_part)
            KratosOA.ExpressionUtils.MapContainerVariableToNodalVariable(physical_space_gradient, physical_element_gradient, nodal_neighbours)

            control_space_gradient = vm_filter.BackwardFilterField(physical_space_gradient)
            control_update = control_space_gradient * (step_size / Kratos.Expression.Utils.NormInf(control_space_gradient))
            physical_update = vm_filter.ForwardFilterField(control_update)

            # Purposefully left out for debugging if required.
            # self.vtu_output.AddContainerExpression("physical_space_gradient", physical_space_gradient)
            # self.vtu_output.AddContainerExpression("control_space_gradient", control_space_gradient)
            # self.vtu_output.AddContainerExpression("control_update", control_update)
            # self.vtu_output.AddContainerExpression("physical_update", physical_update)
            # self.vtu_output.AddContainerExpression("damping_coeffs", vm_filter.GetComponentDataView().GetUnBufferedData()["damping_coefficients"])
            # self.vtu_output.PrintOutput(f"output_{i+1}")

            # update the mesh
            nodal_coords = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.NodalPositionExpressionIO.Read(nodal_coords, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_coords + physical_update, Kratos.Configuration.Initial)
            Kratos.Expression.NodalPositionExpressionIO.Write(nodal_coords + physical_update, Kratos.Configuration.Current)

            vm_filter.Update()

        with kratos_unittest.WorkFolderScope(".", __file__):
            self.vtu_output.PrintOutput(f"output_{ref_file[:-4]}")
            params = Kratos.Parameters("""{
                "reference_file_name"   : "explicit_filter_reference_1.vtu.orig",
                "output_file_name"      : "explicit_filter_reference.vtu",
                "remove_output_file"    : true,
                "comparison_type"       : "deterministic"
            }""")
            params["reference_file_name"].SetString(ref_file)
            params["output_file_name"].SetString(f"output_{ref_file}")
            CompareTwoFilesCheckProcess(params).Execute()


if __name__ == "__main__":
    kratos_unittest.main()