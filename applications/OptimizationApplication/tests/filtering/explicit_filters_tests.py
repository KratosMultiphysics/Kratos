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
        radius_exp = Kratos.TensorAdaptors.DoubleTensorAdaptor(entities, Kratos.DoubleNDData([len(entities)]), copy=False)
        if isinstance(entities, Kratos.NodesArray):
            damping_type = KratosOA.NearestNodeExplicitDamping
        elif isinstance(entities, Kratos.ConditionsArray):
            damping_type = KratosOA.NearestConditionExplicitDamping
        elif isinstance(entities, Kratos.ElementsArray):
            damping_type = KratosOA.NearestElementExplicitDamping

        radius_exp.data[:] = 2.0
        vm_filter.SetRadius(radius_exp)
        damping = damping_type(self.model, Kratos.Parameters("""{"damping_function_type": "cosine"}"""), 3)
        damping.SetRadius(radius_exp)
        vm_filter.SetDamping(damping)

        damping.Update()
        vm_filter.Update()

        unfiltered_field = Kratos.TensorAdaptors.DoubleTensorAdaptor(entities, Kratos.DoubleNDData([len(entities), 3], 2.0), copy=False)
        filtered_data = vm_filter.ForwardFilterField(unfiltered_field)
        self.assertAlmostEqual(numpy.linalg.norm(filtered_data.data), math.sqrt(12 * len(entities)), 12)

        # check with damping
        damping = damping_type(self.model, Kratos.Parameters("""{"damping_function_type": "cosine", "damped_model_part_settings": { "test.fixed": [true, true, false] }}"""), 3)
        damping.SetRadius(radius_exp)
        vm_filter.SetDamping(damping)
        damping.Update()

        physical_sensitivity_numpy_array = []
        for entity in entities:
            physical_sensitivity_numpy_array.append([entity.Id, entity.Id + 1, entity.Id + 2])
        physical_sensitivity_numpy_array = numpy.array(physical_sensitivity_numpy_array, dtype=numpy.float64)

        physical_sensitivity_field = Kratos.TensorAdaptors.DoubleTensorAdaptor(entities, Kratos.DoubleNDData(physical_sensitivity_numpy_array), copy=False)
        control_update = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_sensitivity_field)
        control_update.data = physical_sensitivity_field.data * -1.0
        vm_filter.Update()

        control_sensitivity_field = vm_filter.BackwardFilterField(physical_sensitivity_field)
        physical_update = vm_filter.ForwardFilterField(control_update)

        self.assertAlmostEqual(
            numpy.dot(physical_sensitivity_field.data.ravel(), physical_update.data.ravel()),
            numpy.dot(control_sensitivity_field.data.ravel(), control_update.data.ravel()), 9)

        integration_weights = Kratos.TensorAdaptors.DoubleTensorAdaptor(entities, Kratos.DoubleNDData([len(entities), 3], 2.0), copy=False)
        vm_filter.GetIntegrationWeights(integration_weights)
        integrated_physical_sensitivity_field = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_sensitivity_field)
        integrated_physical_sensitivity_field.data = physical_sensitivity_field.data * integration_weights.data

        temp = vm_filter.BackwardFilterIntegratedField(integrated_physical_sensitivity_field)
        self.assertAlmostEqual(numpy.linalg.norm(temp.data - control_sensitivity_field.data), 0.0, 9)

    def __RunMatrixTest(self, damping: DampingType) -> None:
        model_part = self.model_part
        vm_filter = KratosOA.NodeExplicitFilterUtils(model_part, "linear", 1000, 0)

        radius_exp = Kratos.TensorAdaptors.DoubleTensorAdaptor(model_part.Nodes, Kratos.DoubleNDData([len(model_part.Nodes)], 0.7), copy=False)
        vm_filter.SetRadius(radius_exp)

        damping.SetRadius(radius_exp)
        vm_filter.SetDamping(damping)
        damping.Update()
        vm_filter.Update()

        physical_sensitivity_numpy_array = []
        for entity in model_part.Nodes:
            physical_sensitivity_numpy_array.append(entity.Id)
        physical_sensitivity_numpy_array = numpy.array(physical_sensitivity_numpy_array, dtype=numpy.float64)

        physical_sensitivity_field = Kratos.TensorAdaptors.DoubleTensorAdaptor(model_part.Nodes, Kratos.DoubleNDData(physical_sensitivity_numpy_array), copy=False)
        control_update = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_sensitivity_field)
        control_update.data = physical_sensitivity_field.data * -1.0

        # compute nodal area
        Kratos.CalculateNonHistoricalNodalAreaProcess(model_part).Execute()
        nodal_area_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, Kratos.NODAL_AREA)
        nodal_area_ta.CollectData()

        integrated_physical_sensitivity_field = Kratos.TensorAdaptors.DoubleTensorAdaptor(physical_sensitivity_field)
        integrated_physical_sensitivity_field.data[:] = physical_sensitivity_field.data[:] * nodal_area_ta.data

        A = Kratos.Matrix()
        damping.CalculateMatrix(A, 0)
        D = numpy.array(A)
        vm_filter.CalculateMatrix(A)
        A = numpy.array(A)

        p =  control_update.data
        df_dx = physical_sensitivity_field.data

        # test forward filtering
        x = D @ A @ p
        x_vm = vm_filter.ForwardFilterField(control_update).data
        self.assertVectorAlmostEqual(x, x_vm)

        # now test the backward filtering
        df_dp = A.T @ D @ df_dx
        df_dp_vm = vm_filter.BackwardFilterField(physical_sensitivity_field).data
        self.assertVectorAlmostEqual(df_dp, df_dp_vm)

        self.assertAlmostEqual(df_dx.dot(x), df_dp.dot(p))

        # now test the integrated backward filtering
        int_df_dx = integrated_physical_sensitivity_field.data
        df_dp_1 = A.T @ D @ (int_df_dx / nodal_area_ta.data)
        df_dp_1_vm = vm_filter.BackwardFilterIntegratedField(integrated_physical_sensitivity_field).data
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

        cls.initial_nodal_pos = Kratos.TensorAdaptors.NodePositionTensorAdaptor(cls.model_part.Nodes, Kratos.Configuration.Initial)
        cls.initial_nodal_pos.CollectData()

        cls.filter_data = ComponentDataView("test", cls.optimization_problem)
        cls.filter_data.SetDataBuffer(1)

        cls.vtu_output = Kratos.VtuOutput(cls.model_part, binary_output=Kratos.VtuOutput.ASCII, precision=6)

    def setUp(self) -> None:
        Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.initial_nodal_pos, Kratos.Configuration.Initial, copy=False).StoreData()
        Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.initial_nodal_pos, Kratos.Configuration.Current, copy=False).StoreData()

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

    def test_FilterSigmoidalNodeCloudMesh(self):
        self.__RunTestCase("sigmoidal", "cosine", "explicit_filter_reference_sigmoidal_cloud_mesh.vtu", True)

    def __RunTestCase(self, filter_function_type: str, damping_function_type: str, ref_file: str, node_cloud_mesh=False) -> None:
        settings = Kratos.Parameters("""{
            "filter_type"               : "explicit_filter",
            "filter_function_type"      : "linear",
            "max_nodes_in_filter_radius": 100000,
            "node_cloud_mesh": false,
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
        settings["node_cloud_mesh"].SetBool(node_cloud_mesh)
        vm_filter = FilterFactory(self.model, "test", KratosOA.SHAPE, Kratos.Globals.DataLocation.NodeHistorical, settings)
        vm_filter.SetComponentDataView(ComponentDataView("test", self.optimization_problem))
        vm_filter.Initialize()

        step_size = 5e-2
        for i in range(10):
            Kratos.NormalCalculationUtils().CalculateNormalsInElements(self.model_part, Kratos.NORMAL)
            element_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.NORMAL)
            element_ta.CollectData()

            domain_size_ta = Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor(self.model_part.Elements, Kratos.TensorAdaptors.GeometryMetricsTensorAdaptor.DomainSize)
            domain_size_ta.CollectData()

            physical_element_gradient = Kratos.TensorAdaptors.DoubleTensorAdaptor(element_ta)
            physical_element_gradient.data = element_ta.data * domain_size_ta.data[:, None] # row wise scaling
            physical_space_gradient = KratosOA.OptimizationUtils.MapContainerDataToNodalData(physical_element_gradient, self.model_part.Nodes)

            control_space_gradient = vm_filter.BackwardFilterField(physical_space_gradient)
            control_update = Kratos.TensorAdaptors.DoubleTensorAdaptor(control_space_gradient)
            control_update.data = control_space_gradient.data * (step_size / numpy.max(numpy.abs(control_space_gradient.data)))
            physical_update = vm_filter.ForwardFilterField(control_update)

            # Purposefully left out for debugging if required.
            # self.vtu_output.AddContainerExpression("physical_space_gradient", physical_space_gradient)
            # self.vtu_output.AddContainerExpression("control_space_gradient", control_space_gradient)
            # self.vtu_output.AddContainerExpression("control_update", control_update)
            # self.vtu_output.AddContainerExpression("physical_update", physical_update)
            # self.vtu_output.AddContainerExpression("damping_coeffs", vm_filter.GetComponentDataView().GetUnBufferedData()["damping_coefficients"])
            # self.vtu_output.PrintOutput(f"output_{i+1}")

            # update the mesh
            nodal_coords_initial = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Initial)
            nodal_coords_initial.CollectData()
            nodal_coords_initial.data += physical_update.data
            nodal_coords_initial.StoreData()

            nodal_coords_current = Kratos.TensorAdaptors.NodePositionTensorAdaptor(self.model_part.Nodes, Kratos.Configuration.Current)
            nodal_coords_current.CollectData()
            nodal_coords_current.data += physical_update.data
            nodal_coords_current.StoreData()

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