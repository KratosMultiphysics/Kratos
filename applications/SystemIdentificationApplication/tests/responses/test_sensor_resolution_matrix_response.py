import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import AddMaskStatusController
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetMaskStatusControllers
from KratosMultiphysics.SystemIdentificationApplication.responses.sensor_resolution_matrix_response import SensorResolutionMatrixResponse

class TestSensorResolutionMatrixResponse(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        """
          (0,0)
            1------2------3------4
            |      |      |      |
            |   1  |  2   |  3   |
            |      |      |      |
            5------6------7------8
            |      |      |      |
            |   4  |  5   |  6   |
            |      |      |      |
            9-----10-----11-----12
                                (6,4)
        """

        cls.model = Kratos.Model()
        cls.domain_model_part = cls.model.CreateModelPart("mask")
        cls.domain_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.domain_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(2, 2.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(3, 4.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(4, 6.0, 0.0, 0.0)
        cls.domain_model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(6, 2.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(7, 4.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(8, 6.0, 2.0, 0.0)
        cls.domain_model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(10, 2.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(11, 4.0, 4.0, 0.0)
        cls.domain_model_part.CreateNewNode(12, 6.0, 4.0, 0.0)

        prop = cls.domain_model_part.CreateNewProperties(1)

        cls.domain_model_part.CreateNewElement("Element2D4N", 1, [1, 2, 6, 5], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 2, [2, 3, 7, 6], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 3, [3, 4, 8, 7], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 4, [5, 6, 10, 9], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 5, [6, 7, 11, 10], prop)
        cls.domain_model_part.CreateNewElement("Element2D4N", 6, [7, 8, 12, 11], prop)

        for node in cls.domain_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [1, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [3, 1, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [3, 3, 0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [5, 3, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.optimization_problem = OptimizationProblem()

        cls.sensor_model_part = cls.model.CreateModelPart("sensors")
        cls.sensors = CreateSensors(cls.sensor_model_part, cls.domain_model_part, parameters)

        cls.sensor_group_data = ComponentDataView("sensors", cls.optimization_problem)

        SetSensors(cls.sensor_group_data, cls.sensors)

        for sensor in cls.sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_STATUS, (sensor.GetNode().Id % 2) / 2)

        elem_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.domain_model_part.Elements, Kratos.PRESSURE)

        elem_ta.data[:] = np.array([1, 1, 0, 1, 1, 0], dtype=np.float64)
        cls.sensors[0].AddTensorAdaptor("mask_exp", elem_ta.Clone())

        elem_ta.data[:] = np.array([1, 1, 0, 1, 0, 1], dtype=np.float64)
        cls.sensors[1].AddTensorAdaptor("mask_exp", elem_ta.Clone())

        elem_ta.data[:] = np.array([0, 0, 0, 1, 0, 0], dtype=np.float64)
        cls.sensors[2].AddTensorAdaptor("mask_exp", elem_ta.Clone())

        elem_ta.data[:] = np.array([0, 0, 1, 0, 1, 1], dtype=np.float64)
        cls.sensors[3].AddTensorAdaptor("mask_exp", elem_ta.Clone())

        # add the mask status controller
        cls.sensor_mask_status = KratosSI.SensorMaskStatus(cls.sensor_model_part, [sensor.GetTensorAdaptor("mask_exp").Clone() for sensor in cls.sensors], 0)
        AddMaskStatusController(cls.sensor_group_data, "mask_exp", cls.sensor_mask_status)

    def test_CalculateValue1(self):
        self.sensor_model_part.GetNode(1).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosSI.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosSI.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1
        """
        self.__UpdateStatusDependents()

        params = Kratos.Parameters("""{
            "sensor_group_name"   : "sensors",
            "sensor_mask_name"    : "mask_exp",
            "step_size"           : 1.5,
            "filter_radius"       : 0.0,
            "filter_function_type": "linear",
            "node_cloud_mesh"     : false,
            "max_items_in_bucket" : 10,
            "store_filter_matrix" : false,
            "echo_level"          : 0
        }""")
        response = SensorResolutionMatrixResponse("test1", self.model, params, self.optimization_problem)
        response.Initialize()

        mask_statuses = np.array(response.mask_status.GetMaskStatuses(), dtype=np.float64)
        result = mask_statuses @ mask_statuses.T
        self.assertAlmostEqual(response.CalculateValue(), 0.5 * np.linalg.norm(1.5 * result - np.identity(mask_statuses.shape[0]), ord="fro") ** 2)

    def test_CalculateValue2(self):
        self.sensor_model_part.GetNode(1).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(2).SetValue(KratosSI.SENSOR_STATUS, 0)
        self.sensor_model_part.GetNode(3).SetValue(KratosSI.SENSOR_STATUS, 1)
        self.sensor_model_part.GetNode(4).SetValue(KratosSI.SENSOR_STATUS, 0)

        """
        status  1   0   1   0
        mask    m1  m2  m3  m4
        e1      1   1   0   0
        e2      1   1   0   0
        e3      0   0   0   1
        e4      1   1   1   0
        e5      1   0   0   1
        e6      0   1   0   1
        """
        self.__UpdateStatusDependents()

        params = Kratos.Parameters("""{
            "sensor_group_name"   : "sensors",
            "sensor_mask_name"    : "mask_exp",
            "step_size"           : 1.5,
            "filter_radius"       : 10.5,
            "filter_function_type": "linear",
            "node_cloud_mesh"     : false,
            "max_items_in_bucket" : 10,
            "store_filter_matrix" : false,
            "echo_level"          : 0
        }""")
        response = SensorResolutionMatrixResponse("test1", self.model, params, self.optimization_problem)
        response.Initialize()
        filter = response.utils.GetFilter()

        mask_statuses = np.array(response.mask_status.GetMaskStatuses(), dtype=np.float64)
        mask_statuses = mask_statuses @ mask_statuses.T
        result = np.empty(mask_statuses.shape, dtype=np.float64)

        for col in range(mask_statuses.shape[0]):
            ta = Kratos.TensorAdaptors.DoubleTensorAdaptor(response.mask_status.GetMaskContainer(), Kratos.DoubleNDData(mask_statuses[col, :], copy=False), copy=False)
            filtered_ta = filter.ForwardFilterField(filter.BackwardFilterField(ta))
            result[:, col] = filtered_ta.data[:]
        self.assertAlmostEqual(response.CalculateValue(), 0.5 * np.linalg.norm(1.5 * result - np.identity(mask_statuses.shape[0]), ord="fro") ** 2)

    def test_CalculateGradient1(self):
        params = Kratos.Parameters("""{
            "sensor_group_name"   : "sensors",
            "sensor_mask_name"    : "mask_exp",
            "step_size"           : 1.5,
            "filter_radius"       : 1e-8,
            "filter_function_type": "linear",
            "node_cloud_mesh"     : false,
            "max_items_in_bucket" : 10,
            "store_filter_matrix" : false,
            "echo_level"          : 0
        }""")
        response = SensorResolutionMatrixResponse("test4", self.model, params, self.optimization_problem)
        response.Initialize()

        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 4) / 4)
        self.__UpdateStatusDependents()

        ref_value = response.CalculateValue()
        analytical_gradient = Kratos.TensorAdaptors.VariableTensorAdaptor(self.sensor_model_part.Nodes, Kratos.PRESSURE)
        cta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([analytical_gradient], perform_collect_data_recursively=False, perform_store_data_recursively=False, copy=False)
        response.CalculateGradient({KratosSI.SENSOR_STATUS: cta})

        masks = np.array(response.mask_status.GetMasks(), dtype=np.float64)

        mask_statuses = np.array(response.mask_status.GetMaskStatuses(), dtype=np.float64)
        mat_R = 1.5 * mask_statuses @ mask_statuses.T
        self.assertAlmostEqual(ref_value, 0.5 * np.linalg.norm(mat_R - np.identity(mat_R.shape[0]), ord="fro") ** 2)

        for i, node in enumerate(self.sensor_model_part.Nodes):
            mat_dR_di = np.outer(masks[:, i], masks[:, i]) * 2.0 * 1.5 * node.GetValue(KratosSI.SENSOR_STATUS)
            value = np.sum((mat_R - np.identity(mat_R.shape[0])) * mat_dR_di)
            self.assertAlmostEqual(value, analytical_gradient.data[i], 5)

    def test_CalculateGradient2(self):
        params = Kratos.Parameters("""{
            "sensor_group_name"   : "sensors",
            "sensor_mask_name"    : "mask_exp",
            "step_size"           : 1.5,
            "filter_radius"       : 10.5,
            "filter_function_type": "linear",
            "node_cloud_mesh"     : false,
            "max_items_in_bucket" : 10,
            "store_filter_matrix" : false,
            "echo_level"          : 0
        }""")
        response = SensorResolutionMatrixResponse("test4", self.model, params, self.optimization_problem)
        response.Initialize()

        for node in self.sensor_model_part.Nodes:
            node.SetValue(KratosSI.SENSOR_STATUS, (node.Id % 4) / 4)
        self.__UpdateStatusDependents()

        ref_value = response.CalculateValue()
        analytical_gradient = Kratos.TensorAdaptors.VariableTensorAdaptor(self.sensor_model_part.Nodes, Kratos.PRESSURE)
        cta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([analytical_gradient], perform_collect_data_recursively=False, perform_store_data_recursively=False, copy=False)
        response.CalculateGradient({KratosSI.SENSOR_STATUS: cta})

        delta = 1e-8
        for i, node in enumerate(self.sensor_model_part.Nodes):
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) + delta)
            self.__UpdateStatusDependents()
            fd_sensitivity = (response.CalculateValue() - ref_value) / delta
            node.SetValue(KratosSI.SENSOR_STATUS, node.GetValue(KratosSI.SENSOR_STATUS) - delta)
            self.assertAlmostEqual(fd_sensitivity, analytical_gradient.data[i], 5)

    def __UpdateStatusDependents(self) -> None:
        for mask_status_controller in GetMaskStatusControllers(self.sensor_group_data, "mask_exp"):
            mask_status_controller.Update()

if __name__ == '__main__':
    UnitTest.main()