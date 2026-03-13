import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import SetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class TestSensorUtils(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop).SetValue(Kratos.PRESSURE, 10.0)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop).SetValue(Kratos.PRESSURE, 100.0)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])

        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)

        nodal_ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(cls.model_part.Nodes, Kratos.DISPLACEMENT)
        nodal_ta.CollectData()

        elemental_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.model_part.Elements, Kratos.PRESSURE)
        elemental_ta.CollectData()

        for sensor in cls.sensors:
            ta = nodal_ta.Clone()
            ta.data[:] *= sensor.GetNode().Id
            sensor.AddTensorAdaptor("nodal_ta", ta)

            ta = elemental_ta.Clone()
            ta.data[:] *= sensor.GetNode(). Id ** 2
            sensor.AddTensorAdaptor("element_ta", ta)

    def test_CreateSensorView(self):
        n_nodes = self.model_part.NumberOfNodes()
        for sensor in self.sensors:
            nodal_sensor_view = KratosSI.Sensors.SensorView(sensor, "nodal_ta")
            self.assertAlmostEqual((3 * n_nodes * (n_nodes + 1) / 2 + 3 * n_nodes) * sensor.GetNode().Id, np.sum(nodal_sensor_view.GetTensorAdaptor().data))

            element_sensor_view = KratosSI.Sensors.SensorView(sensor, "element_ta")
            self.assertAlmostEqual(110 * sensor.GetNode().Id ** 2, np.sum(element_sensor_view.GetTensorAdaptor().data))

    def test_SensorViewAuxiliaryTensorAdaptors(self):
        for sensor in self.sensors:
            nodal_sensor_view = KratosSI.Sensors.SensorView(sensor, "nodal_ta")
            ta = nodal_sensor_view.GetTensorAdaptor().Clone()
            ta.data[:] *= sensor.GetNode().Id
            nodal_sensor_view.AddAuxiliaryTensorAdaptor("mapped_id", ta)

            ta = nodal_sensor_view.GetTensorAdaptor().Clone()
            ta.data[:] *= sensor.GetNode().Id ** 2
            nodal_sensor_view.AddAuxiliaryTensorAdaptor("mapped_id_square", ta)

            self.assertIn("mapped_id", nodal_sensor_view.GetAuxiliarySuffixes())
            self.assertIn("mapped_id_square", nodal_sensor_view.GetAuxiliarySuffixes())
            self.assertEqual(np.linalg.norm(nodal_sensor_view.GetAuxiliaryTensorAdaptor("mapped_id").data - nodal_sensor_view.GetTensorAdaptor().data[:] * sensor.GetNode().Id), 0.0)
            self.assertEqual(np.linalg.norm(nodal_sensor_view.GetAuxiliaryTensorAdaptor("mapped_id_square").data - nodal_sensor_view.GetTensorAdaptor().data[:] * sensor.GetNode().Id ** 2), 0.0)

    def test_SetGetSensors(self):
        opt_prob = OptimizationProblem()
        data = ComponentDataView("test", opt_prob)

        SetSensors(data, self.sensors)

        list_of_sensors = GetSensors(data)
        for i, sensor in enumerate(list_of_sensors):
            self.assertEqual(sensor.GetName(), self.sensors[i].GetName())

if __name__ == '__main__':
    UnitTest.main()