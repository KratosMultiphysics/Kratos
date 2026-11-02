from math import sqrt
from random import Random
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.StructuralMechanicsApplication as KratosStruct
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors
class TestDisplacementSensor(UnitTest.TestCase):
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

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

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
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "direction"    : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)
        cls.ref_values = [7/3, 3, (7/3 + 10/3)/sqrt(2), (3 + 4)/sqrt(2)]

    def test_SensorsOnNodes(self):
        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [1.0, 0.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [1.0, 1.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [0.0, 1.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        sensors = CreateSensors(self.model.CreateModelPart("SensorsOnNodes"), self.model_part, parameters)
        for sensor, ref_node_id in zip(sensors, [1, 2, 3, 4]):
            self.assertAlmostEqual(sensor.CalculateValue(self.model_part), self.model_part.GetNode(ref_node_id).GetSolutionStepValue(Kratos.DISPLACEMENT_X))

    def test_SensorsOnEdges(self):
        parameters = [
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_1",
                "value"        : 0,
                "location"     : [0.5, 0.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_2",
                "value"        : 0,
                "location"     : [1.0, 0.5, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_3",
                "value"        : 0,
                "location"     : [0.5, 0.5, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_4",
                "value"        : 0,
                "location"     : [0.5, 1.0, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "displacement_sensor",
                "name"         : "disp_x_5",
                "value"        : 0,
                "location"     : [0.0, 0.5, 0.0],
                "direction"    : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        sensors = CreateSensors(self.model.CreateModelPart("SensorsOnEdges"), self.model_part, parameters)
        for sensor, (ref_node_id_1, ref_node_id_2) in zip(sensors, [(1, 2), (2, 3), (2, 4), (3, 4), (1, 4)]):
            ref_value = (self.model_part.GetNode(ref_node_id_1).GetSolutionStepValue(Kratos.DISPLACEMENT_X) + self.model_part.GetNode(ref_node_id_2).GetSolutionStepValue(Kratos.DISPLACEMENT_X)) / 2.0
            self.assertAlmostEqual(sensor.CalculateValue(self.model_part), ref_value)

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 7)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            element: Kratos.Element = self.model_part.GetElement(sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for j, node in enumerate(element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[j * 6])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[j * 6 + 1])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


class TestStrainSensorShell(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.adjoint_model_part = cls.model.CreateModelPart("TestAdjoint")

        cls.model_part.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = 1e-10

        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ROTATION)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_ROTATION)

        cls.adjoint_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.adjoint_model_part.AddNodalSolutionStepVariable(Kratos.ROTATION)
        cls.adjoint_model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
        cls.adjoint_model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_ROTATION)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("ShellThinElement3D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 2, [2, 3, 4], prop)

        # set the constitutive laws
        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "auxiliary_files/shell_material_properties.json"}}""")
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(cls.model_part, cls.model_part.Elements, False)

        Kratos.ConnectivityPreserveModeler().GenerateModelPart(cls.model_part, cls.adjoint_model_part, "AdjointFiniteDifferencingShellThinElement3D3N")

        for element in cls.model_part.Elements:
            element.Initialize(cls.model_part.ProcessInfo)

        for element in cls.adjoint_model_part.Elements:
            element.Initialize(cls.adjoint_model_part.ProcessInfo)

        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.ROTATION_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Z, cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
            node.SetSolutionStepValue(Kratos.ROTATION, [node.Id * 2, node.Id * 2 + 1, node.Id * 2 + 2])

        parameters = [
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_x_1",
                "value"          : 0,
                "location"       : [0.3333333333333, 0.3333333333333, 0.0],
                "strain_type"    : "strain_xx",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_x_2",
                "value"          : 0,
                "location"       : [0.6666666666667, 0.6666666666667, 0.0],
                "strain_type"    : "strain_xx",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_y_1",
                "value"          : 0,
                "location"       : [0.3333333333333, 0.3333333333333, 0.0],
                "strain_type"    : "strain_yy",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{

                "type"           : "strain_sensor",
                "name"           : "strain_y_2",
                "value"          : 0,
                "location"       : [0.6666666666667, 0.6666666666667, 0.0],
                "strain_type"    : "strain_yy",
                "strain_variable": "SHELL_STRAIN",
                "weight"         : 1.0,
                "variable_data"  : {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)
        cls.ref_values = [0.5, -1.5, 4.5, 0.5]

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 3)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            adjoint_element: Kratos.Element = self.adjoint_model_part.GetElement(sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(adjoint_element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(adjoint_element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)

    def test_CalculateGradientRestarted(self):
        # a restarted run skips element Initialize(), so the FD clone must be initialized with IS_RESTARTED = False
        residual_matrix = Kratos.Matrix(18, 18)
        process_info = self.model_part.ProcessInfo

        def calculate_gradients():
            gradients = []
            for sensor in self.sensors:
                element_id = sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID)
                adjoint_element = self.adjoint_model_part.GetElement(element_id)
                gradient = Kratos.Vector()
                sensor.CalculateGradient(adjoint_element, residual_matrix, gradient, process_info)
                gradients.append(gradient)
            return gradients

        ref_gradients = calculate_gradients()
        process_info[Kratos.IS_RESTARTED] = True
        try:
            restarted_gradients = calculate_gradients()
        finally:
            process_info[Kratos.IS_RESTARTED] = False

        for ref_gradient, restarted_gradient in zip(ref_gradients, restarted_gradients):
            self.assertVectorAlmostEqual(ref_gradient, restarted_gradient, 12)


class TestStrainSensorSolids(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.adjoint_model_part = cls.model.CreateModelPart("TestAdjoint")

        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.model_part.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = 1e-10
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)

        cls.adjoint_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.adjoint_model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)

        cls.model_part.CreateNewNode( 1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode( 2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode( 3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode( 4, 0.0, 1.0, 0.0)
        cls.model_part.CreateNewNode( 5, 0.0, 0.0, 1.0)
        cls.model_part.CreateNewNode( 6, 1.0, 0.0, 1.0)
        cls.model_part.CreateNewNode( 7, 1.0, 1.0, 1.0)
        cls.model_part.CreateNewNode( 8, 0.0, 1.0, 1.0)
        cls.model_part.CreateNewNode( 9, 2.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(10, 2.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(11, 2.0, 0.0, 1.0)
        cls.model_part.CreateNewNode(12, 2.0, 1.0, 1.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("SmallDisplacementElement3D8N", 1, [1, 2,  3, 4, 5,  6,  7, 8], prop)
        cls.model_part.CreateNewElement("SmallDisplacementElement3D8N", 2, [2, 9, 10, 3, 6, 11, 12, 7], prop)

        # set the constitutive laws
        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "auxiliary_files/solid_material_properties.json"}}""")
        Kratos.ReadMaterialsUtility(material_settings, cls.model)

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(cls.model_part, cls.model_part.Elements, False)

        Kratos.ConnectivityPreserveModeler().GenerateModelPart(cls.model_part, cls.adjoint_model_part, "AdjointFiniteDifferencingSmallDisplacementElement3D8N")

        for element in cls.model_part.Elements:
            element.Initialize(cls.model_part.ProcessInfo)

        for element in cls.adjoint_model_part.Elements:
            element.Initialize(cls.adjoint_model_part.ProcessInfo)

        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Z, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_X, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Y, cls.model_part)
        Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Z, cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, [node.Id, node.Id + 1, node.Id + 2])
            node.SetValue(Kratos.DISPLACEMENT, node.GetSolutionStepValue(Kratos.DISPLACEMENT))

        parameters = [
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_x_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_xx",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_x_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_xx",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_y_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_yy",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_y_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_yy",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_xy_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_xy",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_xy_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_xy",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_xz_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_xz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_xz_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_xz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_yz_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_yz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_yz_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_yz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_zz_1",
                "value"          : 0,
                "location"       : [0.5, 0.5, 0.5],
                "strain_type"    : "strain_zz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }"""),
            Kratos.Parameters("""{
                "type"           : "strain_sensor",
                "name"           : "strain_yz_2",
                "value"          : 0,
                "location"       : [1.5, 0.5, 0.5],
                "strain_type"    : "strain_zz",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR",
                "weight"         : 1.0,
                "variable_data"  : {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)
        cls.ref_values = [0, 6.0, 2.0, 1.0, 1.0, 3.5, 2.0, 4.5, 3.0, 2.0, 4.0, 3.0]

    def tearDown(self):
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.DISPLACEMENT, node.GetValue(Kratos.DISPLACEMENT))

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 3)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(24, 24)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            adjoint_element: Kratos.Element = self.adjoint_model_part.GetElement(sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(adjoint_element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for j, node in enumerate(adjoint_element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[j * 3], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[j * 3 + 1], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Z, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Z) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[j * 3 + 2], 4)
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Z, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Z) - delta)

    def test_CalculateGradientRestarted(self):
        # a restarted run skips element Initialize(), so the FD clone must be initialized with IS_RESTARTED = False
        residual_matrix = Kratos.Matrix(24, 24)
        process_info = self.model_part.ProcessInfo

        def calculate_gradients():
            gradients = []
            for sensor in self.sensors:
                element_id = sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID)
                adjoint_element = self.adjoint_model_part.GetElement(element_id)
                gradient = Kratos.Vector()
                sensor.CalculateGradient(adjoint_element, residual_matrix, gradient, process_info)
                gradients.append(gradient)
            return gradients

        ref_gradients = calculate_gradients()
        process_info[Kratos.IS_RESTARTED] = True
        try:
            restarted_gradients = calculate_gradients()
        finally:
            process_info[Kratos.IS_RESTARTED] = False

        for ref_gradient, restarted_gradient in zip(ref_gradients, restarted_gradients):
            self.assertVectorAlmostEqual(ref_gradient, restarted_gradient, 12)


class TestStrainSensorThreadSafety(UnitTest.TestCase):
    """Strain-sensor gradients assembled by the parallel block builder with the adjoint scheme must not depend on the
    number of threads. StrainSensor::CalculateGradient used to perturb the shared nodal DISPLACEMENT in place, so
    elements assembled at the same time by other threads read perturbed values (and concurrent +=/-= could leave the
    nodes shifted). Every element carries a strain sensor, and the element ids are shuffled so that neighbouring
    elements fall into different chunks of the builder's parallel loop (schedule(guided, 512)) and run concurrently."""

    @staticmethod
    def _AssembleAdjointRhs(num_threads: int, num_assemblies: int):
        # before the response is created: MeasurementResidualResponseFunction sizes its per-thread buffers
        Kratos.ParallelUtilities.SetNumThreads(num_threads)
        model = Kratos.Model()
        model_part = model.CreateModelPart("Test")
        adjoint_model_part = model.CreateModelPart("TestAdjoint")
        for mp in (model_part, adjoint_model_part):
            mp.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
            mp.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
            mp.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
            mp.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = 1e-8

        domain = Kratos.Hexahedra3D8(
            Kratos.Node(1, 0.0, 0.0, 0.0), Kratos.Node(2, 4.0, 0.0, 0.0),
            Kratos.Node(3, 4.0, 1.0, 0.0), Kratos.Node(4, 0.0, 1.0, 0.0),
            Kratos.Node(5, 0.0, 0.0, 1.0), Kratos.Node(6, 4.0, 0.0, 1.0),
            Kratos.Node(7, 4.0, 1.0, 1.0), Kratos.Node(8, 0.0, 1.0, 1.0))
        mesher_parameters = Kratos.Parameters("""{
            "element_name": "SmallDisplacementElement3D4N", "condition_name": "SurfaceCondition3D3N",
            "create_skin_sub_model_part": false, "number_of_divisions": 10}""")
        mesh = model.CreateModelPart("Mesh")
        mesh.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        Kratos.StructuredMeshGeneratorProcess(domain, mesh, mesher_parameters).Execute()
        for node in mesh.Nodes:
            model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
        properties = model_part.CreateNewProperties(1)
        new_ids = list(range(1, mesh.NumberOfElements() + 1))
        Random(0).shuffle(new_ids)
        for new_id, element in zip(new_ids, mesh.Elements):
            model_part.CreateNewElement(
                "SmallDisplacementElement3D4N", new_id, [n.Id for n in element.GetGeometry()], properties)
        properties[Kratos.YOUNG_MODULUS] = 30e9
        properties[Kratos.POISSON_RATIO] = 0.2
        properties[Kratos.CONSTITUTIVE_LAW] = KratosStruct.LinearElastic3DLaw()
        for node in model_part.Nodes:      # smooth, non-rigid field
            node.SetSolutionStepValue(
                Kratos.DISPLACEMENT, [1e-3 * node.X * node.Y, 2e-3 * node.Z * node.X ** 2, -1e-3 * node.Y * node.Z])

        Kratos.ConnectivityPreserveModeler().GenerateModelPart(
            model_part, adjoint_model_part, "AdjointFiniteDifferencingSmallDisplacementElement3D4N")
        for mp in (model_part, adjoint_model_part):
            for element in mp.Elements:
                element.Initialize(mp.ProcessInfo)
        for variable in (KratosStruct.ADJOINT_DISPLACEMENT_X, KratosStruct.ADJOINT_DISPLACEMENT_Y,
                         KratosStruct.ADJOINT_DISPLACEMENT_Z):
            Kratos.VariableUtils().AddDof(variable, adjoint_model_part)

        sensor_parameters = []
        for element in model_part.Elements:
            centre = element.GetGeometry().Center()
            sensor_parameters.append(Kratos.Parameters(f"""{{
                "type": "strain_sensor", "name": "strain_{element.Id}", "value": 0,
                "location": [{centre.X}, {centre.Y}, {centre.Z}], "strain_type": "strain_xx",
                "strain_variable": "GREEN_LAGRANGE_STRAIN_TENSOR", "weight": 1.0, "variable_data": {{}}}}"""))
        sensors = CreateSensors(model.CreateModelPart("Sensors"), model_part, sensor_parameters)
        response = KratosSI.Responses.MeasurementResidualResponseFunction(1.0)
        for sensor in sensors:
            sensor.GetNode().SetValue(KratosSI.SENSOR_MEASURED_VALUE, 1e-3)
            sensor.GetNode().SetValue(KratosSI.SENSOR_NORMALIZATION_FACTOR, 1.0)
            response.AddSensor(sensor)
        response.Initialize()
        response.CalculateValue(model_part)

        displacement_before = [list(node.GetSolutionStepValue(Kratos.DISPLACEMENT)) for node in model_part.Nodes]
        scheme = Kratos.ResidualBasedAdjointStaticScheme(response)
        builder = Kratos.ResidualBasedBlockBuilderAndSolver(Kratos.SkylineLUFactorizationSolver())
        space = Kratos.SparseSpace()
        A, Dx, b = space.CreateEmptyMatrixPointer(), space.CreateEmptyVectorPointer(), space.CreateEmptyVectorPointer()
        scheme.Initialize(adjoint_model_part)
        builder.SetUpDofSet(scheme, adjoint_model_part)
        builder.SetUpSystem(adjoint_model_part)
        builder.ResizeAndInitializeVectors(scheme, A, Dx, b, adjoint_model_part)
        assemblies = []
        for _ in range(num_assemblies):
            rhs = Kratos.Vector(builder.GetEquationSystemSize())
            space.SetToZeroVector(rhs)
            builder.BuildRHS(scheme, adjoint_model_part, rhs)
            assemblies.append(
                {(dof.Id(), dof.GetVariable().Name()): rhs[dof.EquationId] for dof in builder.GetDofSet()})
        displacement_after = [list(node.GetSolutionStepValue(Kratos.DISPLACEMENT)) for node in model_part.Nodes]
        # FD swaps node copies into the element geometries: they must hold the model part's own nodes again
        for node in adjoint_model_part.Nodes:
            node.SetValue(Kratos.NODAL_MASS, float(node.Id))
        original_nodes_restored = all(
            node.GetValue(Kratos.NODAL_MASS) == float(node.Id)
            for element in adjoint_model_part.Elements for node in element.GetGeometry())
        return assemblies, displacement_before, displacement_after, original_nodes_restored

    def test_CalculateGradientThreadCount(self):
        initial_threads = Kratos.ParallelUtilities.GetNumThreads()
        try:
            serial, _, _, restored_serial = self._AssembleAdjointRhs(1, 1)
            parallel, before, after, restored_parallel = self._AssembleAdjointRhs(max(8, initial_threads), 3)
        finally:
            Kratos.ParallelUtilities.SetNumThreads(initial_threads)
        self.assertTrue(restored_serial and restored_parallel)
        self.assertEqual(before, after)      # the shared nodal displacements are untouched
        reference = serial[0]
        scale = max(abs(v) for v in reference.values())
        self.assertGreater(scale, 0.0)
        for assembly in parallel:
            self.assertEqual(reference.keys(), assembly.keys())
            worst = max(abs(reference[k] - assembly[k]) for k in reference)
            self.assertLess(worst, 1e-10 * scale)


if __name__ == '__main__':
    UnitTest.main()
