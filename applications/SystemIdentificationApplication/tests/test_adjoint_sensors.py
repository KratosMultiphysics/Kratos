from math import sqrt
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


class TestTemperatureSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.sensor_model_part = cls.model.CreateModelPart("SensorModelPart")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)

        prop = cls.model_part.CreateNewProperties(1)

        cls.model_part.CreateNewElement("Element2D3N", 1, [1, 2, 4], prop)
        cls.model_part.CreateNewElement("Element2D3N", 2, [2, 3, 4], prop)

        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.TEMPERATURE, node.Id + 2)

        parameters = [
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_1",
                "value"        : 0,
                "location"     : [0.3333333333333, 0.3333333333333, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_2",
                "value"        : 0,
                "location"     : [0.6666666666667, 0.6666666666667, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_3",
                "value"        : 0,
                "location"     : [0.5, 0.5, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = CreateSensors(cls.sensor_model_part, cls.model_part, parameters)
        cls.ref_values = [13/3, 5, 5]

    def test_SensorsOnNodes(self):
        parameters = [
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_1",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_2",
                "value"        : 0,
                "location"     : [1.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_3",
                "value"        : 0,
                "location"     : [1.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_4",
                "value"        : 0,
                "location"     : [0.0, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        sensors = CreateSensors(self.model.CreateModelPart("SensorsOnNodes"), self.model_part, parameters)
        for sensor, ref_node_id in zip(sensors, [1, 2, 3, 4]):
            self.assertAlmostEqual(sensor.CalculateValue(self.model_part), self.model_part.GetNode(ref_node_id).GetSolutionStepValue(Kratos.TEMPERATURE))

    def test_SensorsOnEdges(self):
        parameters = [
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_1",
                "value"        : 0,
                "location"     : [0.5, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_2",
                "value"        : 0,
                "location"     : [1.0, 0.5, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_3",
                "value"        : 0,
                "location"     : [0.5, 0.5, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_4",
                "value"        : 0,
                "location"     : [0.5, 1.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "temperature_sensor",
                "name"         : "temp_x_5",
                "value"        : 0,
                "location"     : [0.0, 0.5, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        sensors = CreateSensors(self.model.CreateModelPart("SensorsOnEdges"), self.model_part, parameters)
        for sensor, (ref_node_id_1, ref_node_id_2) in zip(sensors, [(1, 2), (2, 3), (2, 4), (3, 4), (1, 4)]):
            ref_value = (self.model_part.GetNode(ref_node_id_1).GetSolutionStepValue(Kratos.TEMPERATURE) + self.model_part.GetNode(ref_node_id_2).GetSolutionStepValue(Kratos.TEMPERATURE)) / 2.0
            self.assertAlmostEqual(sensor.CalculateValue(self.model_part), ref_value)

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 7)

    def test_PartialSensitivity(self):
        sensitivity_matrix = Kratos.Matrix(3, 18)
        sensitivity_gradient = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            element: Kratos.Element = self.model_part.GetElement(sensor.GetNode().GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculatePartialSensitivity(element, Kratos.TEMPERATURE, sensitivity_matrix, sensitivity_gradient, self.model_part.ProcessInfo)
            for j, node in enumerate(element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.TEMPERATURE, node.GetSolutionStepValue(Kratos.TEMPERATURE) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, -1* sensitivity_gradient[j]) # -1 because structural responses are designed opposite to fluid ones and need a sign change
                node.SetSolutionStepValue(Kratos.TEMPERATURE, node.GetSolutionStepValue(Kratos.TEMPERATURE) - delta)


if __name__ == '__main__':
    UnitTest.main()