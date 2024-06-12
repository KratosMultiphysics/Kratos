from math import sqrt, pi
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.StructuralMechanicsApplication as KratosStruct
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
class TestDisplacementSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
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

        cls.sensors = GetSensors(cls.model_part, parameters)
        cls.ref_values = [7/3, 3, (7/3 + 10/3)/sqrt(2), (3 + 4)/sqrt(2)]

    def test_CalculateValue(self):
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 7)

    def test_CalculateGradient(self):
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = 1e-5

            element: Kratos.Element = self.model_part.GetElement(sensor.GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_X, node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) - delta)

                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) + delta)
                perturbed_value = sensor.CalculateValue(self.model_part)
                sensitivity = (perturbed_value - ref_value) / delta
                self.assertAlmostEqual(sensitivity, response_sensitivities[i * 6 + 1])
                node.SetSolutionStepValue(Kratos.DISPLACEMENT_Y, node.GetSolutionStepValue(Kratos.DISPLACEMENT_Y) - delta)


class TestStrainSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.model_part.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = 1e-10
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ROTATION)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_ROTATION)

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

        for element in cls.model_part.Elements:
            element.Initialize(cls.model_part.ProcessInfo)

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

        cls.sensors = GetSensors(cls.model_part, parameters)
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

            element: Kratos.Element = self.model_part.GetElement(sensor.GetValue(KratosSI.SENSOR_ELEMENT_ID))
            sensor.CalculateGradient(element, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
            for i, node in enumerate(element.GetGeometry()):
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


class TestEigenvalueSensor(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("Test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ROTATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION_MOMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_DISPLACEMENT)
        cls.model_part.AddNodalSolutionStepVariable(KratosStruct.ADJOINT_ROTATION)
        cls.delta = 1e4
        cls.model_part.ProcessInfo[KratosStruct.PERTURBATION_SIZE] = cls.delta
        cls.model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 3)

        cls.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        cls.model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        cls.model_part.CreateNewNode(5, 0.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(6, 1.0, 2.0, 0.0)
        cls.model_part.CreateNewNode(7, 1.0, 3.0, 0.0)
        cls.model_part.CreateNewNode(8, 0.0, 3.0, 0.0)
        cls.model_part.CreateNewNode(9, 0.0, 4.0, 0.0)
        cls.model_part.CreateNewNode(10, 1.0, 4.0, 0.0)

        for node in cls.model_part.Nodes:
            Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_X, Kratos.REACTION_X, cls.model_part)
            Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Y, Kratos.REACTION_Y, cls.model_part)
            Kratos.VariableUtils().AddDof(Kratos.DISPLACEMENT_Z, Kratos.REACTION_Z, cls.model_part)
            Kratos.VariableUtils().AddDof(Kratos.ROTATION_X, Kratos.REACTION_MOMENT_X, cls.model_part)
            Kratos.VariableUtils().AddDof(Kratos.ROTATION_Y, Kratos.REACTION_MOMENT_Y, cls.model_part)
            Kratos.VariableUtils().AddDof(Kratos.ROTATION_Z, Kratos.REACTION_MOMENT_Z, cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_X,  cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Y, cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_DISPLACEMENT_Z, cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_X,  cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Y,  cls.model_part)
            # Kratos.VariableUtils().AddDof(KratosStruct.ADJOINT_ROTATION_Z, cls.model_part)

        prop1 = cls.model_part.CreateNewProperties(1)
        # materials
        cl1 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[1].SetValue(Kratos.CONSTITUTIVE_LAW, cl1)
        cls.model_part.GetProperties()[1].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[1].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[1].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[1].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[1].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop2 = cls.model_part.CreateNewProperties(2)
        # materials
        cl2 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[2].SetValue(Kratos.CONSTITUTIVE_LAW, cl2)
        cls.model_part.GetProperties()[2].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[2].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[2].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[2].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[2].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop3 = cls.model_part.CreateNewProperties(3)
        # materials
        cl3 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[3].SetValue(Kratos.CONSTITUTIVE_LAW, cl3)
        cls.model_part.GetProperties()[3].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[3].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[3].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[3].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[3].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop4 = cls.model_part.CreateNewProperties(4)
        # materials
        cl4 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[4].SetValue(Kratos.CONSTITUTIVE_LAW, cl4)
        cls.model_part.GetProperties()[4].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[4].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[4].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[4].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[4].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop5 = cls.model_part.CreateNewProperties(5)
        # materials
        cl5 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[5].SetValue(Kratos.CONSTITUTIVE_LAW, cl5)
        cls.model_part.GetProperties()[5].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[5].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[5].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[5].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[5].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop6 = cls.model_part.CreateNewProperties(6)
        # materials
        cl6 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[6].SetValue(Kratos.CONSTITUTIVE_LAW, cl6)
        cls.model_part.GetProperties()[6].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[6].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[6].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[6].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[6].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop7 = cls.model_part.CreateNewProperties(7)
        # materials
        cl7 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[7].SetValue(Kratos.CONSTITUTIVE_LAW, cl7)
        cls.model_part.GetProperties()[7].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[7].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[7].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[7].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[7].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        prop8 = cls.model_part.CreateNewProperties(8)
        # materials
        cl8 = Kratos.KratosGlobals.GetConstitutiveLaw("LinearElasticPlaneStress2DLaw").Clone()
        cls.model_part.GetProperties()[8].SetValue(Kratos.CONSTITUTIVE_LAW, cl8)
        cls.model_part.GetProperties()[8].SetValue(Kratos.YOUNG_MODULUS, 3.0e+10)
        cls.model_part.GetProperties()[8].SetValue(Kratos.POISSON_RATIO, 0.29)
        cls.model_part.GetProperties()[8].SetValue(Kratos.DENSITY, 7800)
        cls.model_part.GetProperties()[8].SetValue(Kratos.THICKNESS, 0.1)
        cls.model_part.GetProperties()[8].SetValue(Kratos.COMPUTE_LUMPED_MASS_MATRIX, True)

        # create elements
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 1, [1, 2, 4], prop1)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 2, [2, 3, 4], prop2)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 3, [4, 3, 6], prop3)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 4, [4, 6, 5], prop4)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 5, [5, 6, 8], prop5)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 6, [6, 7, 8], prop6)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 7, [8, 7, 10], prop7)
        cls.model_part.CreateNewElement("ShellThinElement3D3N", 8, [8, 10, 9], prop8)

        # fixing BCs
        n1: Kratos.Node = cls.model_part.GetNode(1)
        n1.Fix(Kratos.DISPLACEMENT_X)
        n1.Fix(Kratos.DISPLACEMENT_Y)
        n1.Fix(Kratos.DISPLACEMENT_Z)
        n1.Fix(Kratos.ROTATION_X)
        n1.Fix(Kratos.ROTATION_Y)
        n1.Fix(Kratos.ROTATION_Z)  
        # n1.Fix(KratosStruct.ADJOINT_DISPLACEMENT_X)
        # n1.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Y)
        # n1.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n1.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n1.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n1.Fix(KratosStruct.ADJOINT_ROTATION_Z)        
        
        n2: Kratos.Node = cls.model_part.GetNode(2)
        n2.Fix(Kratos.DISPLACEMENT_Z)
        n2.Fix(Kratos.ROTATION_X)
        n2.Fix(Kratos.ROTATION_Y)
        n2.Fix(Kratos.ROTATION_Z) 
        # n2.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n2.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n2.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n2.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n3: Kratos.Node = cls.model_part.GetNode(3)
        n3.Fix(Kratos.DISPLACEMENT_Z)
        n3.Fix(Kratos.ROTATION_X)
        n3.Fix(Kratos.ROTATION_Y)
        n3.Fix(Kratos.ROTATION_Z) 
        # n3.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n3.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n3.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n3.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n4: Kratos.Node = cls.model_part.GetNode(4)
        n4.Fix(Kratos.DISPLACEMENT_Z)
        n4.Fix(Kratos.ROTATION_X)
        n4.Fix(Kratos.ROTATION_Y)
        n4.Fix(Kratos.ROTATION_Z) 
        # n4.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n4.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n4.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n4.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n5: Kratos.Node = cls.model_part.GetNode(5)
        n5.Fix(Kratos.DISPLACEMENT_Z)
        n5.Fix(Kratos.ROTATION_X)
        n5.Fix(Kratos.ROTATION_Y)
        n5.Fix(Kratos.ROTATION_Z) 
        # n5.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n5.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n5.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n5.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n6: Kratos.Node = cls.model_part.GetNode(6)
        n6.Fix(Kratos.DISPLACEMENT_Z)
        n6.Fix(Kratos.ROTATION_X)
        n6.Fix(Kratos.ROTATION_Y)
        n6.Fix(Kratos.ROTATION_Z) 
        # n6.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n6.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n6.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n6.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n7: Kratos.Node = cls.model_part.GetNode(7)
        n7.Fix(Kratos.DISPLACEMENT_Z)
        n7.Fix(Kratos.ROTATION_X)
        n7.Fix(Kratos.ROTATION_Y)
        n7.Fix(Kratos.ROTATION_Z) 
        # n7.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n7.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n7.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n7.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n8: Kratos.Node = cls.model_part.GetNode(8)
        n8.Fix(Kratos.DISPLACEMENT_Z)
        n8.Fix(Kratos.ROTATION_X)
        n8.Fix(Kratos.ROTATION_Y)
        n8.Fix(Kratos.ROTATION_Z) 
        # n8.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n8.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n8.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n8.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n9: Kratos.Node = cls.model_part.GetNode(9)
        n9.Fix(Kratos.DISPLACEMENT_Z)
        n9.Fix(Kratos.ROTATION_X)
        n9.Fix(Kratos.ROTATION_Y)
        n9.Fix(Kratos.ROTATION_Z) 
        # n9.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n9.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n9.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n9.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        n10: Kratos.Node = cls.model_part.GetNode(10)
        n10.Fix(Kratos.DISPLACEMENT_Z)
        n10.Fix(Kratos.ROTATION_X)
        n10.Fix(Kratos.ROTATION_Y)
        n10.Fix(Kratos.ROTATION_Z) 
        # n10.Fix(KratosStruct.ADJOINT_DISPLACEMENT_Z)
        # n10.Fix(KratosStruct.ADJOINT_ROTATION_X)
        # n10.Fix(KratosStruct.ADJOINT_ROTATION_Y)
        # n10.Fix(KratosStruct.ADJOINT_ROTATION_Z) 

        # cls.adjoint_model_part = cls.model.CreateModelPart("AdjointStructure")
        # KratosOA.OptimizationUtils.CopySolutionStepVariablesList(cls.adjoint_model_part, cls.model_part)
        # connectivity_preserve_modeller = Kratos.ConnectivityPreserveModeler()
        # connectivity_preserve_modeller.GenerateModelPart(cls.model_part, cls.adjoint_model_part, "AdjointFiniteDifferencingShellThinElement3D3N")

        #solve
        cls.eigensolver_settings = Kratos.Parameters("""
            {
                "eigen_system_settings" : {
                        "solver_type"       : "spectra_sym_g_eigs_shift",
                        "number_of_eigenvalues": 6,
                        "max_iteration": 1000,
                        "echo_level": 0
                    }
            }
            """)
        from KratosMultiphysics import eigen_solver_factory
        cls.solution_scheme = KratosStruct.EigensolverDynamicScheme()
        cls.eigen_linear_solver = eigen_solver_factory.ConstructSolver(cls.eigensolver_settings["eigen_system_settings"])
        cls.builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(cls.eigen_linear_solver)
        # cls.strategy = KratosStruct.EigensolverStrategy(cls.model_part, solution_scheme, builder_and_solver, 0.0, 1.0)
        # cls.strategy.Solve()

        parameters = [
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_1",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_2",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_3",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_4",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_5",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }"""),
            Kratos.Parameters("""{

                "type"         : "eigenvalue_sensor",
                "name"         : "eigenvalue_6",
                "value"        : 0,
                "location"     : [0.0, 0.0, 0.0],
                "weight"       : 1.0,
                "variable_data": {}
            }""")
        ]

        cls.sensors = GetSensors(cls.model_part, parameters)
        #cls.ref_values = [sqrt(5841715.635361228) / (2 * pi), sqrt(8784068.422526835) / (2 * pi)] # 1 element 3 nodes in Hz
        #cls.ref_values = [5841715.635361228, 8784068.422526835] # 1 element 3 nodes in lambda
        
        #cls.ref_values = [sqrt(2679590.868398367) / (2 * pi), sqrt(4345883.404962238) / (2 * pi), 
        #                  sqrt(5678340.053215294)/(2*pi), sqrt(5678351.763857347)/(2*pi)] # 2 elements 4 nodes in Hz
        
        #cls.ref_values = [sqrt(1017115.9373842935) / (2 * pi), sqrt(1626959.4806632798) / (2 * pi)] #4 elements 6 nodes in Hz
        #cls.ref_values = [1017115.9373842935, 1626959.4806632798] #4 elements 6 nodes in lambda
        
        #cls.ref_values = [516066.37176830304, 811665.5479598092, 4082724.576648915, 5468681.136183668] #6 elements 8 nodes in lambda

        #cls.ref_values = [309629.86927694024, 482679.56889471033, 2608512.25072338, 
        #                  3736893.4374975483, 5775693.198165001, 7302176.233452227] #8 elements 10 nodes in lambda
        cls.ref_values = [sqrt(309629.86927694024)/(2*pi), sqrt(482679.56889471033)/(2*pi), 
                          sqrt(2608512.25072338)/(2*pi), sqrt(3736893.4374975483)/(2*pi), 
                          sqrt(5775693.198165001)/(2*pi), sqrt(7302176.233452227)/(2*pi)] #8 elements 10 nodes in Hz
        
    def test_CalculateValue(self):
        self.strategy = KratosStruct.EigensolverStrategy(self.model_part, self.solution_scheme, self.builder_and_solver, 0.0, 1.0)
        self.strategy.Solve()
        values = [sensor.CalculateValue(self.model_part) for sensor in self.sensors]
        self.assertVectorAlmostEqual(values, self.ref_values, 8)

    def test_CalculateGradient(self):
        # def append_line_to_file(file_name, line):
        #     with open(file_name, 'a') as file:
        #         file.write(line + '\n')
        # file_name = "output.txt"
        residual_matrix = Kratos.Matrix(18, 18)
        response_sensitivities = Kratos.Vector()
        for i, sensor in enumerate(self.sensors):
            ref_value = self.ref_values[i]
            delta = self.delta
            
            for j, element in enumerate(self.model_part.Elements):
                self.strategy = KratosStruct.EigensolverStrategy(self.model_part, self.solution_scheme, self.builder_and_solver, 0.0, 1.0)
                self.strategy.Solve()
                sensor.CalculatePartialSensitivity(element, Kratos.YOUNG_MODULUS, residual_matrix, response_sensitivities, self.model_part.ProcessInfo)
                print("sensitivity of eigenvalue ", i+1, " wrt element ",j+1, " is ", response_sensitivities[0])
                # young_modulus = self.model_part.GetProperties()[j+1].GetValue(Kratos.YOUNG_MODULUS)

                # self.model_part.GetProperties()[j+1].SetValue(Kratos.YOUNG_MODULUS, young_modulus+delta)
                
                # self.strategy = KratosStruct.EigensolverStrategy(self.model_part, self.solution_scheme, self.builder_and_solver, 0.0, 1.0)
                # self.strategy.Solve()
                # perturbed_value = sensor.CalculateValue(self.model_part)
                
                # sensitivity = (perturbed_value - ref_value) / delta
                # print("sensitivity by finite difference is ", sensitivity)
                # self.model_part.GetProperties()[j+1].SetValue(Kratos.YOUNG_MODULUS, young_modulus)
                #line = "sensitivity of eigenvalue "+ str(i+1)+ " wrt element "+str(j+1)+ " is "+ response_sensitivities[0]
                #append_line_to_file(file_name, line)
                #print("sensitivity of eigenvalue ", i+1, " wrt element ",j+1, " is ", response_sensitivities[0])
                # self.assertAlmostEqual(sensitivity, response_sensitivities[0],10)

                

if __name__ == '__main__':
    UnitTest.main()