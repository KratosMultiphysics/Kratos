
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.thickness.shell_thickness_control import ShellThicknessControl

class TestShellThicknessControl(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("shell")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": ["shell"],
            "physical_thicknesses"       : [0.01, 0.02],
            "filter_settings": {
                "filter_type"  : "implicit_filter",
                "filter_radius": 0.2
            },
            "thickness_projection_settings": {
                "type": "sigmoidal_projection",
                "beta_value": 25
            }
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.thickness_control = ShellThicknessControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.thickness_control)

        with kratos_unittest.WorkFolderScope(".", __file__):
            Kratos.ModelPartIO("Thick_2x2_Shell", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
        material_params = Kratos.Parameters("""{
            "properties": [
                {
                    "model_part_name": "shell",
                    "properties_id": 1,
                    "Material": {
                        "Variables": {
                            "THICKNESS": 0.015
                        },
                        "Tables": {}
                    }
                }
            ]
        }""")
        Kratos.ReadMaterialsUtility(cls.model).ReadMaterials(material_params)

        cls.thickness_control.Initialize()
        cls.initial_thickness = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(cls.model_part.Elements, Kratos.THICKNESS)
        cls.initial_thickness.CollectData()

    def setUp(self) -> None:
        self.initial_thickness.StoreData()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Thick_2x2_Shell.time")

    # def test_GetControlField(self):
    #     control_field = self.thickness_control.GetControlField()
    #     self.assertAlmostEqual(numpy.linalg.norm(control_field.data, ord=numpy.inf), 0.35, 4)

    # def test_GetPhysicalField(self):
    #     thickness_field = self.thickness_control.GetPhysicalField()
    #     self.assertAlmostEqual(numpy.linalg.norm(thickness_field.data), 0.03, 4)

    # def test_MapGradient(self):
    #     physical_gradient = self.thickness_control.GetEmptyField()
    #     for element in physical_gradient.GetContainer():
    #         element.SetValue(KratosOA.THICKNESS_SENSITIVITY, element.GetGeometry().DomainSize())
    #     Kratos.TensorAdaptors.VariableTensorAdaptor(physical_gradient, KratosOA.THICKNESS_SENSITIVITY, copy=False).CollectData()
    #     mapped_gradient = self.thickness_control.MapGradient({Kratos.THICKNESS: physical_gradient})
    #     self.assertAlmostEqual(numpy.linalg.norm(mapped_gradient.data), 0.5625, 4)

    # def test_Update(self):
    #     update_field = self.thickness_control.GetEmptyField()
    #     update_field.data[:] = 0.75
    #     self.thickness_control.Update(update_field)
    #     control_field = self.thickness_control.GetControlField()
    #     thickness_field = self.thickness_control.GetPhysicalField()
    #     self.assertAlmostEqual(numpy.linalg.norm(control_field.data, ord=numpy.inf), 0.75, 4)
    #     self.assertAlmostEqual(numpy.linalg.norm(thickness_field.data, ord=numpy.inf), 0.019999999998375886, 10)

    def test_AdaptiveBeta(self):
        parameters = Kratos.Parameters("""{
            "controlled_model_part_names": ["shell"],
            "physical_thicknesses"       : [0.01, 0.02],
            "filter_settings": {
                "filter_type"  : "implicit_filter",
                "filter_radius": 0.2
            },
            "thickness_projection_settings": {
                "type": "adaptive_sigmoidal_projection",
                "initial_value": 0.01,
                "max_value"    : 30,
                "increase_fac" : 1.05,
                "update_period": 3
            }
        }""")
        thickness_control = ShellThicknessControl("test_adaptive", self.model, parameters, self.optimization_problem)
        self.optimization_problem.AddComponent(thickness_control)
        thickness_control.Initialize()

        control_field = thickness_control.GetControlField()
        thickness_control.Update(control_field)
        for i in range(20):
            control_field.data *= 1.2
            self.assertTrue(thickness_control.Update(control_field), msg=f"Failed at iteration = {i} with update = {control_field.data}, control_field = {thickness_control.GetControlField().data}")
            self.optimization_problem.AdvanceStep()

        self.assertAlmostEqual(thickness_control.thickness_projection.beta, 0.014071004226562506)
        self.assertAlmostEqual(numpy.linalg.norm(thickness_control.GetPhysicalField().data), 0.03204326420888321)

if __name__ == "__main__":
    kratos_unittest.main()