
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

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
            "filter_settings": {
                "type" : "implicit",
                "radius": 0.2,
                "linear_solver_settings" : {
                    "solver_type" : "amgcl",
                    "smoother_type":"ilu0",
                    "krylov_type": "gmres",
                    "coarsening_type": "aggregation",
                    "max_iteration": 200,
                    "provide_coordinates": false,
                    "gmres_krylov_space_dimension": 100,
                    "verbosity" : 0,
                    "tolerance": 1e-7,
                    "scaling": false,
                    "block_size": 1,
                    "use_block_matrices_if_possible" : true,
                    "coarse_enough" : 5000
                }
            },
            "initial_physical_thickness":0.015,
            "physical_thicknesses": [0.01,0.02]
        }""")

        cls.optimization_problem = OptimizationProblem()
        cls.thickness_control = ShellThicknessControl("test", cls.model, parameters, cls.optimization_problem)
        cls.optimization_problem.AddComponent(cls.thickness_control)

        with kratos_unittest.WorkFolderScope("../../control/thickness", __file__):
            Kratos.ModelPartIO("Thick_2x2_Shell", Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(cls.model_part)
            material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": "materials_2D.json"}} """)
            Kratos.ReadMaterialsUtility(material_settings, cls.model)

        cls.thickness_control.Initialize()

    @classmethod
    def tearDownClass(cls):
        with kratos_unittest.WorkFolderScope(".", __file__):
            DeleteFileIfExisting("Thick_2x2_Shell.time")

    def test_GetControlField(self):
        control_field = self.thickness_control.GetControlField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.5, 4)

    def test_GetPhysicalField(self):
        thickness_field = self.thickness_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(thickness_field), 0.03, 4)

    def test_MapGradient(self):
        physical_gradient = self.thickness_control.GetEmptyField()
        for element in physical_gradient.GetContainer():
            element.SetValue(KratosOA.THICKNESS_SENSITIVITY, element.GetGeometry().Area())
        Kratos.Expression.VariableExpressionIO.Read(physical_gradient, KratosOA.THICKNESS_SENSITIVITY)
        mapped_gradient = self.thickness_control.MapGradient({Kratos.THICKNESS: physical_gradient})
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(mapped_gradient), 0.5625, 4)

    def test_Update(self):
        update_field = self.thickness_control.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(update_field, 0.75)
        self.thickness_control.Update(update_field)
        control_field = self.thickness_control.GetControlField()
        thickness_field = self.thickness_control.GetPhysicalField()
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(control_field), 0.75, 4)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(thickness_field), 0.019999962733607157, 10)

if __name__ == "__main__":
    kratos_unittest.main()