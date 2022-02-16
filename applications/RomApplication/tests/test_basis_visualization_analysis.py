import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.RomApplication.basis_visualization_analysis import BasisVisualizationAnalysis

if KratosMultiphysics.kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication"):
    import KratosMultiphysics.StructuralMechanicsApplication


class TestBasisVisualizationAnalysis(KratosUnittest.TestCase):
    def setUp(self):
        self.files_to_remove = []

    def tearDown(self):
        for f in self.files_to_remove:
            KratosMultiphysics.kratos_utilities.DeleteFileIfExisting(f)

    def testSingleDOF(self):
        with KratosUnittest.WorkFolderScope("thermal_static_test_files", __file__):
            params = KratosMultiphysics.Parameters("""
            {
                "solver_settings": {
                    "model_part_name": "main_model_part",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "Square_Radiation_Stationary"
                    },
                    "rom_parameters_file" : "RomParameters.json"
                }
            }
            """)
            model = KratosMultiphysics.Model()

            analysis = BasisVisualizationAnalysis(model, params)
            analysis.Run()

            self.assertIn(KratosMultiphysics.TEMPERATURE, analysis._GetSolver().GetVariables())

            self.files_to_remove.append(os.path.abspath("Square_Radiation_Stationary_basis.post.bin"))
            self.files_to_remove.append(os.path.abspath("thermal_static_test_files.post.lst"))

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testMultipleDOF(self):
        with KratosUnittest.WorkFolderScope("structural_dynamic_test_files", __file__):
            params = KratosMultiphysics.Parameters("""
            {
                "solver_settings": {
                    "model_part_name": "main_model_part",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "Structure_Dynamic_2D"
                    },
                    "rom_parameters_file" : "RomParameters.json"
                }
            }
            """)
            model = KratosMultiphysics.Model()
            analysis = BasisVisualizationAnalysis(model, params)
            analysis.Run()

            self.assertIn(analysis._GetSolver().GetVariables(), KratosMultiphysics.StructuralMechanicsApplication.DISPLACEMENT_X)
            self.assertIn(analysis._GetSolver().GetVariables(), KratosMultiphysics.StructuralMechanicsApplication.DISPLACEMENT_Y)

            self.files_to_remove.append(os.path.abspath("Structure_Dynamic_2D.post.bin"))
            self.files_to_remove.append(os.path.abspath("Structure_Dynamic_2D.post.lst"))

if __name__ == '__main__':
    KratosUnittest.main()
