import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.RomApplication.basis_visualization_analysis import BasisVisualizationAnalysis


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
                    "model_part_name": "ThermalModelPart",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "Square_Radiation_Stationary"
                    },
                    "rom_parameters_file" : "RomParameters.json"
                },
                "problem_data": {
                    "parallel_type": "OpenMP",
                    "start_time": 0.0,
                    "time_step": 1.0,
                    "end_time": 1.0,
                    "echo_level": 0
                }
            }
            """)
            model = KratosMultiphysics.Model()

            analysis = BasisVisualizationAnalysis(model, params)
            analysis.Run()

            # self.files_to_remove.append(os.path.abspath("Square_Radiation_Stationary.post.bin"))

    @KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication")
    def testMultipleDOF(self):
        with KratosUnittest.WorkFolderScope("structural_dynamic_test_files", __file__):
            params = KratosMultiphysics.Parameters("""
            {
                "solver_settings": {
                    "model_part_name": "Parts_Structure",
                    "model_import_settings": {
                        "input_type": "mdpa",
                        "input_filename": "Structure_Dynamic_2D"
                    },
                    "rom_parameters_file" : "RomParameters.json"
                },
                "problem_data": {
                    "parallel_type": "OpenMP",
                    "start_time": 0.0,
                    "time_step": 1.0,
                    "end_time": 1.0,
                    "echo_level": 0
                }
            }
            """)
            model = KratosMultiphysics.Model()

            analysis = BasisVisualizationAnalysis(model, params)
            analysis.Run()

            # self.files_to_remove.append(os.path.abspath("Square_Radiation_Stationary.post.bin"))

if __name__ == '__main__':
    KratosUnittest.main()
