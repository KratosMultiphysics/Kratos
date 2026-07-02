import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import test_helper
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
from helper_utilities import _compare_case_outputs, run_orchestrator

class KratosGeoMechanicsMuskatTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests to test partially saturated groundwater flow in 2D.
    """

    def _assert_fully_saturated_flow(
        self, parent_name, test_name, node_number, expected_value
    ):
        file_path = test_helper.get_file_path(
            os.path.join(parent_name, test_name)
        )
        simulation = test_helper.run_kratos(file_path)
        water_pressure = [
            water_pressure for water_pressure in test_helper.get_water_pressure(simulation)
        ]
        self.assertAlmostEqual(expected_value, water_pressure[node_number - 1], 6)

    # To be changed, node 100 is a random pick and may be on the prescribed boundary. pick a better result, with computed result
    def no_test_fully_saturated_hydrostatic(self):
        self._assert_fully_saturated_flow('Muskat', 'fully_saturated_hydrostatic', 100, 22800.0)


    def no_test_fully_saturated_hydrostatic_cutoff(self):
        self._assert_fully_saturated_flow('Muskat', 'fully_saturated_hydrostatic_cutoff', 100, 0.0)
        

    def no_test_partially_saturated_van_Genuchten_with_2stage_orchestrator(self):
       # stage 1: fully saturated ( linear solution ) equal to test_fully_saturated_hydrostatic_cutoff
       # stage 2: try to reach nonlinear solution with van Genuchten

       parent_name = 'Muskat'
       test_name = 'partially_saturated_van_Genuchten_2stage'
       file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
       # comparison_data = [("test_stage1.post.res", 50.0, 1.0),
       #                    ("test_stage2.post.res", 100.0, 2.0)]

       project_parameters_filename = test_helper.get_file_path(os.path.join(file_path, "ProjectParameters.json"))
       # Parse simulation settings and run simulation
       with open(project_parameters_filename, 'r') as parameter_file:
           project_parameters = Kratos.Parameters(parameter_file.read())

       cwd = os.getcwd()
       os.chdir(file_path)
       run_orchestrator(project_parameters)

       # bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
       # for output_file_name, expected_total_reaction_y, stage_time in comparison_data:
       #     output_file_path = os.path.join(file_path, output_file_name)
       #     total_reaction_y = self.total_reaction_y_from_file(output_file_path, stage_time, bottom_node_ids)
       #     self.assertAlmostEqual(total_reaction_y, expected_total_reaction_y, places=3)

       os.chdir(cwd)


    def test_partially_saturated_van_Genuchten_given_startfield(self):
       parent_name = 'Muskat'
       test_name = 'partially_saturated_van_Genuchten_startfield'
       file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
       # comparison_data = [("test_stage1.post.res", 50.0, 1.0),
       #                    ("test_stage2.post.res", 100.0, 2.0)]

       project_parameters_filename = test_helper.get_file_path(os.path.join(file_path, "ProjectParameters.json"))
       # Parse simulation settings and run simulation
       with open(project_parameters_filename, 'r') as parameter_file:
           project_parameters = Kratos.Parameters(parameter_file.read())

       cwd = os.getcwd()
       os.chdir(file_path)
       run_orchestrator(project_parameters)

       # bottom_node_ids = [1, 2, 6, 11, 17, 25, 34, 46, 59, 75, 90]
       # for output_file_name, expected_total_reaction_y, stage_time in comparison_data:
       #     output_file_path = os.path.join(file_path, output_file_name)
       #     total_reaction_y = self.total_reaction_y_from_file(output_file_path, stage_time, bottom_node_ids)
       #     self.assertAlmostEqual(total_reaction_y, expected_total_reaction_y, places=3)

       os.chdir(cwd)


if __name__ == '__main__':
    KratosUnittest.main()
