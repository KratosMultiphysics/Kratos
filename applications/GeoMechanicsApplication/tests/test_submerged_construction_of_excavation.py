import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper

import os


class KratosGeoMechanicsSubmergedConstructionOfExcavation(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.unsaturated_unit_weight_of_clay = 16.0e3  # N/m3
        self.saturated_unit_weight_of_clay = 18.0e3  # N/m3
        self.saturated_unit_weight_of_sand = 20.0e3  # N/m3
        self.unit_weight_of_water = 10.0e3  # N/m3
        self.model_width = 65.0  # m
        self.excavation_width = 15.0  # m
        self.unsaturated_clay_layer_thickness = 2.0  # m
        self.saturated_clay_layer_thickness = 18.0  # m
        self.saturated_sand_layer_thickness = 30.0  # m
        self.excavated_middle_clay_layer_thickness = 8.0  # m
        self.excavated_lower_clay_layer_thickness = 10.0  # m
        self.weight_of_diaphragm_wall = 10.0e3  # N / m
        self.length_of_diaphragm_wall = 30.0  # m
        self.distributed_surface_load = 5.0e3  # N / m / m
        self.load_edge_length = 5.0  # m

    def total_reaction_y_from_output_data(self, output_data, time, node_ids):
        reactions = GiDOutputFileReader.nodal_values_at_time("REACTION", time, output_data, node_ids=node_ids)
        return sum([reaction[1] for reaction in reactions])

    def calculate_weight_of_all_soil(self):
        return self.model_width * (self.unsaturated_clay_layer_thickness * self.unsaturated_unit_weight_of_clay +
                                   self.saturated_clay_layer_thickness * self.saturated_unit_weight_of_clay +
                                   self.saturated_sand_layer_thickness * self.saturated_unit_weight_of_sand)

    def calculate_weight_of_excavated_clay_upper_right(self):
        return self.unsaturated_clay_layer_thickness * self.excavation_width * self.unsaturated_unit_weight_of_clay

    def calculate_weight_of_excavated_clay_middle_right(self):
        return self.excavated_middle_clay_layer_thickness * self.excavation_width * self.saturated_unit_weight_of_clay

    def calculate_weight_of_excavated_clay_lower_right(self):
        return self.excavated_lower_clay_layer_thickness * self.excavation_width * self.saturated_unit_weight_of_clay

    def calculate_weight_of_water_after_second_excavation(self):
        return self.excavated_middle_clay_layer_thickness * self.excavation_width * self.unit_weight_of_water

    def calculate_weight_of_water_after_third_excavation(self):
        return (self.excavated_middle_clay_layer_thickness + self.excavated_lower_clay_layer_thickness) * self.excavation_width * self.unit_weight_of_water

    def calculate_weight_of_diaphragm_wall(self):
        return self.weight_of_diaphragm_wall * self.length_of_diaphragm_wall

    def calculate_total_vertical_surface_load(self):
        return self.distributed_surface_load * self.load_edge_length

    def run_simulation(self, sub_directory_name):
        project_path = test_helper.get_file_path(os.path.join("submerged_construction_of_excavation", sub_directory_name))
        project_parameters_filenames = ["1_Initial_stage.json", "2_Null_step.json", "3_Wall_installation.json", "4_First_excavation.json", "5_Strut_installation.json", "6_Second_excavation.json", "7_Third_excavation.json"]

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for filename in project_parameters_filenames:
                with open(os.path.join("..", "common", filename), "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        # Check vertical reaction forces in the initial stage
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(os.path.join(project_path, "1_Initial_stage.post.res"))
        time = -1.0
        bottom_node_ids = [1, 3, 10, 25, 44, 68, 96, 132, 175, 218, 268, 318, 378, 434, 499, 573, 649, 732, 819, 909, 999, 1102, 1206, 1315, 1434, 1560, 1695, 1828, 1966, 2111, 2259, 2406, 2558, 2725, 2896, 3059, 3239, 3422, 3621, 3818, 4023, 4233, 4455, 4739, 5151, 5563, 5959, 6281, 6607, 6924, 7237, 7541, 7831]
        rel_tolerance = 0.005
        expected_total_weight = self.calculate_weight_of_all_soil()
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_weight, places=None, delta=rel_tolerance*expected_total_weight)

        # Check vertical reaction forces in the null step
        output_data = output_reader.read_output_from(os.path.join(project_path, "2_Null_step.post.res"))
        time = 0.0
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_weight, places=None, delta=rel_tolerance*expected_total_weight)

        # Check vertical reaction forces after installing the diaphragm wall
        output_data = output_reader.read_output_from(os.path.join(project_path, "3_Wall_installation.post.res"))
        time = 1.0
        expected_total_weight += self.calculate_weight_of_diaphragm_wall()
        expected_total_vertical_reaction = expected_total_weight + self.calculate_total_vertical_surface_load()
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_vertical_reaction, places=None, delta=rel_tolerance*expected_total_vertical_reaction)

        # Check vertical reaction forces after the first excavation
        output_data = output_reader.read_output_from(os.path.join(project_path, "4_First_excavation.post.res"))
        time = 2.0
        expected_total_weight -= self.calculate_weight_of_excavated_clay_upper_right()
        expected_total_vertical_reaction = expected_total_weight + self.calculate_total_vertical_surface_load()
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_vertical_reaction, places=None, delta=rel_tolerance*expected_total_vertical_reaction)

        # Check vertical reaction forces after strut installation (no changes with respect to the previous stage)
        output_data = output_reader.read_output_from(os.path.join(project_path, "5_Strut_installation.post.res"))
        time = 3.0
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_vertical_reaction, places=None, delta=rel_tolerance*expected_total_vertical_reaction)

        # Check vertical reaction forces after the second excavation
        output_data = output_reader.read_output_from(os.path.join(project_path, "6_Second_excavation.post.res"))
        time = 4.0
        expected_total_weight -= self.calculate_weight_of_excavated_clay_middle_right()
        expected_total_vertical_reaction = expected_total_weight + self.calculate_total_vertical_surface_load() + self.calculate_weight_of_water_after_second_excavation()
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_vertical_reaction, places=None, delta=rel_tolerance*expected_total_vertical_reaction)

        # Check vertical reaction forces after the third excavation
        output_data = output_reader.read_output_from(os.path.join(project_path, "7_Third_excavation.post.res"))
        time = 5.0
        expected_total_weight -= self.calculate_weight_of_excavated_clay_lower_right()
        expected_total_vertical_reaction = expected_total_weight + self.calculate_total_vertical_surface_load() + self.calculate_weight_of_water_after_third_excavation()
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_vertical_reaction, places=None, delta=rel_tolerance*expected_total_vertical_reaction)

    def test_simulation_with_linear_elastic_materials(self):
        self.run_simulation("linear_elastic")

    @KratosUnittest.skip("Test case is being developed")
    def test_simulation_with_Mohr_Coulomb(self):
        self.run_simulation("Mohr_Coulomb")

    @KratosUnittest.skip("Test case is being developed")
    def test_simulation_with_hardening_soil(self):
        self.run_simulation("hardening_soil")


if __name__ == "__main__":
    KratosUnittest.main()
