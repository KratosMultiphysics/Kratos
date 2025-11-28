import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
import test_helper

import os


class KratosGeoMechanicsSubmergedConstructionOfExcavation(KratosUnittest.TestCase):
    STAGE = {
        "initial_stage": {"time": -1.0, "name": "1_Initial_stage"},
        "null_step": {"time": 0.0, "name": "2_Null_step"},
        "wall_installation": {"time": 1.0, "name": "3_Wall_installation"},
        "first_excavation": {"time": 2.0, "name": "4_First_excavation"},
        "strut_installation": {"time": 3.0, "name": "5_Strut_installation"},
        "second_excavation": {"time": 4.0, "name": "6_Second_excavation"},
        "third_excavation": {"time": 5.0, "name": "7_Third_excavation"},
    }

    BOTTOM_NODE_IDS = [1, 3, 10, 25, 44, 68, 96, 132, 175, 218, 268, 318, 378, 434, 499, 573, 649, 732, 819, 909, 999,
                       1102, 1206, 1315, 1434, 1560, 1695, 1828, 1966, 2111, 2259, 2406, 2558, 2725, 2896, 3059, 3239,
                       3422, 3621, 3818, 4023, 4233, 4455, 4739, 5151, 5563, 5959, 6281, 6607, 6924, 7237, 7541, 7831]

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
        reactions = GiDOutputFileReader.nodal_values_at_time(
            "REACTION", time, output_data, node_ids=node_ids
        )
        return sum([reaction[1] for reaction in reactions])

    def calculate_weight_of_all_soil(self):
        return self.model_width * (
            self.unsaturated_clay_layer_thickness * self.unsaturated_unit_weight_of_clay
            + self.saturated_clay_layer_thickness * self.saturated_unit_weight_of_clay
            + self.saturated_sand_layer_thickness * self.saturated_unit_weight_of_sand
        )

    def calculate_weight_of_excavated_clay_upper_right(self):
        return (
            self.unsaturated_clay_layer_thickness
            * self.excavation_width
            * self.unsaturated_unit_weight_of_clay
        )

    def calculate_weight_of_excavated_clay_middle_right(self):
        return (
            self.excavated_middle_clay_layer_thickness
            * self.excavation_width
            * self.saturated_unit_weight_of_clay
        )

    def calculate_weight_of_excavated_clay_lower_right(self):
        return (
            self.excavated_lower_clay_layer_thickness
            * self.excavation_width
            * self.saturated_unit_weight_of_clay
        )

    def calculate_weight_of_water_after_second_excavation(self):
        return (
            self.excavated_middle_clay_layer_thickness
            * self.excavation_width
            * self.unit_weight_of_water
        )

    def calculate_weight_of_water_after_third_excavation(self):
        return (
            (
                self.excavated_middle_clay_layer_thickness
                + self.excavated_lower_clay_layer_thickness
            )
            * self.excavation_width
            * self.unit_weight_of_water
        )

    def calculate_weight_of_diaphragm_wall(self):
        return self.weight_of_diaphragm_wall * self.length_of_diaphragm_wall

    def calculate_total_vertical_surface_load(self):
        return self.distributed_surface_load * self.load_edge_length

    def check_vertical_reaction(self, output_reader, project_path, stage, expected):
        filename = f"{stage['name']}.post.res"
        output_data = output_reader.read_output_from(
            os.path.join(project_path, filename)
        )
        rel_tolerance = 0.005
        self.assertAlmostEqual(
            self.total_reaction_y_from_output_data(
                output_data, stage["time"], self.BOTTOM_NODE_IDS
            ),
            expected,
            places=None,
            delta=rel_tolerance * expected,
        )

    def run_simulation(self, sub_directory_name):
        project_path = test_helper.get_file_path(
            os.path.join("submerged_construction_of_excavation", sub_directory_name)
        )

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for stage in self.STAGE.values():
                param_file = os.path.join("..", "common", f"{stage['name']}.json")
                with open(param_file, "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        output_reader = GiDOutputFileReader()
        expected_total_weight = self.calculate_weight_of_all_soil()
        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["initial_stage"],
            expected_total_weight,
        )

        self.check_vertical_reaction(
            output_reader, project_path, self.STAGE["null_step"], expected_total_weight
        )

        expected_total_weight += self.calculate_weight_of_diaphragm_wall()
        expected_total_vertical_reaction = (
            expected_total_weight + self.calculate_total_vertical_surface_load()
        )
        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["wall_installation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_upper_right()
        expected_total_vertical_reaction = (
            expected_total_weight + self.calculate_total_vertical_surface_load()
        )
        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["first_excavation"],
            expected_total_vertical_reaction,
        )

        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["strut_installation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_middle_right()
        expected_total_vertical_reaction = (
            expected_total_weight
            + self.calculate_total_vertical_surface_load()
            + self.calculate_weight_of_water_after_second_excavation()
        )
        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["second_excavation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_lower_right()
        expected_total_vertical_reaction = (
            expected_total_weight
            + self.calculate_total_vertical_surface_load()
            + self.calculate_weight_of_water_after_third_excavation()
        )
        self.check_vertical_reaction(
            output_reader,
            project_path,
            self.STAGE["third_excavation"],
            expected_total_vertical_reaction,
        )

    def test_simulation_with_linear_elastic_materials(self):
        self.run_simulation("linear_elastic")

if __name__ == "__main__":
    KratosUnittest.main()
