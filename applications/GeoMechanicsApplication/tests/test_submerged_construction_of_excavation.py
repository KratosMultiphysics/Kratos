import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (GiDOutputFileReader, )
import test_helper

import os
from pathlib import Path
if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils
def get_wall_node_ids():
    return [4768,
            4782,
            4794,
            4811,
            4830,
            4845,
            4863,
            4878,
            4896,
            4912,
            4933,
            4946,
            4965,
            4981,
            4995,
            5011,
            5028,
            5039,
            5059,
            5074,
            5091,
            5113,
            5131,
            5147,
            5162,
            5184,
            5201,
            5220,
            5236,
            5253,
            5269,
            5286,
            5301,
            5322,
            5341,
            5362,
            5377,
            5395,
            5413,
            5432,
            5449,
            5467,
            5486,
            5499,
            5521,
            5540,
            5559,
            5579,
            5596,
            5613,
            5632,
            5654,
            5672,
            5688,
            5707,
            5732,
            5749,
            5766,
            5783,
            5803,
            5819,
            5837,
            5857,
            5877,
            5897,
            5916,
            5939,
            5957,
            5977,
            5991,
            6007,
            6029,
            6042,
            6060,
            6078,
            6096,
            6118,
            6128,
            6142,
            6161,
            6173,
            6193,
            6210,
            6223,
            6238,
            6249,
            6268,
            6289,
            6305,
            6317,
            6333,
            6356,
            6370,
            6390,
            6405,
            6425,
            6441,
            6456,
            6473,
            6488,
            6505,
            6524,
            6539,
            6559,
            6582,
            6600,
            6619,
            6635,
            6654,
            6675,
            6692,
            6712,
            6728,
            6748,
            6762,
            6778,
            6797,
            6821,
            6837,
            6851,
            6867,
            6887,
            6903,
            6922,
            6944,
            6961,
            6980,
            6993,
            7014,
            7037,
            7055,
            7075,
            7089,
            7105,
            7123,
            7139,
            7162,
            7175,
            7199,
            7218,
            7236,
            7252,
            7272,
            7289,
            7306,
            7327,
            7342,
            7364,
            7387,
            7402,
            7417,
            7434,
            7451,
            7469,
            7495,
            7510,
            7530,
            7549,
            7561,
            7578,
            7597,
            7617,
            7636,
            7653,
            7668,
            7687,
            7706,
            7721,
            7741,
            7764,
            7787,
            7807,
            7824,
            7849,
            7868,
            7886,
            7909,
            7928,
            7952,
            7971,
            7994,
            8006,
            8027,
            8044,
            8062,
            8082,
            8106,
            8121,
            8141,
            8157,
            8175,
            8188,
            8205,
            8230,
            8250,
            8268,
            8289,
            8305,
            8321,
            8334,
            8351,
            8378,
            8399,
            8419,
            8432,
            8448,
            8464,
            8479,
            8501,
            8522,
            8539,
            8552,
            8567,
            8585,
            8600,
            8616,
            8634,
            8655,
            8673,
            8687,
            8698,
            8713,
            8725,
            8749,
            8764,
            8778,
            8790,
            8806,
            8822,
            8836,
            8852,
            8869,
            8882,
            8898,
            8913,
            8928,
            8940,
            8956,
            8968,
            8979,
            8988]

def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1):
    words = line.split(',')
    return (-1.0 * float(words[index_of_x]), float(words[index_of_y]))


def extract_y_and_moment_from_line(line):
    words = line.split(',')
    return (float(words[11]), float(words[4]))
    # return _extract_x_and_y_from_line(line, index_of_x=11, index_of_y=4)

def extract_y_and_axial_force_from_line(line):
    words = line.split(',')
    return (float(words[5]), float(words[4]))

def extract_y_and_shear_force_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=8, index_of_y=4)


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
        self.weight_of_diaphragm_wall = 10.0e3  # N/m
        self.length_of_diaphragm_wall = 30.0  # m
        self.distributed_surface_load = 5.0e3  # N/m/m
        self.load_edge_length = 5.0  # m

        self.stages_info = {"initial_stage": {"end_time": -1.0, "base_name": "1_Initial_stage"},
                            "null_step": {"end_time": 0.0, "base_name": "2_Null_step"},
                            "wall_installation": {"end_time": 1.0, "base_name": "3_Wall_installation"},
                            "first_excavation": {"end_time": 2.0, "base_name": "4_First_excavation"},
                            "strut_installation": {"end_time": 3.0, "base_name": "5_Strut_installation", },
                            "second_excavation": {"end_time": 4.0, "base_name": "6_Second_excavation"},
                            "third_excavation": {"end_time": 5.0, "base_name": "7_Third_excavation"}, }

        self.bottom_node_ids = [1, 3, 10, 25, 44, 68, 96, 132, 175, 218, 268, 318, 378, 434, 499, 573, 649, 732, 819,
                                909, 999, 1102, 1206, 1315, 1434, 1560, 1695, 1828, 1966, 2111, 2259, 2406, 2558, 2725,
                                2896, 3059, 3239, 3422, 3621, 3818, 4023, 4233, 4455, 4739, 5151, 5563, 5959, 6281,
                                6607, 6924, 7237, 7541, 7831, ]

    def total_reaction_y_from_output_data(self, output_data, time, node_ids):
        reactions = GiDOutputFileReader.nodal_values_at_time("REACTION", time, output_data, node_ids=node_ids)
        return sum([reaction[1] for reaction in reactions])

    def calculate_weight_of_all_soil(self):
        return self.model_width * (
                self.unsaturated_clay_layer_thickness * self.unsaturated_unit_weight_of_clay +
                self.saturated_clay_layer_thickness * self.saturated_unit_weight_of_clay +
                self.saturated_sand_layer_thickness * self.saturated_unit_weight_of_sand)

    def calculate_weight_of_excavated_clay_upper_right(self):
        return (self.unsaturated_clay_layer_thickness * self.excavation_width * self.unsaturated_unit_weight_of_clay)

    def calculate_weight_of_excavated_clay_middle_right(self):
        return (self.excavated_middle_clay_layer_thickness * self.excavation_width * self.saturated_unit_weight_of_clay)

    def calculate_weight_of_excavated_clay_lower_right(self):
        return (self.excavated_lower_clay_layer_thickness * self.excavation_width * self.saturated_unit_weight_of_clay)

    def calculate_weight_of_water_after_second_excavation(self):
        return (self.excavated_middle_clay_layer_thickness * self.excavation_width * self.unit_weight_of_water)

    def calculate_weight_of_water_after_third_excavation(self):
        return ((
                            self.excavated_middle_clay_layer_thickness + self.excavated_lower_clay_layer_thickness) *
                self.excavation_width * self.unit_weight_of_water)

    def calculate_weight_of_diaphragm_wall(self):
        return self.weight_of_diaphragm_wall * self.length_of_diaphragm_wall

    def calculate_total_vertical_surface_load(self):
        return self.distributed_surface_load * self.load_edge_length

    def check_vertical_reaction(self, project_path, stage, expected_total_vertical_reaction):
        filename = f"{stage['base_name']}.post.res"
        output_data = GiDOutputFileReader().read_output_from(os.path.join(project_path, filename))
        rel_tolerance = 0.005
        self.assertAlmostEqual(
            self.total_reaction_y_from_output_data(output_data, stage["end_time"], self.bottom_node_ids),
            expected_total_vertical_reaction, places=None, delta=rel_tolerance * expected_total_vertical_reaction, )

    def run_simulation_and_checks(self, sub_directory_name):
        project_path = test_helper.get_file_path(
            os.path.join("submerged_construction_of_excavation", sub_directory_name))

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for stage in self.stages_info.values():
                with open(Path("..") / "common" / f"{stage['base_name']}.json", "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        if test_helper.want_test_plots():
            structural_stages=[stage for stage in self.stages_info.values() if
                               stage['base_name'] != "1_Initial_stage" and stage['base_name'] != "2_Null_step"
            ]
            titles = self.get_names_of_structural_stages(structural_stages)
            y_coords_per_stage = self.get_y_coords_per_stage(project_path, structural_stages)

            axial_force_plot_label = "Axial force"
            axial_force_kratos_label = "AXIAL_FORCE"
            axial_force_collections = self.get_variable_collections_per_stage(axial_force_kratos_label, axial_force_plot_label,
                                                                              project_path, structural_stages, y_coords_per_stage, extract_y_and_axial_force_from_line)

            plot_utils.make_sub_plots(
                axial_force_collections,
                Path(project_path) / f"axial_forces_all_stages.svg",
                xlabel="Axial Force [kN]",
                ylabel="y [m]",
                titles=titles
            )

            shear_force_plot_label = "Shear force"
            shear_force_kratos_label = "SHEAR_FORCE"
            shear_force_collections = self.get_variable_collections_per_stage(shear_force_kratos_label, shear_force_plot_label,
                                                                              project_path, structural_stages, y_coords_per_stage, extract_y_and_shear_force_from_line)
            plot_utils.make_sub_plots(
                shear_force_collections,
                Path(project_path) / f"shear_forces_all_stages.svg",
                xlabel="Shear Force [kN]",
                ylabel="y [m]",
                titles=titles
            )

            bending_moment_plot_label = "Bending moment"
            bending_moment_kratos_label = "BENDING_MOMENT"
            bending_moment_collections = self.get_variable_collections_per_stage(bending_moment_kratos_label, bending_moment_plot_label,
                                                                              project_path, structural_stages, y_coords_per_stage, extract_y_and_moment_from_line)
            plot_utils.make_sub_plots(
                bending_moment_collections,
                Path(project_path) / f"bending_moments_all_stages.svg",
                xlabel="Bending moment [kNm]",
                ylabel="y [m]",
                titles=titles
            )

        expected_total_weight = self.calculate_weight_of_all_soil()
        self.check_vertical_reaction(project_path, self.stages_info["initial_stage"], expected_total_weight, )

        self.check_vertical_reaction(project_path, self.stages_info["null_step"], expected_total_weight)

        expected_total_weight += self.calculate_weight_of_diaphragm_wall()
        expected_total_vertical_reaction = (expected_total_weight + self.calculate_total_vertical_surface_load())
        self.check_vertical_reaction(project_path, self.stages_info["wall_installation"],
                                     expected_total_vertical_reaction, )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_upper_right()
        expected_total_vertical_reaction = (expected_total_weight + self.calculate_total_vertical_surface_load())
        self.check_vertical_reaction(project_path, self.stages_info["first_excavation"],
                                     expected_total_vertical_reaction, )

        self.check_vertical_reaction(project_path, self.stages_info["strut_installation"],
                                     expected_total_vertical_reaction, )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_middle_right()
        expected_total_vertical_reaction = (
                expected_total_weight + self.calculate_total_vertical_surface_load() +
                self.calculate_weight_of_water_after_second_excavation())
        self.check_vertical_reaction(project_path, self.stages_info["second_excavation"],
                                     expected_total_vertical_reaction, )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_lower_right()
        expected_total_vertical_reaction = (
                expected_total_weight + self.calculate_total_vertical_surface_load() +
                self.calculate_weight_of_water_after_third_excavation())
        self.check_vertical_reaction(project_path, self.stages_info["third_excavation"],
                                     expected_total_vertical_reaction, )

    def get_variable_collections_per_stage(self, axial_force_kratos_label, axial_force_plot_label,
                                           project_path,
                                           structural_stages, y_coords_per_stage, data_point_extractor):
        axial_force_collections = []
        for stage, y_coords in zip(structural_stages, y_coords_per_stage):
            data_series_collection = []
            output_data_wall = GiDOutputFileReader().read_output_from(
                os.path.join(project_path, f"{stage['base_name']}.post.res"))
            axial_forces = GiDOutputFileReader.nodal_values_at_time(axial_force_kratos_label, stage['end_time'],
                                                                    output_data_wall,
                                                                    node_ids=get_wall_node_ids())
            axial_forces = [axial_force / 1000 for axial_force in axial_forces]  # Convert to kN

            data = []
            for force, y in zip(axial_forces, y_coords):
                data.append((force, y))

            data_series_collection.append(plot_utils.DataSeries(
                data,
                f"{axial_force_plot_label} [Kratos]",
                line_style="-",
                marker=".",
            ))
            path_to_comparison_file = Path(project_path) / "comparison_data" / f"{stage['base_name']}_comparison.csv"
            comparison_axial_force = test_helper.get_data_points_from_file(
                path_to_comparison_file,
                data_point_extractor
            )
            data_series_collection.append(
                plot_utils.DataSeries(comparison_axial_force, f"{axial_force_plot_label} [Comparison]", marker="1")
            )
            axial_force_collections.append(data_series_collection)
        return axial_force_collections

    def get_y_coords_per_stage(self, project_path, structural_stages):
        y_coords_per_stage = []
        for stage in structural_stages:
            node_ids = get_wall_node_ids()
            coordinates = test_helper.read_coordinates_from_post_msh_file(
                Path(project_path) / f"{stage['base_name']}.post.msh", node_ids=node_ids
            )
            y_coords_per_stage.append([coord[1] for coord in coordinates])
        return y_coords_per_stage

    def get_names_of_structural_stages(self, structural_stages):
        return [stage['base_name'] for stage in structural_stages]

    def test_simulation_with_linear_elastic_materials(self):
        self.run_simulation_and_checks("linear_elastic")


if __name__ == "__main__":
    KratosUnittest.main()
