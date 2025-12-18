import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
from KratosMultiphysics.GeoMechanicsApplication.unit_conversions import unit_to_k_unit
import test_helper

import os
import json
from pathlib import Path

if test_helper.want_test_plots():
    import KratosMultiphysics.GeoMechanicsApplication.geo_plot_utilities as plot_utils


def get_wall_node_ids():
    return [
        4768,
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
        8988,
    ]

def get_interface_soil_node_ids():
    return [
        9424,
        9421,
        9416,
        9413,
        9406,
        9400,
        9396,
        9391,
        9386,
        9382,
        9377,
        9371,
        9365,
        9359,
        9355,
        9350,
        9345,
        9340,
        9333,
        9330,
        9408,
        9483,
        9479,
        9474,
        9469,
        9463,
        9459,
        9454,
        9449,
        9445,
        9442,
        9438,
        9431,
        9425,
        9419,
        9414,
        9407,
        9402,
        9395,
        9388,
        9378,
        9374,
        9370,
        9364,
        9358,
        9351,
        9344,
        9338,
        9331,
        9487,
        9484,
        9478,
        9473,
        9468,
        9464,
        9458,
        9452,
        9448,
        9443,
        9435,
        9433,
        9427,
        9420,
        9411,
        9403,
        9394,
        9389,
        9383,
        9375,
        9369,
        9363,
        9353,
        9346,
        9339,
        9334,
        9329,
        9485,
        9481,
        9476,
        9470,
        9465,
        9460,
        9455,
        9451,
        9447,
        9440,
        9436,
        9429,
        9423,
        9417,
        9412,
        9405,
        9399,
        9393,
        9385,
        9379,
        9372,
        9366,
        9360,
        9354,
        9348,
        9342,
        9336,
        9328,
        9486,
        9480,
        9475,
        9471,
        9466,
        9462,
        9457,
        9453,
        9444,
        9439,
        9434,
        9430,
        9426,
        9418,
        9410,
        9404,
        9397,
        9392,
        9387,
        9381,
        9373,
        9367,
        9362,
        9356,
        9349,
        9343,
        9337,
        9332,
        9488,
        9482,
        9477,
        9472,
        9467,
        9461,
        9456,
        9450,
        9446,
        9441,
        9437,
        9432,
        9428,
        9422,
        9415,
        9409,
        9401,
        9398,
        9390,
        9384,
        9380,
        9376,
        9368,
        9361,
        9357,
        9352,
        9347,
        9341,
        9335,
        9335,
        9491,
        9568,
        9566,
        9564,
        9560,
        9557,
        9554,
        9551,
        9548,
        9545,
        9542,
        9539,
        9535,
        9533,
        9529,
        9526,
        9523,
        9520,
        9517,
        9514,
        9510,
        9507,
        9505,
        9503,
        9500,
        9498,
        9496,
        9494,
        9492,
        9489,
        9565,
        9562,
        9559,
        9556,
        9552,
        9549,
        9546,
        9543,
        9540,
        9537,
        9534,
        9531,
        9528,
        9525,
        9522,
        9518,
        9515,
        9512,
        9509,
        9506,
        9504,
        9502,
        9501,
        9499,
        9497,
        9495,
        9493,
        9490,
        9567,
        9563,
        9561,
        9558,
        9555,
        9553,
        9550,
        9547,
        9544,
        9541,
        9538,
        9536,
        9532,
        9530,
        9527,
        9524,
        9521,
        9519,
        9516,
        9513,
        9511,
        9508,
    ]

def _extract_x_and_y_from_line(line, index_of_x=0, index_of_y=1, x_transform=None):
    words = line.split(",")
    x_ = float(words[index_of_x])
    if x_transform:
        x_ = x_transform(x_)
    y_ = float(words[index_of_y])

    return x_, y_


def extract_bending_moment_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=4, index_of_y=1)


def extract_axial_force_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=2, index_of_y=1)


def extract_shear_force_and_y_from_line(line):
    # The shear force sign in the comparison data is opposite to the Kratos sign
    return _extract_x_and_y_from_line(
        line, index_of_x=3, index_of_y=1, x_transform=lambda x: -x
    )


def extract_normal_traction_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=2, index_of_y=1)


def extract_shear_traction_and_y_from_line(line):
    return _extract_x_and_y_from_line(line, index_of_x=4, index_of_y=1)


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

        self.stages_info = {
            "initial_stage":      {"end_time": -1.0, "base_name": "1_Initial_stage"},
            "null_step":          {"end_time":  0.0, "base_name": "2_Null_step"},
            "wall_installation":  {"end_time":  1.0, "base_name": "3_Wall_installation"},
            "first_excavation":   {"end_time":  2.0, "base_name": "4_First_excavation"},
            "strut_installation": {"end_time":  3.0, "base_name": "5_Strut_installation"},
            "second_excavation":  {"end_time":  4.0, "base_name": "6_Second_excavation"},
            "third_excavation":   {"end_time":  5.0, "base_name": "7_Third_excavation"},
        }

        self.bottom_node_ids = [
            1,
            3,
            10,
            25,
            44,
            68,
            96,
            132,
            175,
            218,
            268,
            318,
            378,
            434,
            499,
            573,
            649,
            732,
            819,
            909,
            999,
            1102,
            1206,
            1315,
            1434,
            1560,
            1695,
            1828,
            1966,
            2111,
            2259,
            2406,
            2558,
            2725,
            2896,
            3059,
            3239,
            3422,
            3621,
            3818,
            4023,
            4233,
            4455,
            4739,
            5151,
            5563,
            5959,
            6281,
            6607,
            6924,
            7237,
            7541,
            7831,
        ]

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

    def check_vertical_reaction(
        self, project_path, stage, expected_total_vertical_reaction
    ):
        filename = f"{stage['base_name']}.post.res"
        output_data = GiDOutputFileReader().read_output_from(
            os.path.join(project_path, filename)
        )
        rel_tolerance = 0.005
        self.assertAlmostEqual(
            self.total_reaction_y_from_output_data(
                output_data, stage["end_time"], self.bottom_node_ids
            ),
            expected_total_vertical_reaction,
            places=None,
            delta=rel_tolerance * expected_total_vertical_reaction,
        )

    def run_simulation_and_checks(self, sub_directory_name):
        project_path = test_helper.get_file_path(
            os.path.join("submerged_construction_of_excavation", sub_directory_name)
        )

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for stage in self.stages_info.values():
                with open(
                    Path("..") / "common" / f"{stage['base_name']}.json", "r"
                ) as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        if test_helper.want_test_plots():
            self.create_wall_plots(project_path)
            self.create_interface_plots(project_path)

        expected_total_weight = self.calculate_weight_of_all_soil()
        self.check_vertical_reaction(
            project_path,
            self.stages_info["initial_stage"],
            expected_total_weight,
        )

        self.check_vertical_reaction(
            project_path, self.stages_info["null_step"], expected_total_weight
        )

        expected_total_weight += self.calculate_weight_of_diaphragm_wall()
        expected_total_vertical_reaction = (
            expected_total_weight + self.calculate_total_vertical_surface_load()
        )
        self.check_vertical_reaction(
            project_path,
            self.stages_info["wall_installation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_upper_right()
        expected_total_vertical_reaction = (
            expected_total_weight + self.calculate_total_vertical_surface_load()
        )
        self.check_vertical_reaction(
            project_path,
            self.stages_info["first_excavation"],
            expected_total_vertical_reaction,
        )

        self.check_vertical_reaction(
            project_path,
            self.stages_info["strut_installation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_middle_right()
        expected_total_vertical_reaction = (
            expected_total_weight
            + self.calculate_total_vertical_surface_load()
            + self.calculate_weight_of_water_after_second_excavation()
        )
        self.check_vertical_reaction(
            project_path,
            self.stages_info["second_excavation"],
            expected_total_vertical_reaction,
        )

        expected_total_weight -= self.calculate_weight_of_excavated_clay_lower_right()
        expected_total_vertical_reaction = (
            expected_total_weight
            + self.calculate_total_vertical_surface_load()
            + self.calculate_weight_of_water_after_third_excavation()
        )
        self.check_vertical_reaction(
            project_path,
            self.stages_info["third_excavation"],
            expected_total_vertical_reaction,
        )

    def get_structural_stages(self):
        return [
            stage
            for stage in self.stages_info.values()
            if stage["base_name"] != "1_Initial_stage"
            and stage["base_name"] != "2_Null_step"
        ]

    def create_interface_plots(self, project_path):
        structural_stages = self.get_structural_stages()

        plot_titles = [stage["base_name"] for stage in structural_stages]

        normal_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        normal_traction_plot_label = "Normal traction"
        normal_traction_collections = self.get_variable_collections_per_stage(
            normal_traction_kratos_label,
            normal_traction_plot_label,
            project_path,
            structural_stages,
            "left_interface",
            extract_normal_traction_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            normal_traction_collections,
            Path(project_path) / "normal_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel="Normal Traction [kN/m^2]",
            ylabel="y [m]",
            )

        shear_traction_plot_label = "Shear traction"
        shear_traction_kratos_label = "GEO_EFFECTIVE_TRACTION_VECTOR"
        shear_traction_collections = self.get_variable_collections_per_stage(
            shear_traction_kratos_label,
            shear_traction_plot_label,
            project_path,
            structural_stages,
            "left_interface",
            extract_shear_traction_and_y_from_line,
        )
        # sort on column with y
        #for collection in shear_traction_collections:
        #    print(f"{collection[0].label}: {collection[0].data_points}")
        #shear_traction_collections.sort(key=lambda i: (i[1], i[0]))
        plot_utils.make_sub_plots(
            shear_traction_collections,
            Path(project_path) / "shear_tractions_all_stages.svg",
            titles=plot_titles,
            xlabel="Shear Traction [kN/m^2]",
            ylabel="y [m]",
            )

    def create_wall_plots(self, project_path):
        structural_stages = self.get_structural_stages()

        plot_titles = [stage["base_name"] for stage in structural_stages]

        axial_force_plot_label = "Axial force"
        axial_force_kratos_label = "AXIAL_FORCE"
        axial_force_collections = self.get_variable_collections_per_stage(
            axial_force_kratos_label,
            axial_force_plot_label,
            project_path,
            structural_stages,
            "wall",
            extract_axial_force_and_y_from_line,
        )

        plot_utils.make_sub_plots(
            axial_force_collections,
            Path(project_path) / "axial_forces_all_stages.svg",
            titles=plot_titles,
            xlabel="Axial Force [kN/m]",
            ylabel="y [m]",
        )

        shear_force_plot_label = "Shear force"
        shear_force_kratos_label = "SHEAR_FORCE"
        shear_force_collections = self.get_variable_collections_per_stage(
            shear_force_kratos_label,
            shear_force_plot_label,
            project_path,
            structural_stages,
            "wall",
            extract_shear_force_and_y_from_line,
        )
        plot_utils.make_sub_plots(
            shear_force_collections,
            Path(project_path) / "shear_forces_all_stages.svg",
            titles=plot_titles,
            xlabel="Shear Force [kN/m]",
            ylabel="y [m]",
        )

        bending_moment_plot_label = "Bending moment"
        bending_moment_kratos_label = "BENDING_MOMENT"
        bending_moment_collections = self.get_variable_collections_per_stage(
            bending_moment_kratos_label,
            bending_moment_plot_label,
            project_path,
            structural_stages,
            "wall",
            extract_bending_moment_and_y_from_line,
        )
        plot_utils.make_sub_plots(
            bending_moment_collections,
            Path(project_path) / "bending_moments_all_stages.svg",
            titles=plot_titles,
            xlabel="Bending moment [kNm/m]",
            ylabel="y [m]",
        )

    def get_variable_collections_per_stage(
        self,
        kratos_variable_label,
        variable_plot_label,
        project_path,
        structural_stages,
        object_name,
        data_point_extractor,
    ):
        node_ids = []
        if object_name.__contains__("wall"):
            node_ids = get_wall_node_ids()
        else:
            node_ids = get_interface_soil_node_ids()

        # Since the coordinates do not change between stages, we base them on the first stage
        y_coords = self.get_y_coords(
            project_path, structural_stages[0]["base_name"], node_ids
        )

        variable_data_collections = []
        for stage in structural_stages:
            output_data = GiDOutputFileReader().read_output_from(
                os.path.join(project_path, f"{stage['base_name']}.post.res")
            )
            if object_name.__contains__("wall"):
                variable_kratos_data = GiDOutputFileReader.nodal_values_at_time(
                    kratos_variable_label,
                    stage["end_time"],
                    output_data,
                    node_ids=node_ids,
                )
            else:
                json_data = self.read_json_output(project_path, stage)
                variable_kratos_data = []
                for node_label in [f"NODE_{node_id}" for node_id in node_ids]:
                    variable_kratos_data.append(json_data[node_label][kratos_variable_label][0][1])

            variable_kratos_data = [
                unit_to_k_unit(value) for value in variable_kratos_data
            ]

            data_series_collection = []
            data_series_collection.append(
                plot_utils.DataSeries(
                    zip(variable_kratos_data, y_coords),
                    f"{variable_plot_label} [Kratos]",
                    line_style="-",
                    marker=".",
                )
            )
            path_to_comparison_file = (
                Path(project_path)
                / "comparison_data"
                / f"{stage['base_name']}_{object_name}_comparison.csv"
            )
            comparison_variable = test_helper.get_data_points_from_file(
                path_to_comparison_file, data_point_extractor
            )
            data_series_collection.append(
                plot_utils.DataSeries(
                    comparison_variable,
                    f"{variable_plot_label} [Comparison]",
                    marker="1",
                )
            )
            variable_data_collections.append(data_series_collection)
        return variable_data_collections

    def read_json_output(self, project_path, stage):
        with open(
                os.path.join(project_path, f"{stage['base_name']}_interface_output.json"), "r"
        ) as output_file:
            return json.load(output_file)

    def get_y_coords(self, project_path, stage_base_name, node_ids):
        coordinates = test_helper.read_coordinates_from_post_msh_file(
            Path(project_path) / f"{stage_base_name}.post.msh", node_ids=node_ids
        )
        return [coord[1] for coord in coordinates]

    def test_simulation_with_linear_elastic_materials(self):
        self.run_simulation_and_checks("linear_elastic")


if __name__ == "__main__":
    KratosUnittest.main()
