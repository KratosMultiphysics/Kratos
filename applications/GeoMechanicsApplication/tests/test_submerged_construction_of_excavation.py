import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
import test_helper

import os


class KratosGeoMechanicsSubmergedConstructionOfExcavation(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.unsaturated_unit_weight_of_clay = 16.0e3  # N/m3
        self.saturated_unit_weight_of_clay = 18.0e3  # N/m3
        self.saturated_unit_weight_of_sand = 20.0e3  # N/m3
        self.model_width = 65.0  # m
        self.unsaturated_clay_layer_thickness = 2.0  # m
        self.saturated_clay_layer_thickness = 18.0  # m
        self.saturated_sand_layer_thickness = 30.0  # m
        self.weight_of_diaphragm_wall = 10.0e3  # N / m
        self.length_of_diaphragm_wall = 30.0  # m

    def total_reaction_y_from_output_data(self, output_data, time, node_ids):
        reactions = test_helper.GiDOutputFileReader.nodal_values_at_time("REACTION", time, output_data, node_ids=node_ids)
        return sum([reaction[1] for reaction in reactions])

    def calculate_weight_of_all_soil(self):
        return self.model_width * (self.unsaturated_clay_layer_thickness * self.unsaturated_unit_weight_of_clay +
                                   self.saturated_clay_layer_thickness * self.saturated_unit_weight_of_clay +
                                   self.saturated_sand_layer_thickness * self.saturated_unit_weight_of_sand)

    def calculate_weight_of_diaphragm_wall(self):
        return self.weight_of_diaphragm_wall * self.length_of_diaphragm_wall

    def test_run_simulation(self):
        project_path = test_helper.get_file_path("submerged_construction_of_excavation")
        project_parameters_filenames = ["1_Initial_stage.json", "2_Null_step.json", "3_Wall_installation.json"]

        with context_managers.set_cwd_to(project_path):
            model = Kratos.Model()
            for filename in project_parameters_filenames:
                with open(filename, "r") as f:
                    stage_parameters = Kratos.Parameters(f.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

        # Check vertical reaction forces in the initial stage
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(os.path.join(project_path, "1_Initial_stage.post.res"))
        time = -1.0
        bottom_node_ids = [1, 4, 10, 25, 43, 69, 97, 133, 174, 218, 269, 318, 377, 434, 500, 572, 648, 732, 820, 909, 1000, 1101, 1206, 1316, 1433, 1561, 1695, 1829, 1965, 2110, 2258, 2406, 2559, 2725, 2896, 3058, 3238, 3423, 3620, 3818, 4022, 4233, 4455, 4739, 5151, 5563, 5959, 6281, 6607, 6924, 7236, 7541, 7831]
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
        self.assertAlmostEqual(self.total_reaction_y_from_output_data(output_data, time, bottom_node_ids),
                               expected_total_weight, places=None, delta=rel_tolerance*expected_total_weight)


if __name__ == "__main__":
    KratosUnittest.main()
