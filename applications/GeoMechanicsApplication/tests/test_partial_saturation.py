import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsPartialSaturation(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def __test_saturated_below_phreatic_level_pw(self, test_name):
        n_stages = 2

        # get the parameter file names for all stages
        file_path = test_helper.get_file_path(os.path.join('test_partially_saturated', test_name))
        parameter_file_names = [os.path.join(file_path, f'ProjectParameters_stage{i+1}.json') for i in
                                range(n_stages)]

        # set stage parameters
        parameters_stages = []
        initial_directory = os.getcwd()
        os.chdir(file_path)
        for parameter_file_name in parameter_file_names:
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages.append(Kratos.Parameters(parameter_file.read()))

        model = Kratos.Model()

        # run stages and get water pressure/displacement results per stage
        stage_water_pressure = []
        coords = []
        for stage_parameters in parameters_stages:
            stage = analysis.GeoMechanicsAnalysis(model, stage_parameters)
            stage.Run()
            stage_water_pressure.append(test_helper.get_water_pressure(stage))
            coords.append(test_helper.get_nodal_coordinates(stage))

        os.chdir(initial_directory)
        # get y coords of all the nodes
        y_coords = [coord[1] for coord in coords[0]]

        # calculate water pressure analytical solution for all stages and calculate the error
        rel_p_stage  = [self.__compute_hydrostatic_water_pressure(y_coord) for y_coord in y_coords]
        
        errors_stage = [actual_pressure - expected_pressure for actual_pressure, expected_pressure in
                        zip(stage_water_pressure[1], rel_p_stage)]
        rmse_stages = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below accuracy
        accuracy = 1.0e-7
        self.assertLess(rmse_stages, accuracy)

    def test_saturated_below_phreatic_level_pw_triangle3N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle3N')
        
    def test_saturated_below_phreatic_level_pw_triangle6N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle6N')
        
    def test_saturated_below_phreatic_level_upw_difforder_triangle6n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_difforder_triangle6n')
        
    def test_saturated_below_phreatic_level_upw_smallstrain_triangle3n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_smallstrain_triangle3n')
        
    def test_saturated_below_phreatic_level_upw_smallstrain_triangle6n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_smallstrain_triangle6n')

    def __compute_hydrostatic_water_pressure(self, y_coord):
        water_density = 1019.367991845056
        gravity = -9.81
        phreatic_level = -2.0
        result = gravity * water_density * (phreatic_level - y_coord)
        return min([result, 0.0])



    def define_expected_water_pressures(self):
        self.expected_water_pressures = [{"output_filename": "saturatedbelowphreaticlevel.post.res",
                                          "time": 1.0,
                                          "expected_values": [
                                              {"node": 1, "WATER_PRESSURE": -10000.0},
                                              {"node": 4, "WATER_PRESSURE": -10000.0},
                                              {"node": 5, "WATER_PRESSURE": 0.0},
                                              {"node": 18, "WATER_PRESSURE": 0.0}
                                          ]},
                                         {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                          "time": 4.0,
                                          "expected_values": [
                                              {"node": 1, "WATER_PRESSURE": -20000.0},
                                              {"node": 4, "WATER_PRESSURE": -20000.0},
                                              {"node": 5, "WATER_PRESSURE": -10000.0},
                                              {"node": 18, "WATER_PRESSURE": -10000.0}
                                          ]},
                                         {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                          "time": 7.0,
                                          "expected_values": [
                                              {"node": 1, "WATER_PRESSURE": -30000.0},
                                              {"node": 4, "WATER_PRESSURE": -30000.0},
                                              {"node": 5, "WATER_PRESSURE": -20000.0},
                                              {"node": 18, "WATER_PRESSURE": -20000.0}
                                          ]},
                                         {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                          "time": 10.0,
                                          "expected_values": [
                                              {"node": 1, "WATER_PRESSURE": -40000.0},
                                              {"node": 4, "WATER_PRESSURE": -40000.0},
                                              {"node": 5, "WATER_PRESSURE": -30000.0},
                                              {"node": 18, "WATER_PRESSURE": -30000.0}
                                          ]}]

    def define_expected_bishop_coefficients(self):
        self.expected_bishop_coefficients = [{"output_filename": "saturatedbelowphreaticlevel.post.res",
                                              "time": 1.0,
                                              "expected_values": [
                                                  {"element": 1, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 2, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 3, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 4, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 5, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]}
                                              ]},
                                             {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                              "time": 4.0,
                                              "expected_values": [
                                                  {"element": 1, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 2, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 3, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 4, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 5, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]}
                                              ]},
                                             {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                              "time": 7.0,
                                              "expected_values": [
                                                  {"element": 1, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 2, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 3, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 4, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]},
                                                  {"element": 5, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]}
                                              ]},
                                             {"output_filename": "saturatedbelowphreaticlevel.post.res",
                                              "time": 10.0,
                                              "expected_values": [
                                                  {"element": 1, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 2, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 3, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 4, "BISHOP_COEFFICIENT": [1.0, 1.0, 1.0, 1.0]},
                                                  {"element": 5, "BISHOP_COEFFICIENT": [0.0, 0.0, 0.0, 0.0]}
                                              ]}]
    def check_water_pressures(self, test_path):
        reader = test_helper.GiDOutputFileReader()

        for item in self.expected_water_pressures:
            time = item["time"]
            node_ids = [sub_item["node"] for sub_item in item["expected_values"]]
            expected_water_pressures = [sub_item["WATER_PRESSURE"] for sub_item in item["expected_values"]]

            actual_data = reader.read_output_from(os.path.join(test_path, item["output_filename"]))
            actual_water_pressures = reader.nodal_values_at_time("WATER_PRESSURE", time, actual_data, node_ids)

            self.assertEqual(len(actual_water_pressures), len(expected_water_pressures))
            for actual_water_pressure, expected_water_pressure in zip(actual_water_pressures, expected_water_pressures):
                self.assertAlmostEqual(actual_water_pressure, expected_water_pressure, 1)

    def check_Bishop_coefficients(self, test_path):
        reader = test_helper.GiDOutputFileReader()

        for item in self.expected_bishop_coefficients:
            time = item["time"]
            expected_bishop_coefficients = [sub_item["BISHOP_COEFFICIENT"] for sub_item in item["expected_values"]]

            actual_data = reader.read_output_from(os.path.join(test_path, item["output_filename"]))
            actual_bishop_coefficients = reader.element_integration_point_values_at_time("BISHOP_COEFFICIENT", time, actual_data, [1,2,3,4,5], [0,1,2,3])
            self.assertEqual(len(actual_bishop_coefficients), len(expected_bishop_coefficients))
            for actual_bishop_coefficient, expected_bishop_coefficient in zip(actual_bishop_coefficients, expected_bishop_coefficients):
                self.assertVectorAlmostEqual(actual_bishop_coefficient, expected_bishop_coefficient)

    def test_rising_water_quad4N(self):
        test_path = test_helper.get_file_path(os.path.join('test_partially_saturated', 'rising_water_pw_quad4N'))

        # The expected values are analytical results
        self.define_expected_water_pressures()
        self.define_expected_bishop_coefficients()
        test_helper.run_kratos(test_path)
        self.check_water_pressures(test_path)
        self.check_Bishop_coefficients(test_path)

if __name__ == '__main__':
    KratosUnittest.main()
