import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

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
        rel_p_stage  = [self.__compute_hydrostatic_water_pressure(y_coord, -2.0) for y_coord in y_coords]
        
        errors_stage = [actual_pressure - expected_pressure for actual_pressure, expected_pressure in
                        zip(stage_water_pressure[1], rel_p_stage)]
        rmse_stages = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below accuracy
        accuracy = 1.0e-3
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

    def test_climbing_falling_phreatic_level_upw_smallstrain_quad4n(self):
        # only waterpressures below phreatic level are checked with an analytical solution.
        # values above phreatic level give suction of an unchecked amount.
        file_path = test_helper.get_file_path(os.path.join('test_partially_saturated', 'test_rising_falling_phreatic_level_pw_quad4N'))
        simulation = test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(os.path.join(file_path, 'rising_falling_phreatic_level_pw_quad4n.post.res'))
        coords = test_helper.get_nodal_coordinates(simulation)
        times = [1.0, 5.0, 9.0, 13.0, 17.0, 21.0, 25.0, 29.0]
        water_levels = [-4.0, -3.0, -2.0, -1.0, -2.0, -3.0, -4.0, -5.0]
        for time, water_level in zip(times, water_levels):
            water_pressures = reader.nodal_values_at_time('WATER_PRESSURE', time, output_data)
            negative_water_pressures = [min([water_pressure, 0.0]) for water_pressure in water_pressures]
            analytical_water_pressures = [self.__compute_hydrostatic_water_pressure(coord[1], water_level) for coord in coords]
            self.assertVectorAlmostEqual(negative_water_pressures, analytical_water_pressures, places=None, msg=f"water pressures at time {time}", delta=10.0)

    def __compute_hydrostatic_water_pressure(self, y_coord, phreatic_level):
        water_weight = -10000.0
        result = water_weight * (phreatic_level - y_coord)
        return min([result, 0.0])

if __name__ == '__main__':
    KratosUnittest.main()
