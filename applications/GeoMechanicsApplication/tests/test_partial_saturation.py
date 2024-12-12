# from KratosMultiphysics import * as Kratos

import sys
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsPartialSaturation(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def __test_saturated_below_phreatic_level_pw(self, test_name):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        from math import fabs
        #from analytical_solutions import calculate_column_saturated_under_phreatic_level
        
        # define number of stages
        n_stages = 2

        # get the parameter file names for all stages
        #test_name = test_name
        file_path = test_helper.get_file_path(os.path.join('.', 'test_partially_saturated', test_name))
        parameter_file_names = [os.path.join(file_path, 'ProjectParameters_stage' + str(i + 1) + '.json') for i in
                                range(n_stages)]

        # set stage parameters
        parameters_stages = [None] * n_stages
        os.chdir(file_path)
        for idx, parameter_file_name in enumerate(parameter_file_names):
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages[idx] = Kratos.Parameters(parameter_file.read())

        model = Kratos.Model()
        stages = [analysis.GeoMechanicsAnalysis(model, stage_parameters) for stage_parameters in parameters_stages]

        # run stages and get water pressure/displacement results per stage
        stage_water_pressure = [None] * n_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            stage_water_pressure[idx] = test_helper.get_water_pressure(stage)

        # get y coords of all the nodes
        coords = test_helper.get_nodal_coordinates(stages[1])
        y_coords = [coord[1] for coord in coords]
        
        

        # calculate water pressure analytical solution for all stages and calculate the error
        #rmse_stages = [None] * (n_stages - 1)
        rel_p_stage  = [self.__calculate_column_saturated_under_phreatic_level(y_coord) for y_coord in y_coords]
        print(rel_p_stage)
        
        errors_stage = [stage_water_pressure[1][node_idx] - rel_p for node_idx, rel_p in
                         enumerate(rel_p_stage)]
        rmse_stages = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 1.0e-7
        self.assertLess(rmse_stages, accuracy)

    
    def test_saturated_below_phreatic_level_pw_triangle3N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle3N')
        
    def test_saturated_below_phreatic_level_pw_triangle6N(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_pw_triangle6N')
        
    def test_saturated_below_phreatic_level_upw_difforder_triangle6n(self):
        self.__test_saturated_below_phreatic_level_pw('test_saturated_below_phreatic_level_upw_difforder_triangle6n')


    def __calculate_column_saturated_under_phreatic_level(self, y_coord):
        water_density = 1019.367991845056
        gravity = 9.81
        phreatic_level = -2.0
        result = gravity * water_density * (y_coord - phreatic_level)
        result = min([result, 0.0])
        return result
        return 0.0
        
        
        
if __name__ == '__main__':
    KratosUnittest.main()
