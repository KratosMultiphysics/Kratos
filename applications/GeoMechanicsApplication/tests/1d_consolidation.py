# from KratosMultiphysics import * as Kratos

import sys
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsBenchmarkSet1(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass
        
    def test_1d_consolidation(self):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        from math import fabs
        from analytical_solutions import calculate_pore_pressure_1D_consolidation, calculate_degree_of_1D_consolidation
        
        # define number of stages
        n_stages = 11

        # get the parameter file names for all stages
        test_name = '1D-Consolidation_all_stages'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
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
        stage_displacement   = [None] * n_stages
        for idx, stage in enumerate(stages):
            stage.Run()
            stage_water_pressure[idx] = test_helper.get_water_pressure(stage)
            displacements = test_helper.get_displacement(stage)
            stage_displacement[idx]   = [displacement[1] for displacement in displacements] 

        # get y coords of all the nodes
        coords = test_helper.get_nodal_coordinates(stages[0])
        y_coords = [coord[1] + 1 for coord in coords]
        
        t_vs = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
        # calculate water pressure analytical solution for all stages and calculate the error
        sample_height = 1
        rmse_stages = [None] * (n_stages - 1)
        for idx, t_v in enumerate(t_vs):
            rel_p_stage  = [calculate_pore_pressure_1D_consolidation(y_coord, sample_height, t_v) * -1 for y_coord in y_coords]
            errors_stage = [stage_water_pressure[idx + 1][node_idx] - rel_p for node_idx, rel_p in
                            enumerate(rel_p_stage)]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for rmse_stage in rmse_stages:
            self.assertLess(rmse_stage, accuracy)

        # calculate the degree of consolidation analytical solution for all stages and calculate the error
        # Veruijt's notations
        q = -1.0
        K = 3.33e2
        G = 5.0e2
        m_v = 1/(K + 4/3 * G)
        delta_h8 = (-1)*m_v*sample_height*q
        beta = 0.5e-9
        n = 0.3
        delta_h0 = (-1)*m_v*sample_height*q*n*beta/(m_v+n*beta)
        for idx, t_v in enumerate(t_vs):
            rel_displacement = [(stage_displacement[idx + 1][node_idx]-delta_h0)/(delta_h8-delta_h0) for node_idx, y_coord in
                            enumerate(y_coords) if fabs(y_coord - sample_height) < 0.001]
            degree_of_consolidation = [calculate_degree_of_1D_consolidation(sample_height, t_v) * -1 
                            for y_coord in y_coords if fabs(y_coord - sample_height) < 0.001]
            errors_stage = [rel_displacement[node_idx] - degree_of_consolidation[node_idx] for node_idx in range(len(degree_of_consolidation))]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for rmse_stage in rmse_stages:
            self.assertLess(rmse_stage, accuracy)

if __name__ == '__main__':
    KratosUnittest.main()
