# from KratosMultiphysics import * as Kratos

import sys
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

import matplotlib.pyplot as plt
import numpy as np


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
        
    @KratosUnittest.skip("Test test_benchmark1_1 skipped temporary.")
    def test_benchmark1_1(self):
        """
        In this benchmark Bia-axial shearing test conditions are tested in the Kratos-Geomechanics application.
        In that way, this example can be used to verify elastic deformation with a linear elastic model.
        :return:
        """
        
        test_name = 'Biaxialshearstresswithlinearelasticmodel'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))

        simulation = test_helper.run_kratos(file_path)

        strain_tensor = test_helper.get_green_lagrange_strain_tensor(simulation)

        # check strain_tensor on all gauss points
        for i in range(len(strain_tensor)):
            for j in range(len(strain_tensor[i])):
                self.assertAlmostEqual(-0.00125, strain_tensor[i][j][0,1])

    #@KratosUnittest.skip("Test skipped, because it is unreliable. To be enabled after investigation of the root cause.")
    def test_benchmark1_4(self):
        """
        test 1D consolidation on elastic soil.

        :return:
        """
        from math import fabs
        from analytical_solutions import calculate_pore_pressure_1D_consolidation, calculate_total_deformation_1D_consolidation
        
        print("**************test_benchmark1_4******************************************")

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

        print("**************Creating pictures******************************************")
        # creating a water pressure picture in svg format
        plt.figure(figsize=(12, 10))
        plt.ylabel('height', fontsize=18)
        plt.title('Water pressure', fontsize=22)
        for i, curve in enumerate(stage_water_pressure):
            x_values = y_coords.copy()
            y_values = curve.copy()
            x_values, y_values = zip(*sorted(zip(x_values, y_values)))
            plt.plot(y_values, x_values, color=np.random.rand(3, ), marker="o", label=f"Stage #{i}", )

        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.savefig("water_pressure.svg", format="svg")
        plt.show
        
        plt.figure(figsize=(12, 10))
        plt.ylabel('height', fontsize=18)
        plt.title('Displacement', fontsize=22)
        for i, curve in enumerate(stage_displacement):
            x_values = y_coords.copy()
            y_values = curve.copy()
            x_values, y_values = zip(*sorted(zip(x_values, y_values)))
            plt.plot(y_values, x_values, color=np.random.rand(3, ), marker="o", label=f"Stage #{i}", )

        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.savefig("displacement.svg", format="svg")
        plt.show       

        print("**************Saved the pictures******************************************")

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

        # calculate total deformation analytical solution for all stages and calculate the error
        delta_h0 = 0
        delta_h8 = 0.001
        for idx, t_v in enumerate(t_vs):
            rel_disp   = [(stage_displacement[idx + 1][node_idx]-delta_h0)/(delta_h8-delta_h0) for node_idx, y_coord in
                            enumerate(y_coords) if fabs(y_coord - sample_height) < 0.001]
            tot_defomr = [calculate_total_deformation_1D_consolidation(sample_height, t_v) * -1 
                            for y_coord in y_coords if fabs(y_coord - sample_height) < 0.001]
            errors_stage = [rel_disp[node_idx] - tot_defomr[node_idx] for node_idx in range(len(tot_defomr))]
            rmse_stages[idx] = (sum([error ** 2 for error in errors_stage]) / len(errors_stage)) ** 0.5

        # assert if average error in all stages is below 1 percent
        accuracy = 0.01
        for rmse_stage in rmse_stages:
            print("rmse_stage",rmse_stage)
            self.assertLess(rmse_stage, accuracy)

    @KratosUnittest.skip("Test test_benchmark1_1 skipped temporary.")
    def test_benchmark1_5(self):
        """
        test point load on circular tunnel with beam elements
        :return:
        """
        import math
        from analytical_solutions import calculate_max_deflections_ring, calculate_bending_moments_ring

        # Calculate
        test_name = 'test_tunnel'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        # get results
        moments = test_helper.get_moment(simulation)
        displacements = test_helper.get_displacement(simulation)

        max_x_disp = max([abs(displacement[0]) for displacement in displacements])
        max_y_disp = max([abs(displacement[1]) for displacement in displacements])
        z_moments = [moment[2] for moment in moments]
        max_moment, min_moment = max(z_moments), min(z_moments)

        # calculate analytical solution
        point_load = -100
        radius = 0.5
        youngs_modulus = 1e6
        m_inertia = 8.333e-8

        eps_h, eps_v = calculate_max_deflections_ring(point_load, radius, youngs_modulus, m_inertia)

        hor_deform_analytic = abs(eps_h * radius)
        vert_deform_analytic = abs(eps_v * radius * 2)

        max_moment_analytic = calculate_bending_moments_ring(point_load, radius, 0)
        min_moment_analytic = calculate_bending_moments_ring(point_load, radius, math.pi/2)

        # calculate error between analytical solution and numerical solution
        error_max_moment = abs(max_moment_analytic - max_moment)
        error_min_moment = abs(min_moment_analytic - min_moment)

        error_max_x_disp = abs(hor_deform_analytic - max_x_disp)
        error_max_y_disp = abs(vert_deform_analytic - max_y_disp)

        # assert if error is smaller than the precision
        precision = 0.001
        self.assertLess(error_max_moment, abs(max_moment * precision))
        self.assertLess(error_min_moment, abs(min_moment * precision))

        self.assertLess(error_max_x_disp, abs(max_x_disp * precision))
        self.assertLess(error_max_y_disp, abs(max_y_disp * precision))


if __name__ == '__main__':
    KratosUnittest.main()
