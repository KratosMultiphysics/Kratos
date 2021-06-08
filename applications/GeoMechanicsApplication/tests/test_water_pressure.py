import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsWaterPressureTests(KratosUnittest.TestCase):
    """
    This class contains tests which compare water pressure calculated in kratos with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_inclined_phreatic_line(self):
        """
        test hydrostatic water pressure under an inclided phreatic line ranging over the width of the geometry

        :return:
        """
        test_name = 'test_inclinded_phreatic_line'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        water_pressure = test_helper.get_water_pressure(simulation)

        p_bot_left = water_pressure[0]
        p_bot_middle = water_pressure[47]
        p_bot_right = water_pressure[185]

        p_top_left = water_pressure[186]
        p_top_middle = water_pressure[222]
        p_top_right = water_pressure[250]

        self.assertAlmostEqual(-10000, p_bot_left)
        self.assertAlmostEqual(-7500, p_bot_middle)
        self.assertAlmostEqual(-5000, p_bot_right)

        self.assertAlmostEqual(0, p_top_left)
        self.assertAlmostEqual(0, p_top_middle)
        self.assertAlmostEqual(0, p_top_right)

    def test_inclinded_phreatic_line_smaller_line(self):
        """
        test hydrostatic water pressure under an inclided phreatic line ranging over a part of the geometry

        :return:
        """
        test_name = 'test_inclinded_phreatic_line_smaller_line'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        # Get water pressure in all the nodes
        water_pressure = test_helper.get_water_pressure(simulation)

        # get node information of the bottom nodes in the geometry
        bottom_n_nbrs = [0, 2, 4, 9, 16, 24, 34, 47, 61, 78, 96, 115, 136, 162, 185]
        x_coords_boundary_nodes = [i / 14 for i in range(15)]

        # Set turning points of the phreatic line
        head_values = [1, 1, 0.5, 0.5]
        head_x_coords = [0, 0.3, 0.7, 1]

        # Analytically calculate water pressure at the bottom of the geometry
        dhdx = (head_values[2] - head_values[1]) / (head_x_coords[2] - head_x_coords[1])
        gamma_w = 10000

        p_bottom_analytic = [-head_values[1] * gamma_w if x <= head_x_coords[1]
                             else -head_values[2] * gamma_w if x >= head_x_coords[2]
        else -(dhdx * (x - head_x_coords[1]) + head_values[1]) * gamma_w
                             for x in x_coords_boundary_nodes]

        # assert
        for idx, node_nr in enumerate(bottom_n_nbrs):
            self.assertAlmostEqual(p_bottom_analytic[idx], water_pressure[node_nr], delta=1e-6)

    @KratosUnittest.skip("unit test not implemented")
    def test_interpolate_water_pressure(self):
        pass

    def test_interpolate_water_pressure_inclined(self):
        """
        Test interpolate water pressure with an inclined phreatic line.
        The geometry consists out of 3 blocks. On the top block, an inclined phreatic line is defined. On the bottom
        block, 100 pressure is defined. The water pressure in the middle block is interpolated between the water pressure
        in the top and bottom block
        :return:
        """

        test_name = 'interpolate_line_2'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        water_pressure = test_helper.get_water_pressure(simulation)

        # Node numbers of the left and right boundary of the middle block
        middle_left_n_nbrs = [181, 155, 134, 110, 92, 75, 59, 45]
        middle_right_n_nbrs = [301, 287, 274, 260, 249, 235, 227, 217]

        # Pore pressures at the corners of the middle block
        p_top_middle_left = -10000
        p_bot_middle_left = -100

        p_top_middle_right = -5000
        p_bot_middle_tight = -100

        # Vertical coordinates of the left and right boundary nodes of the middle block
        vert_coord_bound_mid = [-i / 7 * 0.5 + 1.0 for i in range(8)]

        # Analytically calculate pore pressure along left and right boundary of the middle block
        dpdy_left = (p_top_middle_left - p_bot_middle_left) / (vert_coord_bound_mid[0] - vert_coord_bound_mid[-1])
        dpdy_right = (p_top_middle_right - p_bot_middle_tight) / (vert_coord_bound_mid[0] - vert_coord_bound_mid[-1])

        p_middle_left_analytic = [(dpdy_left * (y - vert_coord_bound_mid[0]) + p_top_middle_left)
                                  for y in vert_coord_bound_mid]

        p_middle_right_analytic = [(dpdy_right * (y - vert_coord_bound_mid[0]) + p_top_middle_right)
                                   for y in vert_coord_bound_mid]

        # Assert
        for idx, node_nr in enumerate(middle_left_n_nbrs):
            self.assertAlmostEqual(p_middle_left_analytic[idx], water_pressure[node_nr], delta=1e-6)

        for idx, node_nr in enumerate(middle_right_n_nbrs):
            self.assertAlmostEqual(p_middle_right_analytic[idx], water_pressure[node_nr], delta=1e-6)

    @KratosUnittest.skip("check if test is correct")
    def test_normal_load_triangle_3n_fic(self):
        test_name = 'test_normal_load_triangle_3n_fic'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        n_dim = 2
        self.assert_linear_elastic_saturated_block(simulation, n_dim)

    def test_normal_load_triangle_6n(self):
        test_name = 'test_normal_load_triangle_6n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        n_dim = 2
        self.assert_linear_elastic_saturated_block(simulation, n_dim)

    def test_normal_load_quad_8n(self):
        test_name = 'test_normal_load_quad_8n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        n_dim = 2
        self.assert_linear_elastic_saturated_block(simulation, n_dim)

    def test_normal_load_tetra_10n(self):
        test_name = 'test_normal_load_tetra_10n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        n_dim = 3
        self.assert_linear_elastic_saturated_block(simulation, n_dim)

    def assert_linear_elastic_saturated_block(self, simulation, n_dim):
        """
        Assert results of a linear elastic fully saturated block. The sides of the block can move freely in vertical
        direction and are fixed in horizontal direction. The bottom of the block is fixed. On top of the block, a
        hydrostatic water load is applied with a hydraulic head 1 meter above the top of the soil placed. Results are:
        total stresses, effective stresses, displacements, green langrange strains and water pressure.

        :param simulation: Kratos simulation
        :param top_node_nbrs: node numbers of the nodes at the top of the geometry
        :param n_dim: number of dimensions
        :return:
        """
        # get total stresses
        total_stresses = test_helper.get_total_stress_tensor(simulation)
        total_stresses_xx = [integration_point[0,0] for element in total_stresses for integration_point in element]

        if n_dim >= 2:
            total_stresses_yy = [integration_point[1,1] for element in total_stresses for integration_point in element]
        if n_dim >= 3:
            total_stresses_zz = [integration_point[2,2] for element in total_stresses for integration_point in element]

        # get effective stresses
        efective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        efective_stresses_xx = [integration_point[0,0] for element in efective_stresses for integration_point in element]
        if n_dim >= 2:
            efective_stresses_yy = [integration_point[1,1] for element in efective_stresses
                                    for integration_point in element]
        if n_dim >= 3:
            efective_stresses_zz = [integration_point[2,2] for element in efective_stresses
                                    for integration_point in element]

        # get strains
        green_lagrange_strains = test_helper.get_green_lagrange_strain_tensor(simulation)
        green_lagrange_strains_xx = [integration_point[0,0] for element in green_lagrange_strains for integration_point in
                                     element]
        if n_dim >= 2:
            green_lagrange_strains_yy = [integration_point[1,1] for element in green_lagrange_strains
                                         for integration_point in element]
        if n_dim >= 3:
            green_lagrange_strains_zz = [integration_point[2,2] for element in green_lagrange_strains
                                         for integration_point in element]

        # get displacements
        displacements = test_helper.get_displacement(simulation)
        x_displacements = [displacement[0] for displacement in displacements]
        if n_dim >= 2:
            y_displacements = [displacement[1] for displacement in displacements]
        if n_dim >= 3:
            z_displacements = [displacement[2] for displacement in displacements]

        # get water pressures
        water_pressures = test_helper.get_water_pressure(simulation)

        # get coordinates
        nodal_coordinates = test_helper.get_nodal_coordinates(simulation)
        gauss_coordinates = test_helper.get_gauss_coordinates(simulation)

        # get the coordinates of all gauss points in a list
        gauss_coordinates_x = [integration_point[0] for element in gauss_coordinates for integration_point in
                                         element]
        if n_dim >= 2:
            gauss_coordinates_y = [integration_point[1] for element in gauss_coordinates for integration_point in
                                             element]
        if n_dim >= 3:
            gauss_coordinates_z = [integration_point[2] for element in gauss_coordinates for integration_point in
                                             element]

        # Calculate expected values
        expected_water_pressure = [coord[1] * 1e4 - 2e4 for coord in nodal_coordinates]
        #todo correct this
        expected_displacements_y = [coord[1] * 0.0001667 for coord in nodal_coordinates]

        expected_total_stress_xx = [gauss_coordinate_y * 1e4 - 2e4 for gauss_coordinate_y in gauss_coordinates_y]
        if n_dim >= 2:
            expected_eff_stress_yy = [(1 - gauss_coordinate_y) * 1e4 for gauss_coordinate_y in gauss_coordinates_y]
            expected_strain_yy = [(1-gauss_coordinate_y) * 0.00033333 for gauss_coordinate_y in gauss_coordinates_y]
        if n_dim >= 3:
            expected_total_stress_zz = expected_total_stress_xx

        # Assert integration point information
        for idx, total_stress_xx in enumerate(total_stresses_xx):
            self.assertAlmostEqual(expected_total_stress_xx[idx], total_stress_xx, 4)
            if n_dim >= 2:
                self.assertAlmostEqual(-1e4, total_stresses_yy[idx], 1)
            if n_dim >= 3:
                self.assertAlmostEqual(expected_total_stress_zz[idx], total_stresses_zz[idx], 4)

            self.assertAlmostEqual(0.0, efective_stresses_xx[idx], 4)
            if n_dim >= 2:
                self.assertAlmostEqual(expected_eff_stress_yy[idx], efective_stresses_yy[idx],4)
            if n_dim >= 3:
                self.assertAlmostEqual(0.0, efective_stresses_zz[idx], 4)

            self.assertAlmostEqual(0.0, green_lagrange_strains_xx[idx])
            if n_dim >= 2:
                self.assertAlmostEqual(expected_strain_yy[idx], green_lagrange_strains_yy[idx])
            if n_dim >= 3:
                self.assertAlmostEqual(0.0, green_lagrange_strains_zz[idx])

        # Assert displacements
        for x_displacement in x_displacements:
            self.assertAlmostEqual(0.0, x_displacement)

        for node_idx in range(len(nodal_coordinates)):
            #todo correct expected displacement
            # if n_dim >= 2:
            #     self.assertAlmostEqual(expected_displacements_y[node_idx], y_displacements[node_idx], 6)
            self.assertAlmostEqual(expected_water_pressure[node_idx], water_pressures[node_idx], 6)

        if n_dim >= 3:
            for z_displacement in z_displacements:
                self.assertAlmostEqual(0.0, z_displacement)

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsWaterPressureTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
