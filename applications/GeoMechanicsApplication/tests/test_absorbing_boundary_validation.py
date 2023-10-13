import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsAbsorbingBoundaryColumnValidationTests(KratosUnittest.TestCase):
    """
    This class contains tests which check if the lysmer absorbing boundary works in a 1d column made out of different
    element types
    """

    def setUp(self):
        self.E = 10000
        self.nu = 0.2
        self.rho = 2.65 * 0.7
        self.load = -10
        self.height_column = 10
        self.vp = None
        self.expected_velocity = None

        self.calculate_expected_p_wave_velocity()
        self.calculate_expected_velocity_column()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_absorbing_boundary_on_1d_column_tetra(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of tetrahedrals. The boundary is a 3d3n triangle.

        :return:
        """
        test_name = 'test_lysmer_boundary_column3d_tetra.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        # quarter node, middle node, three quarter node
        node_nbrs_to_assert = [33,54,81]
        direction = 2

        self.run_and_assert_1d_column(file_path, node_nbrs_to_assert, direction)

    def test_absorbing_boundary_on_1d_column_tetra_horizontal(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of tetrahedrals. In this case the wave columns
        moves in a positive y-direction. The boundary is a 3d3n triangle.
        :return:
        """

        test_name = 'test_lysmer_boundary_column3d_tetra_hor.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        # quarter node, middle node, three quarter node
        node_nbrs_to_assert = [33,54,81]
        direction = 1

        self.run_and_assert_1d_column(file_path, node_nbrs_to_assert, direction)

    def test_absorbing_boundary_on_1d_column_hexa(self):
        """
        Tests the lysmer absorbing boundary condition on a column made of hexahedrals. The boundary is a 3d4n rectangle.

        :return:
        """
        test_name = 'test_lysmer_column3d_hexa.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        # quarter node, middle node, three quarter node
        node_nbrs_to_assert = [31, 54, 81]
        direction = 2

        self.run_and_assert_1d_column(file_path, node_nbrs_to_assert, direction)

    def run_and_assert_1d_column(self, file_path, node_nbrs, direction):
        """
        Runs and asserts a dynamic test on a 1d column. This test checks when a p-wave arrives at a certain coordinate
        and the velocity of the corresponding node afterwards

        :param file_path: path of test
        :param node_nbrs: node nbrs to be checked
        :param direction: direction of the wave 0 = x; 1 = y; 2 = z;
        :return:

        """

        # get name of output file
        _, output_file_name = os.path.split(file_path)
        output_file_name = os.path.splitext(output_file_name)[0] + ".post.res"
        output_file_path = os.path.join(file_path,output_file_name)

        # clear old results
        if os.path.exists(output_file_path):
            os.remove(output_file_path)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        # get results from calculation
        coords = test_helper.get_nodal_coordinates(simulation)
        res = test_helper.get_nodal_variable_from_ascii(output_file_path, "DISPLACEMENT")

        # calculate and check velocity and time of arrival at all nodes
        for node_nbr in node_nbrs:
            vert_coord = coords[node_nbr-1][direction]

            # distance between top and vertical coordinate of node
            dist = self.height_column-vert_coord

            # expected time of wave arrival
            expected_ini_time = dist/self.vp

            # find index of expected time of wave arrival in time list
            ini_time_idx = test_helper.find_closest_index_greater_than_value(res["time"], expected_ini_time)

            # calculate velocity after wave arrival
            dt = (res["time"][-1] - res["time"][ini_time_idx])
            velocity_part_two = (res[str(node_nbr)][-1][direction] - res[str(node_nbr)][ini_time_idx][direction])/dt

            # assert velocities
            self.assertAlmostEqual(velocity_part_two, self.expected_velocity, 2)

    def calculate_expected_p_wave_velocity(self):
        """
        Calculates the expected p_wave velocity in a 1d column
        Returns
        -------

        """
        # calculate expected p-wave velocity

        Ec = self.E * (1 - self.nu) / ((1 + self.nu) * (1 - 2 * self.nu))

        # calculate p-wave velocity
        self.vp = (Ec / self.rho) ** 0.5

    def calculate_expected_velocity_column(self):
        """
        Calculates the expected velocity in a 1d column based on the p-wave
        Returns
        -------

        """
        self.expected_velocity = self.load / (self.vp * self.rho)


if __name__ == '__main__':
    KratosUnittest.main()
