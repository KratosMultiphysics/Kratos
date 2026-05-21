import sys
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsElementTypeTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_triangle_3n(self):
        test_name = 'test_triangle_3n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 1, 5]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

        bottom_node_ids = [5, 7, 9]
        self.assertVerticalStressAtBottomNodes(output_data, bottom_node_ids)

    def test_triangle_3n_rebuild_level_0(self):
        test_name = 'test_triangle_3n_rebuild_0'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, "test_triangle_3n.post.res"))
        top_node_nbrs = [0, 1, 5]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

        bottom_node_ids = [5, 7, 9]
        self.assertVerticalStressAtBottomNodes(output_data, bottom_node_ids)

    def test_triangle_6n(self):
        test_name = 'test_triangle_6n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 2, 4, 9, 15]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

    def test_triangle_10n(self):
        test_name = 'test_triangle_10n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [9, 8, 7, 6]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

    def test_triangle_15n(self):
        test_name = 'test_triangle_15n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [12, 11, 10, 9, 8]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

    def test_triangle_3n_fic(self):
        test_name = 'test_triangle_3n_fic'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 1, 5]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

        bottom_node_ids = [5, 7, 9]
        self.assertVerticalStressAtBottomNodes(output_data, bottom_node_ids)

    def test_triangle_6n_fic(self):
        test_name = 'test_triangle_6n_fic'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 2, 4, 9, 15]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

    def test_quad_4n(self):
        test_name = 'test_quad_4n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 1, 5]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

        bottom_node_ids = [5, 7, 9]
        self.assertVerticalStressAtBottomNodes(output_data, bottom_node_ids)

    def test_quad_8n(self):
        test_name = 'test_quad_8n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 2, 4, 8, 12]
        n_dim = 2
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

    def test_tetra_4n(self):
        test_name = 'test_tetra_4n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 2, 9, 3, 6, 13, 8, 14, 20]
        n_dim = 3
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)

        bottom_node_ids = [11, 12, 17, 20, 22, 23, 25, 26, 27]
        self.assertVerticalStressAtBottomNodes(output_data, bottom_node_ids)

    def test_tetra_10n(self):
        test_name = 'test_tetra_10n'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        output_data = self.fetchOutputFromFile(os.path.join(file_path, f"{test_name}.post.res"))
        top_node_nbrs = [0, 3, 9, 28, 53,
                         1, 6, 13, 32,
                         10, 15, 21, 44, 76,
                         27, 34, 41, 64, 94,
                         51, 56, 74, 95, 110]
        n_dim = 3
        self.assert_linear_elastic_block(simulation, output_data, top_node_nbrs, n_dim)


    def assert_linear_elastic_block(self, simulation, output_data, top_node_nbrs, n_dim):
        """
        Assert results of a linear elastic block. The sides of the block can move freely in vertical direction and are
        fixed in horizontal direction. The bottom of the block is fixed. On top of the block, a load of 10kN/m2 is
        placed. Results are: total stresses, effective stresses, displacements and green langrange strains.

        :param simulation: Kratos simulation
        :param top_node_nbrs: node numbers of the nodes at the top of the geometry
        :param n_dim: number of dimensions
        :return:
        """
        total_stresses = test_helper.get_total_stress_tensor(simulation)
        total_stresses_xx = [integration_point[0,0] for element in total_stresses for integration_point in element]
        if n_dim >= 2:
            total_stresses_yy = [integration_point[1,1] for element in total_stresses for integration_point in element]
        if n_dim >= 3:
            total_stresses_zz = [integration_point[2,2] for element in total_stresses for integration_point in element]

        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses for integration_point in element]
        if n_dim >= 2:
            effective_stresses_yy = [integration_point[1,1] for element in effective_stresses for integration_point in element]
        if n_dim >= 3:
            effective_stresses_zz = [integration_point[2,2] for element in effective_stresses for integration_point in element]

        end_time = 1.0
        displacements = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", end_time, output_data)
        x_displacements = [displacement[0] for displacement in displacements]
        if n_dim >= 2:
            y_displacements = [displacement[1] for displacement in displacements]
        if n_dim >= 3:
            z_displacements = [displacement[2] for displacement in displacements]

        green_lagrange_strains = test_helper.get_green_lagrange_strain_tensor(simulation)
        green_lagrange_strains_xx = [integration_point[0,0] for element in green_lagrange_strains for integration_point in
                                     element]
        if n_dim == 2:
            green_lagrange_strains_yy = [integration_point[1,1] for element in green_lagrange_strains for integration_point in
                                         element]
        elif n_dim == 3:
            green_lagrange_strains_yy = [integration_point[1,1] for element in green_lagrange_strains for integration_point in
                                         element]
            green_lagrange_strains_zz = [integration_point[2,2] for element in green_lagrange_strains for integration_point in
                                         element]

        # Assert integration point information
        for idx, total_stress_xx in enumerate(total_stresses_xx):
            self.assertAlmostEqual(0.0, total_stress_xx)
            self.assertAlmostEqual(-1e4, total_stresses_yy[idx])
            if n_dim >= 3:
                self.assertAlmostEqual(0.0, total_stresses_zz[idx])

            self.assertAlmostEqual(0.0, effective_stresses_xx[idx])
            self.assertAlmostEqual(-1e4, effective_stresses_yy[idx])
            if n_dim >= 3:
                self.assertAlmostEqual(0.0, effective_stresses_zz[idx])

            self.assertAlmostEqual(0.0, green_lagrange_strains_xx[idx])
            self.assertAlmostEqual(-0.00033333, green_lagrange_strains_yy[idx])
            if n_dim >= 3:
                self.assertAlmostEqual(0.0, green_lagrange_strains_zz[idx])

        # Assert displacements
        for x_displacement in x_displacements:
            self.assertAlmostEqual(0.0, x_displacement)

        for top_node_nbr in top_node_nbrs:
            self.assertAlmostEqual(-0.00033333, y_displacements[top_node_nbr], 6)

        if n_dim >= 3:
            for z_displacement in z_displacements:
                self.assertAlmostEqual(0.0, z_displacement)

    def fetchOutputFromFile(self, output_file_path):
        reader = test_helper.GiDOutputFileReader()
        return reader.read_output_from(output_file_path)

    def assertVerticalStressAtBottomNodes(self, output_data, bottom_node_ids):
        end_time = 1.0
        nodal_stress_tensors = test_helper.GiDOutputFileReader.nodal_values_at_time("NODAL_CAUCHY_STRESS_TENSOR", end_time, output_data, bottom_node_ids)
        expected_stress_yy = -1e4
        for stress_tensor in nodal_stress_tensors:
            self.assertAlmostEqual(stress_tensor[1], expected_stress_yy)

if __name__ == '__main__':
    KratosUnittest.main()
