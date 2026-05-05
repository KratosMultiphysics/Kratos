import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsPiecewiseLinearMomentCapacityTests(KratosUnittest.TestCase):
    """
    This class contains tests for the piecewise linear moment capacity material
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_piecewise_linear_moment_capacity_tension(self):
        """
        2 element elongation test for piecewise linear moment capacity material.
        """
        project_path = test_helper.get_file_path(os.path.join("piecewise_linear_moment_capacity", "tension"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'piecewise_linear_moment_capacity_mat.post.res')
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [0.81055, -1.08659, 1.20015, 1.5]
        expected_displacements_x = [1.0, 0.75, 1.0, 2.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            section_force = GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 4)
            displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_piecewise_linear_moment_capacity_compression(self):
        """
        2 element compression test for piecewise linear moment capacity material.
        """
        project_path = test_helper.get_file_path(os.path.join("piecewise_linear_moment_capacity", "compression"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'piecewise_linear_moment_capacity_mat.post.res')
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [-0.81055, 1.08659, -1.20015, -1.5]
        expected_displacements_x = [-1.0, -0.75, -1.0, -2.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            section_force = GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 4)
            displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_piecewise_linear_moment_capacity_tension_compression(self):
        """
        Elongation-compression loop test for piecewise linear moment capacity material.
        """
        project_path = test_helper.get_file_path(os.path.join("piecewise_linear_moment_capacity", "tension_compression"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'piecewise_linear_moment_capacity_mat.post.res')
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [0.81055, -1.38546, -1.5, -1.5]
        expected_displacements_x = [1.0, 0.0, -1.0, 0.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            section_force = GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 4)
            displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_piecewise_linear_moment_capacity_compression_tension(self):
        """
        Compression-elongation loop test for piecewise linear moment capacity material.
        """
        project_path = test_helper.get_file_path(os.path.join("piecewise_linear_moment_capacity", "compression_tension"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'piecewise_linear_moment_capacity_mat.post.res')
        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [-0.81055, 1.38546, 1.5, 1.5]
        expected_displacements_x = [-1.0, -0.0, 1.0, 0.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            section_force = GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 4)
            displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)


if __name__ == '__main__':
    KratosUnittest.main()
