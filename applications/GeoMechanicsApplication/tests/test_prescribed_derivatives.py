import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

import test_helper


class KratosGeoMechanicsPrescribedDerivatives(KratosUnittest.TestCase):
    """
    This class contains tests which check if prescribed derivatives are correctly applied to the model and not overwritten
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_prescribed_acceleration(self):
        """
        This test checks if the prescribed acceleration in the x-direction is correctly applied to the model
        """

        test_name = 'prescribed_acceleration'
        file_path = test_helper.get_file_path(os.path.join('prescribed_derivative_tests', test_name))
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, test_name + '.post.res')
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        times = [0.05, 0.1, 0.15, 0.2]
        expected_x_accelerations = [-0.01, -0.005, -0.003, 0.01]
        for expected_x_acceleration, time in zip(expected_x_accelerations, times):
            accelerations = GiDOutputFileReader.nodal_values_at_time(
                "ACCELERATION", time, output_data, node_ids=[1, 2, 3]
            )
            self.assertEqual(len(accelerations), 3)

            for acceleration in accelerations:
                self.assertAlmostEqual(expected_x_acceleration, acceleration[0], 6)
                self.assertAlmostEqual(0.0, acceleration[1], 6)
                self.assertAlmostEqual(0.0, acceleration[2], 6)

    def test_prescribed_velocity(self):
        """
        This test checks if the prescribed velocity in the x-direction is correctly applied to the model
        """

        test_name = 'prescribed_velocity'
        file_path = test_helper.get_file_path(os.path.join('prescribed_derivative_tests', test_name))
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, test_name + '.post.res')
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        times = [0.05, 0.1, 0.15, 0.2]
        expected_velocities = [-0.01, -0.005, -0.003, 0.01]
        for expected_velocity, time in zip(expected_velocities, times):
            velocities = GiDOutputFileReader.nodal_values_at_time(
                "VELOCITY", time, output_data, node_ids=[1, 2, 3]
            )

            self.assertEqual(len(velocities), 3)
            for velocity in velocities:
                self.assertAlmostEqual(expected_velocity, velocity[0], 6)
                self.assertAlmostEqual(0.0, velocity[1], 6)
                self.assertAlmostEqual(0.0, velocity[2], 6)

    def test_prescribed_dt_water_pressure(self):
        """
        This test checks if the prescribed dt_water_pressure is correctly applied to the model
        """

        test_name = 'prescribed_dt_water_pressure'
        file_path = test_helper.get_file_path(os.path.join('prescribed_derivative_tests', test_name))
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, test_name + '.post.res')
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)

        times = [0.05, 0.1, 0.15, 0.2]
        expected_dt_water_pressures = [-0.01, -0.005, -0.003, 0.01]
        for expected_dt_water_pressure, time in zip(expected_dt_water_pressures, times):
            dt_water_pressures = GiDOutputFileReader.nodal_values_at_time(
                "DT_WATER_PRESSURE", time, output_data, node_ids=[1, 2, 3]
            )
            self.assertEqual(len(dt_water_pressures), 3)
            for dt_water_pressure in dt_water_pressures:
                self.assertAlmostEqual(expected_dt_water_pressure, dt_water_pressure, 6)


if __name__ == '__main__':
    KratosUnittest.main()
