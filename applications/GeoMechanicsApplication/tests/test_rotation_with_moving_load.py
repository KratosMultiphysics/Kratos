import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

import test_helper

result_extension = ".post.res"

class KratosGeoMechanicsRotationWithMovingLoadTests(KratosUnittest.TestCase):

    def test_rotation_with_moving_load(self):
        test_name = "test_rotation_with_moving_load"
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        test_helper.run_kratos(file_path)
        res_path = os.path.join(file_path, test_name + result_extension)

        reader = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(res_path)
        rotations = simulation_output["results"]["ROTATION"]

        # Validation for Time = 0.02
        values_time_0_02 = reader.get_values_at_time(1.0, rotations)
        value_for_node = reader.get_value_at_node(1, values_time_0_02)

        # This test is a regression test, if rotation is not added to the
        # newmark upw scheme as a variable that needs to be predicted and
        # updated, the results will be different.
        self.assertAlmostEqual(-0.000177858, value_for_node[2])

    def test_rotation_with_moving_load_constant_system_matrices(self):
        """
        This test is the same as test_rotation_with_moving_load, but the system matrices are indicated as constant in
        the solver settings. Therefore, another builder and solver is used.
        """
        test_name = "test_rotation_with_moving_load_constant_system_matrices"
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        simulation = test_helper.run_kratos(file_path)

        self.assertTrue(isinstance(simulation._GetSolver().builder_and_solver, KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping))

        res_path = os.path.join(file_path, test_name + result_extension)

        reader = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(res_path)
        rotations = simulation_output["results"]["ROTATION"]

        # Validation for Time = 0.02
        values_time_0_02 = reader.get_values_at_time(1.0, rotations)
        value_for_node = reader.get_value_at_node(1, values_time_0_02)

        # This test is a regression test, if rotation is not added to the
        # newmark upw scheme as a variable that needs to be predicted and
        # updated, the results will be different.
        self.assertAlmostEqual(-0.000177858, value_for_node[2])

if __name__ == '__main__':
    KratosUnittest.main()
