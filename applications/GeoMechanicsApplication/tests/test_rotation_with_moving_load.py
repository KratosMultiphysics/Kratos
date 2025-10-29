import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

import test_helper

result_extension = ".post.res"

class KratosGeoMechanicsRotationWithMovingLoadTests(KratosUnittest.TestCase):

    def test_rotation_with_moving_load(self):
        test_name = "test_rotation_with_moving_load"
        file_path = test_helper.get_file_path(test_name)

        # run simulation
        test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        res_path = os.path.join(file_path, test_name + result_extension)
        simulation_output = reader.read_output_from(res_path)

        # Validation for Time = 1.0, validate total and incremental rotations also for time 1.5
        rotations_for_node_1_time_1 = reader.nodal_values_at_time("ROTATION", 1,
                                                                  simulation_output, [1])[0]

        total_rotations_for_node_1_time_1 = reader.nodal_values_at_time("TOTAL_ROTATION", 1,
                                                                        simulation_output, [1])[0]
        total_rotations_for_node_1_time_1_5 = reader.nodal_values_at_time("TOTAL_ROTATION", 1.5,
                                                                          simulation_output, [1])[0]

        incremental_rotations_for_node_1_time_1 = reader.nodal_values_at_time("INCREMENTAL_ROTATION", 1,
                                                                              simulation_output, [1])[0]
        incremental_rotations_for_node_1_time_1_5 = reader.nodal_values_at_time("INCREMENTAL_ROTATION", 1.5,
                                                                                simulation_output, [1])[0]

        # This test is a regression test, if rotation is not added to the
        # newmark upw scheme as a variable that needs to be predicted and
        # updated, the results will be different. TOTAL_ROTATION should be equal to ROTATION
        self.assertAlmostEqual(-0.000177858, rotations_for_node_1_time_1[2])

        self.assertAlmostEqual(-0.000177858, total_rotations_for_node_1_time_1[2])
        self.assertAlmostEqual(-0.000238568, total_rotations_for_node_1_time_1_5[2])

        # Incremental rotations are different from total rotation after the first time step
        self.assertAlmostEqual(-0.000177858, incremental_rotations_for_node_1_time_1[2])
        self.assertAlmostEqual(-6.07101e-05, incremental_rotations_for_node_1_time_1_5[2])

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

        reader = GiDOutputFileReader()
        res_path = os.path.join(file_path, test_name + result_extension)
        simulation_output = reader.read_output_from(res_path)

        rotations_for_node_1_time_1 = reader.nodal_values_at_time("ROTATION", 1,
                                                           simulation_output, [1])[0]

        # This test is a regression test, if rotation is not added to the
        # ResidualBasedBlockBuilderAndSolverWithMassAndDamping as a variable that
        # needs to update the first/second time derivatives, the results will be different.
        self.assertAlmostEqual(-0.000177858, rotations_for_node_1_time_1[2])

if __name__ == '__main__':
    KratosUnittest.main()
