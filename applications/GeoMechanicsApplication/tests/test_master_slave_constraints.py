import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsMasterSlaveConstraints(KratosUnittest.TestCase):
    def check_x_displacements(self, output_data):
        end_time = 1.0
        loaded_node_ids = [6, 7]
        actual_displacements = test_helper.GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", end_time, output_data, node_ids=loaded_node_ids
        )
        expected_x_displacements = [2.0e-2, 2.0e-2]
        for actual_displacement, expected_x_displacement in zip(
            actual_displacements, expected_x_displacements
        ):
            self.assertAlmostEqual(
                actual_displacement[0], expected_x_displacement, places=3
            )

    def check_x_reaction_forces(self, output_data):
        end_time = 1.0
        supported_node_ids = [1, 4]
        actual_reaction_forces = test_helper.GiDOutputFileReader.nodal_values_at_time(
            "REACTION", end_time, output_data, node_ids=supported_node_ids
        )
        expected_x_reaction_forces = [-5000.0, -5000.0]
        for actual_reaction_force, expected_x_reaction_force in zip(
            actual_reaction_forces, expected_x_reaction_forces
        ):
            self.assertAlmostEqual(actual_reaction_force[0], expected_x_reaction_force)

    def test_equal_displacements_constraints_at_middle_section(self):
        test_files_path = test_helper.get_file_path("test_master_slave_constraints")
        test_helper.run_kratos(test_files_path)

        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_files_path, "test_master_slave_constraints.post.res")
        )
        self.check_x_displacements(output_data)
        self.check_x_reaction_forces(output_data)


if __name__ == "__main__":
    KratosUnittest.main()
