import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsMasterSlaveConstraints(KratosUnittest.TestCase):
    def check_x_displacements(self, output_data):
        end_time = 1.0
        loaded_node_ids = [6, 7]
        actual_x_displacements = [
            displacement_vector[0]
            for displacement_vector in GiDOutputFileReader.nodal_values_at_time(
                "DISPLACEMENT", end_time, output_data, node_ids=loaded_node_ids
            )
        ]
        expected_x_displacements = [2.0e-2, 2.0e-2]
        self.assertTrue(
            test_helper.are_iterables_almost_equal(
                expected_x_displacements, actual_x_displacements
            ),
            msg=f"x displacements don't compare equal: expected = {expected_x_displacements}, actual = {actual_x_displacements}",
        )

    def check_displacements_at_tied_nodes(self, output_data):
        end_time = 1.0
        master_node_ids = [2, 3]
        master_node_displacements = (
            GiDOutputFileReader.nodal_values_at_time(
                "DISPLACEMENT", end_time, output_data, node_ids=master_node_ids
            )
        )
        slave_node_ids = [5, 8]
        slave_node_displacements = GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", end_time, output_data, node_ids=slave_node_ids
        )
        self.assertTrue(
            test_helper.are_iterables_almost_equal(
                master_node_displacements, slave_node_displacements
            ),
            msg=f"displacement vectors don't compare equal: master node displacements = {master_node_displacements}, slave node displacements = {slave_node_displacements}",
        )

    def check_x_reaction_forces(self, output_data):
        end_time = 1.0
        supported_node_ids = [1, 4]
        actual_x_reaction_forces = [
            reaction_force_vector[0]
            for reaction_force_vector in GiDOutputFileReader.nodal_values_at_time(
                "REACTION", end_time, output_data, node_ids=supported_node_ids
            )
        ]
        expected_x_reaction_forces = [-5000.0, -5000.0]
        self.assertTrue(
            test_helper.are_iterables_almost_equal(
                expected_x_reaction_forces, actual_x_reaction_forces
            ),
            msg=f"x reaction forces don't compare equal: expected = {expected_x_reaction_forces}, actual = {actual_x_reaction_forces}",
        )

    def test_equal_displacements_constraints_at_middle_section(self):
        test_files_path = test_helper.get_file_path("test_master_slave_constraints")
        test_helper.run_kratos(test_files_path)

        reader = GiDOutputFileReader()
        output_data = reader.read_output_from(
            os.path.join(test_files_path, "test_master_slave_constraints.post.res")
        )
        self.check_x_displacements(output_data)
        self.check_displacements_at_tied_nodes(output_data)
        self.check_x_reaction_forces(output_data)


if __name__ == "__main__":
    KratosUnittest.main()
