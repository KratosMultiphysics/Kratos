import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsTrussBackboneMaterialTests(KratosUnittest.TestCase):
    """
    This class contains tests for the truss backbone material
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_truss_backbone_mat_tension(self):
        """
        2 element elongation test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is elongated from 2 to 3 m, then compressed to 2.75 m, followed by
        elongation to 4 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, reload elastically, then follow the backbone further.
        """
        project_path = test_helper.get_file_path(os.path.join("truss_backbone_material", "tension"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [50., 25., 50, 100]
        expected_displacements_x = [1.0, 0.75, 1.0, 2.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 2)
            displacement = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_truss_backbone_mat_compression(self):
        """
        2 element compression test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is compressed from 2 to 1 m, then elongated to 1.25 m, followed by
        compression to 0 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, reload elastically, then follow the backbone further.
        """
        project_path = test_helper.get_file_path(os.path.join("truss_backbone_material", "compression"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [-50., -25., -50, -100]
        expected_displacements_x = [-1.0, -0.75, -1.0, -2.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 2)
            displacement = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_truss_backbone_mat_tension_compression(self):
        """
        2 element elongation-compression loop test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is elongated from 2 to 3 m, then compressed to 1 m, followed by
        elongation to 4 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, meet the backbone on the other side, unload elastically
        """
        project_path = test_helper.get_file_path(os.path.join("truss_backbone_material", "tension_compression"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        expected_forces_x = [50., -50., -100, 0., 100.]
        expected_displacements_x = [1.0, 0.0, -1.0, 0.0, 1.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 2)
            displacement = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

    def test_truss_backbone_mat_compression_tension(self):
        """
        2 element compression_elongation loop test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is compressed from 2 to 1 m, then elongated to 3 m, followed by
        compression to 0 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, meet the backbone on the other side, unload elastically
        """
        project_path = test_helper.get_file_path(os.path.join("truss_backbone_material", "compression_tension"))
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_name)

        times = [1.0, 2.0, 3.0, 4.0, 5.0]
        expected_forces_x = [-50., 50., 100, 0., -100.]
        expected_displacements_x = [-1.0, 0.0, 1.0, 0.0, -1.0]
        for time, expected_force_x, expected_displacement_x in zip(times, expected_forces_x, expected_displacements_x):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data, [1], [0])[0][0]
            self.assertAlmostEqual(section_force[0], expected_force_x, 2)
            displacement = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data, [3])[0]
            self.assertAlmostEqual(displacement[0], expected_displacement_x, 2)

if __name__ == '__main__':
    KratosUnittest.main()
