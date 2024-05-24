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
        test_name    = r'truss_backbone_material\tension'
        project_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = []
        output_data.append(reader.read_output_from(output_file_name))

        stage_nr = 0
        times    = [1.0, 2.0, 3.0, 4.0]
        expected_forces = [50., 25., 50, 100]
        for time, expected_force in zip(times, expected_forces):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data[stage_nr], [1], [0])[0][0]
            self.assertAlmostEqual(float(section_force[0]), expected_force, 2)

    def test_truss_backbone_mat_compression(self):
        """
        2 element compression test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is compressed from 2 to 1 m, then elongated to 1.25 m, followed by
        compression to 0 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, reload elastically, then follow the backbone further.
        """
        test_name    = r'truss_backbone_material\compression'
        project_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = []
        output_data.append(reader.read_output_from(output_file_name))

        stage_nr = 0
        times    = [1.0, 2.0, 3.0, 4.0]
        expected_forces = [-50., -25., -50, -100]
        for time, expected_force in zip(times, expected_forces):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data[stage_nr], [1], [0])[0][0]
            self.assertAlmostEqual(float(section_force[0]), expected_force, 2)

    def test_truss_backbone_mat_tension_compression(self):
        """
        2 element elongation-compression loop test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is elongated from 2 to 3 m, then compressed to 1 m, followed by
        elongation to 4 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, meet the backbone on the other side, unload elastically
        """
        test_name    = r'truss_backbone_material\tension_compression'
        project_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = []
        output_data.append(reader.read_output_from(output_file_name))

        stage_nr = 0
        times    = [1.0, 2.0, 3.0, 4.0, 5.0]
        expected_forces = [50., -50., -100, 0., 100.]
        for time, expected_force in zip(times, expected_forces):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data[stage_nr], [1], [0])[0][0]
            self.assertAlmostEqual(float(section_force[0]), expected_force, 2)

    def test_truss_backbone_mat_compression_tension(self):
        """
        2 element compression_elongation loop test. The truss backbone material has a linear backbone with a
        stiffness of 100.0 [N/m^2] and an un- and reloading stiffness of 200.0 [N/m^2].
        The 2 element rod is compressed from 2 to 1 m, then elongated to 3 m, followed by
        compression to 0 m. The resulting force displacement curve should initially follow
        the backbone, unload elastically, meet the backbone on the other side, unload elastically
        """
        test_name    = r'truss_backbone_material\compression_tension'
        project_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(project_path)

        # output
        output_file_name = os.path.join(project_path, 'truss_backbone_mat.post.res')
        reader = test_helper.GiDOutputFileReader()
        output_data = []
        output_data.append(reader.read_output_from(output_file_name))

        stage_nr = 0
        times    = [1.0, 2.0, 3.0, 4.0, 5.0]
        expected_forces = [-50., 50., 100, 0., -100.]
        for time, expected_force in zip(times, expected_forces):
            # integration point check in element 1, integration point 1 ( uniform stress and strain so an arbitrary choice )
            section_force = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("FORCE", time, output_data[stage_nr], [1], [0])[0][0]
            self.assertAlmostEqual(float(section_force[0]), expected_force, 2)

if __name__ == '__main__':
    KratosUnittest.main()
