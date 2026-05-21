import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsNormalLoadHexaTests(KratosUnittest.TestCase):
    """
    This class contains tests for applying normal loads on 3D elements.
    """

    def test_hexa_8n_normal_load(self):
        test_name = 'hexa_8n_normal_load'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)

        output_file_path = os.path.join(file_path, test_name + '.post.res')
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)
        total_stress_vectors = \
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time("TOTAL_STRESS_TENSOR", 1.0,
                                                                                     output_data, element_ids=[1])[0]

        self.assertEqual(len(total_stress_vectors), 8)

        for total_stress in total_stress_vectors:
            self.assertAlmostEqual(total_stress[0], 0.0, 6)  # Sxx
            self.assertAlmostEqual(total_stress[1], 0.0, 6)  # Syy

            # Szz, note the negative sign, for a positive NORMAL_CONTACT_STRESS on the top face.
            self.assertAlmostEqual(total_stress[2], -1000.0, 6)

        top_node_nbrs = [5, 6, 7, 8]
        displacements_at_top = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 1.0, output_data,
                                                                                    node_ids=top_node_nbrs)
        bottom_node_nbrs = [1, 2, 3, 4]
        displacements_at_bottom = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 1.0, output_data,
                                                                                       node_ids=bottom_node_nbrs)

        young_modulus = 3e7
        area = 25.0
        length = 5.0
        normal_contact_stress_at_top = -1000.0
        expected_displacement = normal_contact_stress_at_top * area / (young_modulus * area / length)

        for displacement_at_top in displacements_at_top:
            self.assertAlmostEqual(displacement_at_top[2], expected_displacement, 6)
        for displacement_at_bottom in displacements_at_bottom:
            self.assertAlmostEqual(displacement_at_bottom[2], 0.0, 6)


if __name__ == '__main__':
    KratosUnittest.main()
