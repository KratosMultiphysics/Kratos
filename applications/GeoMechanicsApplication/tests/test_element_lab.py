import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsLabElementTests(KratosUnittest.TestCase):
    """
    This class contains some element tests, such as triaxial and oedometer tests
    """
    def test_triaxial_comp_6n(self):
        """
        Drained compression triaxial test on Mohr-Coulomb model with axisymmetric 2D6N elements
        It consistes of two calculation phases:
        1) apply confining stress of -100 kPa
        2) apply deviatoric stress of -200 kPa
        """
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', 'triaxial_comp_6n'))

        n_stages = 2
        run_multiple_stages.run_stages(project_path, n_stages)

        reader = GiDOutputFileReader()

        # Assert
        output_data = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage1.post.res"))
        time = 1.0
        stress_vectors_per_element = reader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data)
        number_of_elements = 2
        self.assertEqual(number_of_elements, len(stress_vectors_per_element))

        number_of_integration_points_per_element = 3
        for element_stress_vectors in stress_vectors_per_element:
            self.assertEqual(number_of_integration_points_per_element, len(element_stress_vectors))
            for stress_vector in element_stress_vectors:
                self.assertAlmostEqual(-100.0, stress_vector[0], 3)  # sigma_xx
                self.assertAlmostEqual(-100.0, stress_vector[1], 3)  # sigma_yy
                self.assertAlmostEqual(-100.0, stress_vector[2], 3)  # sigma_zz

        output_data = reader.read_output_from(os.path.join(project_path, "triaxial_comp_6n_stage2.post.res"))
        time = 1.25
        stress_vectors_per_element = reader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data)
        self.assertEqual(number_of_elements, len(stress_vectors_per_element))

        for element_stress_vectors in stress_vectors_per_element:
            self.assertEqual(number_of_integration_points_per_element, len(element_stress_vectors))
            for stress_vector in element_stress_vectors:
                self.assertAlmostEqual(-100.0, stress_vector[0], 2)  # sigma_xx
                self.assertAlmostEqual(-300.0, stress_vector[1], 2)  # sigma_yy
                self.assertAlmostEqual(-100.0, stress_vector[2], 2)  # sigma_zz


    def test_oedometer_ULFEM(self):
        """
        Oedometer test on a linear elastic model with 2D6N elements
        """
        test_name = 'oedometer_ULFEM'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses for integration_point in element]

        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(0.0,      effective_stress_xx,        3)
            self.assertAlmostEqual(-1000000, effective_stresses_yy[idx], 2)
            self.assertAlmostEqual(0.0,      effective_stresses_zz[idx], 3)

        top_node_nbrs = [1, 2]
        output_file_path = os.path.join(project_path, test_name+'.post.res')
        output_reader = GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)
        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.1, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.00990099, y_displacement, 6)

        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.7, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.0654206, y_displacement, 6)

        displacements = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 1.0, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.0909090909516868, y_displacement, 6)

    def test_oedometer_ULFEM_diff_order(self):
        """
        Oedometer test on a linear elastic model with 2D6N with different order elements
        """
        test_name = 'oedometer_ULFEM_diff_order'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))
        simulation = test_helper.run_kratos(project_path)
        effective_stresses = test_helper.get_cauchy_stress_tensor(simulation)
        displacements = test_helper.get_displacement(simulation)

        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses for integration_point in element]

        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(0.0,   effective_stress_xx,        3)
            self.assertAlmostEqual(-1000, effective_stresses_yy[idx], 3)
            self.assertAlmostEqual(0.0,   effective_stresses_zz[idx], 3)

        y_displacements = [displacement[1] for displacement in displacements]
        top_node_nbrs = [1]
        for top_node_nbr in top_node_nbrs:
            self.assertAlmostEqual(-1e-04, y_displacements[top_node_nbr], 6)

if __name__ == '__main__':
    KratosUnittest.main()
