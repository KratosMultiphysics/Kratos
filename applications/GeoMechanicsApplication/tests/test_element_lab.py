import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsLabElementTests(KratosUnittest.TestCase):
    """
    This class contains some element tests, such as triaxial and oedometer tests
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_triaxial_comp_6n(self):
        """
        Drained compression triaxial test on Mohr-Coulomb model with axisymmetric 2D6N elements
        It consistes of two calculation phases:
        1) apply confining stress of -100 kPa
        2) apply deviatoric stress of -200 kPa
        """
        test_name = 'triaxial_comp_6n'
        project_path = test_helper.get_file_path(os.path.join('test_element_lab', test_name))

        n_stages = 2
        directory_names = test_helper.get_separated_directory_names(project_path, n_stages)

        stages = test_helper.get_separated_stages(directory_names)

        effective_stresses_stages = [None] * n_stages
        # run stages and get results
        for idx, stage in enumerate(stages):
            os.chdir(directory_names[idx])
            stage.Run()
            effective_stresses_stages[idx] = test_helper.get_cauchy_stress_tensor(stage)

        # Assert
        stage_nr = 0
        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(-100.0, effective_stress_xx,        3)
            self.assertAlmostEqual(-100.0, effective_stresses_yy[idx], 3)
            self.assertAlmostEqual(-100.0, effective_stresses_zz[idx], 3)

        stage_nr = 1
        effective_stresses_xx = [integration_point[0,0] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        effective_stresses_yy = [integration_point[1,1] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        effective_stresses_zz = [integration_point[2,2] for element in effective_stresses_stages[stage_nr] for integration_point in element]
        # Assert integration point information
        for idx, effective_stress_xx in enumerate(effective_stresses_xx):
            self.assertAlmostEqual(-100.0, effective_stress_xx,        2)
            self.assertAlmostEqual(-300.0, effective_stresses_yy[idx], 2)
            self.assertAlmostEqual(-100.0, effective_stresses_zz[idx], 2)

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
        output_reader = test_helper.GiDOutputFileReader()
        output_data = output_reader.read_output_from(output_file_path)
        displacements = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.1, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.00990099, y_displacement, 6)

        displacements = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 0.7, output_data, node_ids=top_node_nbrs)
        for displacement in displacements:
            y_displacement = displacement[1]
            self.assertAlmostEqual(-0.0654206, y_displacement, 6)

        displacements = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", 1.0, output_data, node_ids=top_node_nbrs)
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
