import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsDirichletUTests(KratosUnittest.TestCase):
    """
    This class contains a test for displacement Dirichlet boundary condition
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_dirichlet_u(self):
        """
        4 element elongation test in 2 stages. The incremental elastic material should show 
        linear increase in stress_yy over the steps, with continuity over the stages.
        To reach this a saw tooth like table for the prescribed (stage) displacement is given.
        H = 1 [m]
        Young's modulus E = 30E+06 [N/m^2]
        stage 1) elongate 0.05 [m] in 2 steps
        stage 2) elongate another 0.05 [0.05] in 2 steps
        """
        test_name    = 'dirichlet_u'
        project_path = test_helper.get_file_path(test_name)
        n_stages     = 2
        stages       = test_helper.get_stages(project_path,n_stages)

        # name of output file
        output_file_names = [os.path.join(project_path, f'dirichlet_u_stage{i+1}.post.res') for i in
                            range(n_stages)]
        output_data       = []

        # run stages and get results
        for stage, output_file_name in zip(stages, output_file_names):
            stage.Run()
            reader = test_helper.GiDOutputFileReader()
            output_data.append(reader.read_output_from(output_file_name))

        E = 30E+06
        stage_nr = 0
        for time in [0.5, 1.0]:
            # displacement of the central node ( half the prescribed displacement at the top )
            displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            displacement_8_y = displacements_8[1]
            self.assertAlmostEqual(time*0.05/2., displacement_8_y, 2)
            total_displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("TOTAL_DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            total_displacement_8_y = total_displacements_8[1]
            self.assertAlmostEqual(time*0.05/2, total_displacement_8_y, 2)
            incremental_displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("INCREMENTAL_DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            incremental_displacement_8_y = incremental_displacements_8[1]
            self.assertAlmostEqual(0.025/2, incremental_displacement_8_y, 2)

            # integration point check in element 2, integration point 4 ( uniform stress and strain so an arbitrary choice )
            green_lagrange_strains_2_4    = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("GREEN_LAGRANGE_STRAIN_TENSOR", time, output_data[stage_nr], [2], [3])[0][0]
            green_lagrange_strains_2_4_yy = green_lagrange_strains_2_4[1]
            self.assertAlmostEqual(time*0.05, green_lagrange_strains_2_4_yy, 2)
            cauchy_stresses_2_4    = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data[stage_nr], [2], [3])[0][0]
            cauchy_stresses_2_4_yy = cauchy_stresses_2_4[1]
            self.assertAlmostEqual(E*time*0.05, cauchy_stresses_2_4_yy, 2)

        stage_nr = 1
        for time in [1.5, 2.0]:
            # displacement start at 0 at begin of second stage
            displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            displacement_8_y = displacements_8[1]
            self.assertAlmostEqual((time-1.0)*0.05/2, displacement_8_y, 2)
            # total displacement continuous over stages
            total_displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("TOTAL_DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            total_displacement_8_y = total_displacements_8[1]
            self.assertAlmostEqual(time*0.05/2, total_displacement_8_y, 2)
            # incremental displacement every step the same
            incremental_displacements_8  = test_helper.GiDOutputFileReader.nodal_values_at_time("INCREMENTAL_DISPLACEMENT", time, output_data[stage_nr], [8])[0]
            incremental_displacement_8_y = incremental_displacements_8[1]
            self.assertAlmostEqual(0.025/2, incremental_displacement_8_y, 2)
            # integration point check in element 2, integration point 4 ( uniform stress and strain so an arbitrary choice )
            # strains start at 0 at begin of second stage ( like the displacement )
            green_lagrange_strains_2_4 = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("GREEN_LAGRANGE_STRAIN_TENSOR", time, output_data[stage_nr], [2], [3])[0][0]
            green_lagrange_strains_2_4_yy = green_lagrange_strains_2_4[1]
            self.assertAlmostEqual((time-1.0)*0.05, green_lagrange_strains_2_4_yy, 2)
            # stresses are continuous over the stages
            cauchy_stresses_2_4    = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data[stage_nr], [2], [3])[0][0]
            cauchy_stresses_2_4_yy = cauchy_stresses_2_4[1]
            self.assertAlmostEqual(E*time*0.05, cauchy_stresses_2_4_yy, 2)

if __name__ == '__main__':
    KratosUnittest.main()
