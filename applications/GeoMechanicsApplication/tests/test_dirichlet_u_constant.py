import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsDirichletUConstantTests(KratosUnittest.TestCase):
    """
    This class contains a test for displacement constant Dirichlet boundary condition and checks consistency of displacement flavours
    """
    def test_dirichlet_u_constant(self):
        """
        4 element elongation test in 2 stages. The incremental elastic material should show 
        constant total_strain_yy and stress_yy over the steps, the incremental_strain_yy should be 0 from step 2
        H = 1 [m]
        Young's modulus E = 30E+06 [N/m^2]
        compress -0.1 [m] in 1 steps and keep constant the next step
        """
        test_name    = 'dirichlet_u_constant'
        project_path = test_helper.get_file_path(test_name)
        n_stages     = 1
        run_multiple_stages.run_stages(project_path, n_stages)

        output_data = []
        reader = test_helper.GiDOutputFileReader()
        for i in range(n_stages):
            output_data.append(reader.read_output_from(os.path.join(project_path, f'dirichlet_u_constant_stage{i+1}.post.res')))

        E = 10E+03
        stage_nr = 0
        for time in [0.5, 1.0]:
            # displacement of a top node
            displacements_3  = test_helper.GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", time, output_data[stage_nr], [3])[0]
            displacement_3_y = displacements_3[1]
            self.assertAlmostEqual(-0.1, displacement_3_y, 2)
            total_displacements_3  = test_helper.GiDOutputFileReader.nodal_values_at_time("TOTAL_DISPLACEMENT", time, output_data[stage_nr], [3])[0]
            total_displacement_3_y = total_displacements_3[1]
            self.assertAlmostEqual(-0.1, total_displacement_3_y, 2)
            incremental_displacements_3  = test_helper.GiDOutputFileReader.nodal_values_at_time("INCREMENTAL_DISPLACEMENT", time, output_data[stage_nr], [3])[0]
            incremental_displacement_3_y = incremental_displacements_3[1]
            if time == 0.5:
                self.assertAlmostEqual(-0.1, incremental_displacement_3_y, 2)
            else:
                self.assertAlmostEqual(0.0, incremental_displacement_3_y, 2)

            # integration point check in element 1, integration point 4 ( uniform stress and strain so an arbitrary choice )
            green_lagrange_strains_1_4    = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("GREEN_LAGRANGE_STRAIN_TENSOR", time, output_data[stage_nr], [1], [3])[0][0]
            green_lagrange_strains_1_4_yy = green_lagrange_strains_1_4[1]
            self.assertAlmostEqual(-0.1, green_lagrange_strains_1_4_yy, 2)
            cauchy_stresses_1_4    = test_helper.GiDOutputFileReader.element_integration_point_values_at_time("CAUCHY_STRESS_TENSOR", time, output_data[stage_nr], [1], [3])[0][0]
            cauchy_stresses_1_4_yy = cauchy_stresses_1_4[1]
            self.assertAlmostEqual(E*-0.1, cauchy_stresses_1_4_yy, 2)

if __name__ == '__main__':
    KratosUnittest.main()
