import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper


class KratosGeoMechanicsDirichletUConstantTests(KratosUnittest.TestCase):
    """
    This class contains a test for displacement constant Dirichlet boundary condition and checks consistency of displacement flavours
    """

    def test_dirichlet_u_constant(self):
        """
        Single element compression test. The incremental elastic material should show
        constant total_strain_yy and stress_yy over the steps, the incremental_strain_yy should be 0 in step 2 of every stage
        H = 1 [m]
        Young's modulus E = 10E+03 [N/m^2]
        stage 0: compress -0.1 [m] in 1 steps and keep constant the next step
        stage 1: elongate  0.2 [m] in 1 steps and keep constant the next step
        stage 2: compress -0.1 [m] in 1 steps and keep constant the next step reset_displacement=true
        """
        test_name = "dirichlet_u_constant"
        project_path = test_helper.get_file_path(test_name)
        n_stages = 3
        run_multiple_stages.run_stages(project_path, n_stages)

        output_data = []
        reader = GiDOutputFileReader()
        for i in range(n_stages):
            output_data.append(
                reader.read_output_from(
                    os.path.join(
                        project_path, f"dirichlet_u_constant_stage{i+1}.post.res"
                    )
                )
            )

        E = 10e03

        # fmt: off
        expected_results = [(0, 0.5, -0.1, -0.1, -0.1, -0.1, -0.1 * E),
                            (0, 1.0, -0.1, -0.1,  0.0, -0.1, -0.1 * E),
                            (1, 1.5,  0.1,  0.2,  0.2,  0.2,  0.1 * E),
                            (1, 2.0,  0.1,  0.2,  0.0,  0.2,  0.1 * E),
                            (2, 2.5, -0.1, -0.1, -0.1, -0.1,  0.0 * E),
                            (2, 3.0, -0.1, -0.1,  0.0, -0.1,  0.0 * E),
                           ]
        # fmt: on
        for (
            stage_index,
            time,
            expected_total_uy,
            expected_stage_uy,
            expected_incremental_uy,
            expected_stage_eps_yy,
            expected_stage_sig_yy,
        ) in expected_results:
            # displacement of top node 3
            total_displacements_3 = (
                GiDOutputFileReader.nodal_values_at_time(
                    "TOTAL_DISPLACEMENT", time, output_data[stage_index], [3])[0]
            )
            self.assertAlmostEqual(expected_total_uy, total_displacements_3[1], 2)
            displacements_3 = GiDOutputFileReader.nodal_values_at_time(
                "DISPLACEMENT", time, output_data[stage_index], [3])[0]
            self.assertAlmostEqual(expected_stage_uy, displacements_3[1], 2)
            incremental_displacements_3 = (
                GiDOutputFileReader.nodal_values_at_time(
                    "INCREMENTAL_DISPLACEMENT", time, output_data[stage_index], [3])[0]
            )
            self.assertAlmostEqual(expected_incremental_uy, incremental_displacements_3[1], 2)

            # integration point check in element 1, integration point 4 ( uniform stress and strain so an arbitrary choice )
            green_lagrange_strains_1_4 = GiDOutputFileReader.element_integration_point_values_at_time(
                "GREEN_LAGRANGE_STRAIN_TENSOR", time, output_data[stage_index], [1], [3]
            )[0][0]
            self.assertAlmostEqual(expected_stage_eps_yy, green_lagrange_strains_1_4[1], 2)
            cauchy_stresses_1_4 = GiDOutputFileReader.element_integration_point_values_at_time(
                "CAUCHY_STRESS_TENSOR", time, output_data[stage_index], [1], [3]
            )[0][0]
            self.assertAlmostEqual(expected_stage_sig_yy, cauchy_stresses_1_4[1], 2)


if __name__ == "__main__":
    KratosUnittest.main()
