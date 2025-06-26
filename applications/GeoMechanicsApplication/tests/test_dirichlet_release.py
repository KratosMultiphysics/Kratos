import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper


class KratosGeoMechanicsDirichletReleaseTests(KratosUnittest.TestCase):
    """
    This class contains a test for displacement Dirichlet boundary condition which is released in a second stage.
    """

    def test_dirichlet_release(self):
        """
        Single element test in 2 stages. The incremental elastic material should show
        linear increase in stress_yy in the first step. In the second stage the Dirichlet condition
        is removed, meaning the stress_yy should be 0.
        H = 1 [m]
        Young's modulus E = 10000 [N/m^2]
        stage 1) compress with 0.1 [m]
        stage 2) release the prescribed displacement, meaning the stress should be 0
        """
        test_name = "dirichlet_release"
        project_path = test_helper.get_file_path(test_name)
        n_stages = 2
        run_multiple_stages.run_stages(project_path, n_stages)

        # name of output file
        output_file_names = [
            os.path.join(project_path, f"dirichlet_release_stage{i + 1}.post.res")
            for i in range(n_stages)
        ]
        output_data = []

        reader = test_helper.GiDOutputFileReader()
        for output_file_name in output_file_names:
            output_data.append(reader.read_output_from(output_file_name))

        # Expected stress, resulting from the compression, Poisson effects (nu = 0.2) and
        # the strain in the x-direction (partially induced by the water pressure).
        expected_cauchy_stress_yy = -1016.67

        self.check_expected_outputs(
            expected_stage_displacement_and_strain=-0.1,
            expected_stress=expected_cauchy_stress_yy,
            output_data=output_data[0],
            time=1.0,
        )

        self.check_expected_outputs(
            expected_stage_displacement_and_strain=0.1072,
            expected_stress=100.0,
            output_data=output_data[1],
            time=2.0,
        )

    def check_expected_outputs(
        self,
        expected_stage_displacement_and_strain,
        expected_stress,
        output_data,
        time,
    ):
        # displacement of the top node
        displacement_top_node = test_helper.GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", time, output_data, [3]
        )[0]
        self.assertAlmostEqual(
            expected_stage_displacement_and_strain, displacement_top_node[1], 2
        )

        # integration point check in element 1, integration point 4 ( uniform stress and strain so an arbitrary choice )
        green_lagrange_strains_2_4 = (
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time(
                "GREEN_LAGRANGE_STRAIN_TENSOR", time, output_data, [1], [3]
            )[0][0]
        )
        green_lagrange_strains_2_4_yy = green_lagrange_strains_2_4[1]
        self.assertAlmostEqual(
            expected_stage_displacement_and_strain, green_lagrange_strains_2_4_yy, 2
        )
        cauchy_stresses_2_4 = (
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time(
                "CAUCHY_STRESS_TENSOR", time, output_data, [1], [3]
            )[0][0]
        )
        cauchy_stresses_2_4_yy = cauchy_stresses_2_4[1]
        self.assertAlmostEqual(expected_stress, cauchy_stresses_2_4_yy, 2)


if __name__ == "__main__":
    KratosUnittest.main()
