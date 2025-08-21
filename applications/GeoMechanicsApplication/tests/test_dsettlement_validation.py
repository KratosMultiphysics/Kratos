from KratosMultiphysics.GeoMechanicsApplication import run_multiple_stages
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

import test_helper


class KratosGeoMechanicsDSettlementValidationTests(KratosUnittest.TestCase):
    def test_settlement_dry_column(self):
        """
        This test validates the settlement of a dry column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "dry_column_uniform_load"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        original_working_dir = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 5
        project_parameters_filenames = [
            f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)
        ]
        status = run_geo_settlement.run_stages(
            project_path, project_parameters_filenames
        )
        self.assertEqual(status, 0)

        reader = test_helper.GiDOutputFileReader()

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage3.post.res")
        )
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_data, [104]
        )[0][1]
        expected_settlement_after_one_hundred_days = -3.22
        # Within 1% of the expected value
        self.assertTrue((expected_settlement_after_one_hundred_days - actual_settlement_after_one_hundred_days)/expected_settlement_after_one_hundred_days < 0.01)

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage5.post.res")
        )
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_data, [104]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01
        self.assertTrue((expected_settlement_after_ten_thousand_days - actual_settlement_after_ten_thousand_days)/expected_settlement_after_one_hundred_days < 0.01)


        os.chdir(original_working_dir)


if __name__ == "__main__":
    KratosUnittest.main()
