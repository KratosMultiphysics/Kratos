from KratosMultiphysics.GeoMechanicsApplication import run_multiple_stages
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

import test_helper

class KratosGeoMechanicsDSettlementTests(KratosUnittest.TestCase):
    def test_settlement_dry_column(self):
        test_name="dry_column_uniform_load"
        test_relative_path=os.path.join("dsettlement",test_name)
        project_path = test_helper.get_file_path(test_relative_path)
        n_stages     = 5
        cwd=os.getcwd()
        os.chdir(project_path)
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement
        self.project_parameters_filenames = [f"ProjectParameters_stage{i+1}.json" for i in range(n_stages)]
        status = run_geo_settlement.run_stages(project_path, self.project_parameters_filenames)
        self.assertEqual(status, 0)
        reader=test_helper.GiDOutputFileReader()
        output_data=reader.read_output_from(os.path.join(project_path, "stage3.post.res"))
        actual_settlement_after_one_hundred_days=reader.nodal_values_at_time("TOTAL_DISPLACEMENT", 8640000, output_data, [12])[0][1]
        expected_settlement_after_one_hundred_days=-1.75393
        output_data=reader.read_output_from(os.path.join(project_path, "stage5.post.res"))
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, expected_settlement_after_one_hundred_days, 4)
        actual_settlement_after_ten_thousand_days=reader.nodal_values_at_time("TOTAL_DISPLACEMENT", 864000000, output_data, [12])[0][1]

        expected_settlement_after_ten_thousand_days = -6.4576
        self.assertAlmostEqual(actual_settlement_after_ten_thousand_days, expected_settlement_after_ten_thousand_days, 4)

        
        os.chdir(cwd)


if __name__ == '__main__':
    KratosUnittest.main()
            