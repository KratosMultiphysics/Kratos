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

        current_working_directory = os.getcwd()
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
        top_middle_node_id = 104
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_data, [top_middle_node_id]
        )[0][1]
        expected_settlement_after_one_hundred_days = -3.22

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            (
                expected_settlement_after_one_hundred_days
                - actual_settlement_after_one_hundred_days
            )
            / expected_settlement_after_one_hundred_days
            < 0.01
        )

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage5.post.res")
        )
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_data, [top_middle_node_id]
        )[0][1]

        expected_settlement_after_ten_thousand_days = -8.01

        # Assert the value to be within 1% of the analytical solution
        self.assertTrue(
            (
                expected_settlement_after_ten_thousand_days
                - actual_settlement_after_ten_thousand_days
            )
            / expected_settlement_after_one_hundred_days
            < 0.01
        )

        os.chdir(current_working_directory)


    def test_settlement_consolidation_coarse_mesh(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        pore pressure and stresses at specific time but several locations against expected results.
        The expected settlement values are based on regression tests.
        """
        test_name = "coarse_mesh"
        test_root = "dsettlement/consolidation_uniform_load"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        current_working_directory = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 3
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
        
        nodes = [3, 15, 16, 17, 18, 4]
        expected_pressure = [0.0, -18524.3, -24096.0, -33456.2, -42816.6, -50000.0]
        expected_total_stress = [-5332.28, -24858.1, -52279.1, -70030.9, -89823.8, -109627.0]
        
        for i in range(6):
            actual_pressure = reader.nodal_values_at_time("WATER_PRESSURE", 8726.4, output_data, [nodes[i]])[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)
        
            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]])[0][1]
            self.assertAlmostEqual(actual_total_stress, expected_total_stress[i], 4)

        os.chdir(current_working_directory)


    def test_settlement_consolidation_fine_mesh(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        pore pressure and stresses at specific time but several locations against expected results.
        The expected settlement values are based on regression tests.
        """
        test_name = "fine_mesh"
        test_root = "dsettlement/consolidation_uniform_load"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        current_working_directory = os.getcwd()
        os.chdir(project_path)

        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        n_stages = 3
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
        
        nodes = [87, 154, 266, 390, 516, 641]
        expected_pressure = [0.0, -16382.0, -24416.0, -33477.7, -42701.9, -50000.0]
        expected_total_stress = [-10324.4, -30053.2, -50021.0, -70014.2, -89998.3, -110028.0]
        
        for i in range(6):
            actual_pressure = reader.nodal_values_at_time("WATER_PRESSURE", 8726.4, output_data, [nodes[i]])[0]
            self.assertAlmostEqual(actual_pressure, expected_pressure[i], 4)
        
            actual_total_stress = reader.nodal_values_at_time(
                "TOTAL_STRESS_TENSOR", 8726.4, output_data, [nodes[i]])[0][1]
            self.assertAlmostEqual(actual_total_stress, expected_total_stress[i], 4)

        os.chdir(current_working_directory)

    def test_settlement_fully_saturated_column(self):
        """
        This test validates the settlement of a fully saturated column under uniform load.
        The test runs multiple stages of a settlement simulation and checks the
        settlement values at specific times against expected results.
        The expected settlement values are based on an analytical solution.
        """
        test_name = "fully_saturated_column_uniform_load"
        test_root = "dsettlement"
        project_path = test_helper.get_file_path(os.path.join(test_root, test_name))

        cwd = os.getcwd()
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
            os.path.join(project_path, "stage2.post.res")
        )
        actual_settlement_after_one_hundred_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 8640000, output_data, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_one_hundred_days, -1.71094, 4)

        output_data = reader.read_output_from(
            os.path.join(project_path, "stage5.post.res")
        )
        actual_settlement_after_ten_thousand_days = reader.nodal_values_at_time(
            "TOTAL_DISPLACEMENT", 864000000, output_data, [104]
        )[0][1]
        self.assertAlmostEqual(actual_settlement_after_ten_thousand_days, -8.63753, 4)

        os.chdir(cwd)

      
if __name__ == "__main__":
    KratosUnittest.main()
