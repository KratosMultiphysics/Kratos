import os
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper


class KratosGeoMechanicsSettlementWorkflow(KratosUnittest.TestCase):
    """
    This test class is used to check the settlement workflow test, same as test_settlement_workflow.cpp to
    make sure the python workflow yields the same results as the c++ workflow.
    """

    def test_d_settlement_workflow(self):
        test_name = 'test_settlement_workflow'
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_stages(file_path, 4)

        times_to_check = [1.0, 2.0, 3.0, 3.2]

        for i in range(4):
            result_file_name = os.path.join(file_path, f'test_model_stage{i+1}.post.res')
            expected_result_file_name = os.path.join(file_path, f'test_model_stage{i+1}.post.orig.res')

            reader = test_helper.GiDOutputFileReader()

            # These node ids are at the top corner, at the bottom of the excavation,
            # in the middle of the model and at the bottom.
            node_ids = [1, 102, 1085, 1442]

            actual_data = reader.read_output_from(result_file_name)
            actual_nodal_values = reader.nodal_values_at_time("DISPLACEMENT", times_to_check[i], actual_data, node_ids)

            expected_data = reader.read_output_from(expected_result_file_name)
            expected_nodal_values = reader.nodal_values_at_time(
                "DISPLACEMENT", times_to_check[i], expected_data, node_ids)

            self.assertEqual(len(actual_nodal_values), len(expected_nodal_values))
            for (actual_displacement, expected_displacement) in zip(actual_nodal_values, expected_nodal_values):
                self.assertVectorAlmostEqual(actual_displacement, expected_displacement, 3)

            if i > 2:
                self.check_stress_values(expected_result_file_name, i, node_ids, reader, result_file_name,
                                         times_to_check, "TOTAL_STRESS_TENSOR")
                self.check_stress_values(expected_result_file_name, i, node_ids, reader, result_file_name,
                                         times_to_check, "CAUCHY_STRESS_TENSOR")

    def check_stress_values(self, expected_result_file_name, i, node_ids, reader, result_file_name, times_to_check,
                            variable_name):
        actual_data = reader.read_output_from(result_file_name)
        actual_nodal_values = reader.nodal_values_at_time(variable_name, times_to_check[i], actual_data, node_ids)
        expected_data = reader.read_output_from(expected_result_file_name)
        expected_nodal_values = reader.nodal_values_at_time(
            variable_name, times_to_check[i], expected_data, node_ids)
        self.assertEqual(len(actual_nodal_values), len(expected_nodal_values))
        for (actual_total_stress, expected_total_stress) in zip(actual_nodal_values, expected_nodal_values):
            # Although the values are matrices, they are read as lists,
            # meaning we can use assertVectorAlmostEqual
            self.assertVectorAlmostEqual(actual_total_stress, expected_total_stress, 3)


if __name__ == '__main__':
    KratosUnittest.main()
