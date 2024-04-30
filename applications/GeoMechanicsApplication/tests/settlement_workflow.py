# from KratosMultiphysics import * as Kratos

import os
import filecmp
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsSettlementWorkflow(KratosUnittest.TestCase):
    """
    This test class is used to check the settlement workflow test, same as test_settlement_workflow.cpp to
    make sure the python workflow yields the same results as the c++ workflow.
    """

    def test_DSettlement_workflow(self):
        test_name = 'test_settlement_workflow'
        file_path = test_helper.get_file_path(os.path.join('.', test_name, 'python_workflow'))
        test_helper.run_stages(file_path, 4)

        times_to_check = [1.0, 2.0, 3.0, 3.2]

        for i in range(4):
            result_file_name = os.path.join(file_path, 'test_model_stage' + str(i + 1) + '.post.res')
            expected_result_file_name = (test_helper.get_file_path(
                os.path.join('.', test_name, 'common', 'test_model_stage' + str(i + 1) + '.post.orig.res')))

            reader = test_helper.GiDOutputFileReader()

            # These node ids are at the top corner, at the bottom of the excavation,
            # in the middle of the model and at the bottom.
            node_ids = [1, 102, 1085, 1442]

            actual_data = reader.read_output_from(result_file_name)
            actual_nodal_values = reader.nodal_values_at_time("DISPLACEMENT", times_to_check[i], actual_data, [1])

            expected_data = reader.read_output_from(expected_result_file_name)
            expected_nodal_values = reader.nodal_values_at_time("DISPLACEMENT", times_to_check[i], expected_data, [1])

            self.assertEqual(len(actual_nodal_values), len(expected_nodal_values))
            for j in range(len(actual_nodal_values)):
                self.assertAlmostEqual(actual_nodal_values[j][0], expected_nodal_values[j][0], 3)
                self.assertAlmostEqual(actual_nodal_values[j][1], expected_nodal_values[j][1], 3)
                self.assertAlmostEqual(actual_nodal_values[j][2], expected_nodal_values[j][2], 3)




if __name__ == '__main__':
    KratosUnittest.main()
