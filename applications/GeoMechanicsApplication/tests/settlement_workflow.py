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

        for i in range(4):
            result_file_name = os.path.join(file_path, 'test_model_stage' + str(i + 1) + '.post.res')
            expected_result_file_name = test_helper.get_file_path(os.path.join('.', test_name, 'common', 'test_model_stage' + str(i + 1) + '.post.orig.res'))
            self.assertTrue(filecmp.cmp(result_file_name, expected_result_file_name, shallow=False))

if __name__ == '__main__':
    KratosUnittest.main()
