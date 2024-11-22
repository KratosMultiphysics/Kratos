import os
import shutil

import test_helper
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication


class KratosGeoMechanicsSettlementWorkflowCppRoute(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("test_settlement_workflow")
        self.test_path = os.path.join(self.test_root, "cpp")

        try:
            shutil.rmtree(self.test_path)
        except FileNotFoundError:
            pass

        os.makedirs(self.test_path)

        self.project_parameters_filenames = ["ProjectParameters_stage1.json", "ProjectParameters_stage2.json", "ProjectParameters_stage3.json", "ProjectParameters_stage4.json"]
        input_filenames = self.project_parameters_filenames[:] + ["MaterialParameters.json", "test_model.mdpa"]
        for filename in input_filenames:
            shutil.copy(os.path.join(self.test_root, filename), os.path.join(self.test_path, filename))


    def test_d_settlement_workflow(self):
        noop = lambda *args, **kwargs: None
        settlement_api = GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoSettlement()

        times_to_check = [1.0, 2.0, 3.0, 3.2]
        node_ids = [1, 102, 1085, 1442]

        reader = test_helper.GiDOutputFileReader()

        for index, filename in enumerate(self.project_parameters_filenames):
            status = settlement_api.RunStage(self.test_path, filename, noop, noop, noop, noop)
            self.assertEqual(status, 0)

            result_file_name = os.path.join(self.test_path, f'test_model_stage{index+1}.post.res')
            expected_result_file_name = os.path.join(self.test_root, f'test_model_stage{index+1}.post.orig.res')

            actual_data = reader.read_output_from(result_file_name)
            actual_nodal_values = reader.nodal_values_at_time("DISPLACEMENT", times_to_check[index], actual_data, node_ids)

            expected_data = reader.read_output_from(expected_result_file_name)
            expected_nodal_values = reader.nodal_values_at_time(
                "DISPLACEMENT", times_to_check[index], expected_data, node_ids)

            self.assertEqual(len(actual_nodal_values), len(expected_nodal_values))
            for actual_displacement, expected_displacement in zip(actual_nodal_values, expected_nodal_values):
                self.assertVectorAlmostEqual(actual_displacement, expected_displacement, 3)


if __name__ == '__main__':
    KratosUnittest.main()
