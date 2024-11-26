import os
import shutil
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the Kratos applications that we need. We need to reload these later.
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

import test_helper


class KratosGeoMechanicsSettlementWorkflow(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("test_settlement_workflow")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

        try:
            shutil.rmtree(self.test_path)
        except FileNotFoundError:
            pass

        os.makedirs(self.test_path)

        self.number_of_stages = 4
        self.project_parameters_filenames = [f"ProjectParameters_stage{i+1}.json" for i in range(self.number_of_stages)]
        input_filenames = self.project_parameters_filenames[:] + ["MaterialParameters.json", "test_model.mdpa"]

        for filename in input_filenames:
            shutil.copy(os.path.join(self.test_root, filename), os.path.join(self.test_path, filename))

    def get_test_dir_name(self):
        raise RuntimeError("This base class does not provide a generic test directory name")


class KratosGeoMechanicsSettlementWorkflowPyRoute(KratosGeoMechanicsSettlementWorkflow):
    """
    This test class is used to check the settlement workflow test, same as test_settlement_workflow.cpp to
    make sure the python workflow yields the same results as the c++ workflow.
    """
    def get_test_dir_name(self):
        return "python"


    def test_d_settlement_workflow(self):
        test_helper.run_stages(self.test_path, self.number_of_stages)

        times_to_check = [1.0, 2.0, 3.0, 3.2]

        for i in range(self.number_of_stages):
            result_file_name = os.path.join(self.test_path, f'test_model_stage{i+1}.post.res')
            expected_result_file_name = os.path.join(self.test_root, f'test_model_stage{i+1}.post.orig.res')

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
            for actual_displacement, expected_displacement in zip(actual_nodal_values, expected_nodal_values):
                self.assertVectorAlmostEqual(actual_displacement, expected_displacement, 3)

            if i > 2:
                self.check_stress_values(expected_result_file_name, times_to_check[i], node_ids, reader,
                                         result_file_name, "TOTAL_STRESS_TENSOR")
                self.check_stress_values(expected_result_file_name, times_to_check[i], node_ids, reader,
                                         result_file_name, "CAUCHY_STRESS_TENSOR")

    def check_stress_values(self, expected_result_file_name, time_to_check, node_ids, reader, result_file_name,
                            variable_name):
        actual_data = reader.read_output_from(result_file_name)
        actual_nodal_stress_values = reader.nodal_values_at_time(variable_name, time_to_check, actual_data, node_ids)
        expected_data = reader.read_output_from(expected_result_file_name)
        expected_nodal_stress_values = reader.nodal_values_at_time(
            variable_name, time_to_check, expected_data, node_ids)
        self.assertEqual(len(actual_nodal_stress_values), len(expected_nodal_stress_values))
        for actual_total_stress, expected_total_stress in zip(actual_nodal_stress_values, expected_nodal_stress_values):
            # Although the values are matrices, they are read as lists,
            # meaning we can use assertVectorAlmostEqual
            self.assertVectorAlmostEqual(actual_total_stress, expected_total_stress, 3)


class KratosGeoMechanicsSettlementWorkflowCppRoute(KratosGeoMechanicsSettlementWorkflow):
    """
    This test class is used to check the settlement workflow test, same as test_settlement_workflow.cpp to
    make sure the python workflow yields the same results as the c++ workflow.
    """
    def tearDown(self):
        # The `KratosGeoSettlement` instance used by this test removes all registered GeoMechanicsApplication
        # components when it's destroyed. If no action is taken, any following tests will start to fail. It seems
        # that reloading the relevant Kratos applications overcomes this problem.
        importlib.reload(KratosMultiphysics.LinearSolversApplication)
        importlib.reload(KratosMultiphysics.StructuralMechanicsApplication)
        importlib.reload(KratosMultiphysics.GeoMechanicsApplication)


    def get_test_dir_name(self):
        return "cpp"


    def test_d_settlement_workflow(self):
        settlement_api = GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoSettlement()

        no_logging = lambda msg: None
        no_progress_reporting = lambda fraction_done: None
        no_progress_message = lambda msg: None
        dont_cancel = lambda: False

        times_to_check = [1.0, 2.0, 3.0, 3.2]
        node_ids = [1, 102, 1085, 1442]
        reader = test_helper.GiDOutputFileReader()
        for i in range(self.number_of_stages):
            status = settlement_api.RunStage(self.test_path, self.project_parameters_filenames[i], no_logging, no_progress_reporting, no_progress_message, dont_cancel)
            self.assertEqual(status, 0)

            result_file_name = os.path.join(self.test_path, f'test_model_stage{i+1}.post.res')
            expected_result_file_name = os.path.join(self.test_root, f'test_model_stage{i+1}.post.orig.res')

            actual_data = reader.read_output_from(result_file_name)
            actual_nodal_values = reader.nodal_values_at_time("DISPLACEMENT", times_to_check[i], actual_data, node_ids)

            expected_data = reader.read_output_from(expected_result_file_name)
            expected_nodal_values = reader.nodal_values_at_time(
                "DISPLACEMENT", times_to_check[i], expected_data, node_ids)

            self.assertEqual(len(actual_nodal_values), len(expected_nodal_values))
            for actual_displacement, expected_displacement in zip(actual_nodal_values, expected_nodal_values):
                self.assertVectorAlmostEqual(actual_displacement, expected_displacement, 3)

        # Don't rely on the garbage collector to clean up the API object. Make sure it's destructor has run before
        # executing the test case's `tearDown` method (which will reload the relevant Kratos applications)
        del settlement_api


if __name__ == '__main__':
    KratosUnittest.main()
