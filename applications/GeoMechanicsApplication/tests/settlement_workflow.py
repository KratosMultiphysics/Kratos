import os
import shutil
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
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
    def setUp(self):
        super().setUp()

        # The Kratos kernel uses a **static** list of applications that is implicitly shared by all kernel objects.
        # When a kernel object is destroyed, it will deregister all applications. Since the C++ route has its own
        # kernel object, the Python route may no longer work when the C++ route test is run first. To work around
        # this design flaw, explicitly import the GeoMechanicsApplication when a new test is about to be run:
        self.geo_app = GeoMechanicsApplication.KratosGeoMechanicsApplication()
        KratosMultiphysics._ImportApplication(self.geo_app, "KratosGeoMechanicsApplication")

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
    def setUp(self):
        super().setUp()

        self.settlement_api = GeoMechanicsApplication.CustomWorkflowFactory.CreateKratosGeoSettlement()


    def get_test_dir_name(self):
        return "cpp"


    def test_d_settlement_workflow(self):
        noop = lambda *args, **kwargs: None

        times_to_check = [1.0, 2.0, 3.0, 3.2]
        node_ids = [1, 102, 1085, 1442]
        reader = test_helper.GiDOutputFileReader()
        for i in range(self.number_of_stages):
            status = self.settlement_api.RunStage(self.test_path, self.project_parameters_filenames[i], noop, noop, noop, noop)
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


if __name__ == '__main__':
    KratosUnittest.main()
