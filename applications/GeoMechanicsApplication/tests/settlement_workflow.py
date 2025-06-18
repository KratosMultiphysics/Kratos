import os
import shutil
import importlib

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the Kratos applications that we need
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

import test_helper


class KratosGeoMechanicsSettlementWorkflow(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("test_settlement_workflow")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

        shutil.rmtree(self.test_path, ignore_errors=True)

        os.makedirs(self.test_path)

        self.number_of_stages = 4
        self.project_parameters_filenames = [f"ProjectParameters_stage{i+1}.json" for i in range(self.number_of_stages)]
        input_filenames = self.project_parameters_filenames + ["MaterialParameters.json", "test_model.mdpa"]

        for filename in input_filenames:
            shutil.copy(os.path.join(self.test_root, filename), os.path.join(self.test_path, filename))

        # The expected values have been taken from a validated test run
        self.define_expected_displacements()
        self.define_expected_stresses()


    def get_test_dir_name(self):
        raise RuntimeError("This base class does not provide a generic test directory name")


    def define_expected_displacements(self):
        # The selected nodes are at the top corner, at the bottom of the excavation, in the middle of the model,
        # and at the bottom.
        self.expected_displacements = [{"output_filename": "test_model_stage1.post.res",
                                        "time": 1.0,
                                        "expected_values": [
                                            {"node": 1, "DISPLACEMENT": [0.0, -0.123882, 0.0]},
                                            {"node": 102, "DISPLACEMENT": [8.49409e-06, -0.119848, 0.0]},
                                            {"node": 1085, "DISPLACEMENT": [-7.25396e-09, -0.101875, 0.0]},
                                            {"node": 1442, "DISPLACEMENT": [0.0, 0.0, 0.0]}
                                        ]},
                                       {"output_filename": "test_model_stage2.post.res",
                                        "time": 2.0,
                                        "expected_values": [
                                            {"node": 1, "DISPLACEMENT": [0.0, -0.000162389, 0.0]},
                                            {"node": 102, "DISPLACEMENT": [8.50645e-05, -0.000467135, 0.0]},
                                            {"node": 1085, "DISPLACEMENT": [-0.00012122, -0.00186572, 0.0]},
                                            {"node": 1442, "DISPLACEMENT": [0.0, 0.0, 0.0]}
                                        ]},
                                       {"output_filename": "test_model_stage3.post.res",
                                        "time": 3.0,
                                        "expected_values": [
                                            {"node": 1, "DISPLACEMENT": [0.0, -0.000526393, 0.0]},
                                            {"node": 102, "DISPLACEMENT": [-4.47571e-05, -0.000954867, 0.0]},
                                            {"node": 1085, "DISPLACEMENT": [-0.000330983, -0.00372074, 0.0]},
                                            {"node": 1442, "DISPLACEMENT": [0.0, 0.0, 0.0]}
                                        ]},
                                       {"output_filename": "test_model_stage4.post.res",
                                        "time": 3.2,
                                        "expected_values": [
                                            {"node": 1, "DISPLACEMENT": [0.0, 0.0200705, 0.0]},
                                            {"node": 102, "DISPLACEMENT": [-0.000423281, 0.023303, 0.0]},
                                            {"node": 1085, "DISPLACEMENT": [0.000927345, 0.00102933, 0.0]},
                                            {"node": 1442, "DISPLACEMENT": [0.0, 0.0, 0.0]}
                                        ]}]


    def define_expected_stresses(self):
        self.expected_stresses = [{"output_filename": "test_model_stage4.post.res",
                                   "time": 3.2,
                                   "expected_values": [
                                       {"node": 1,
                                        "TOTAL_STRESS_TENSOR": [-7.69856, -0.127088, -1.56513, 0.127088, 0.0, 0.0],
                                        "CAUCHY_STRESS_TENSOR": [-7.69856, -0.127088, -1.56513, 0.127088, 0.0, 0.0]},
                                       {"node": 102,
                                        "TOTAL_STRESS_TENSOR": [-53.8354, -21.2794, -38.6748, 0.755572, 0.0, 0.0],
                                        "CAUCHY_STRESS_TENSOR": [-34.2154, -1.65944, -19.0548, 0.755572, 0.0, 0.0]},
                                       {"node": 1085,
                                        "TOTAL_STRESS_TENSOR": [-99.4165, -142.12, -102.597, -1.27531, 0.0, 0.0],
                                        "CAUCHY_STRESS_TENSOR": [-47.0336, -89.7376, -50.2142, -1.27531, 0.0, 0.0]},
                                       {"node": 1442,
                                        "TOTAL_STRESS_TENSOR": [-245.877, -319.876, -245.89, 1.67524, 0.0, 0.0],
                                        "CAUCHY_STRESS_TENSOR": [-108.537, -182.536, -108.55, 1.67524, 0.0, 0.0]}]
                                 }]


    def check_displacements(self):
        reader = test_helper.GiDOutputFileReader()

        for item in self.expected_displacements:
            time = item["time"]
            node_ids = [sub_item["node"] for sub_item in item["expected_values"]]
            expected_displacements = [sub_item["DISPLACEMENT"] for sub_item in item["expected_values"]]

            actual_data = reader.read_output_from(os.path.join(self.test_path, item["output_filename"]))
            actual_displacements = reader.nodal_values_at_time("DISPLACEMENT", time, actual_data, node_ids)

            self.assertEqual(len(actual_displacements), len(expected_displacements))
            for actual_displacement, expected_displacement in zip(actual_displacements, expected_displacements):
                self.assertVectorAlmostEqual(actual_displacement, expected_displacement, 3)


    def _check_nodal_stresses_with_name(self, output_data, stress_tensor_name, expected_stress_item):
        time = expected_stress_item["time"]
        node_ids = [sub_item["node"] for sub_item in expected_stress_item["expected_values"]]
        expected_nodal_stresses = [sub_item[stress_tensor_name] for sub_item in expected_stress_item["expected_values"]]

        actual_nodal_stresses = test_helper.GiDOutputFileReader.nodal_values_at_time(stress_tensor_name, time, output_data, node_ids)

        self.assertEqual(len(actual_nodal_stresses), len(expected_nodal_stresses))
        for actual_total_stress, expected_total_stress in zip(actual_nodal_stresses, expected_nodal_stresses):
            # Although the values are matrices, they are read as lists, meaning we can use assertVectorAlmostEqual
            self.assertVectorAlmostEqual(actual_total_stress, expected_total_stress, 3)


    def check_nodal_stresses(self):
        reader = test_helper.GiDOutputFileReader()

        for item in self.expected_stresses:
            actual_data = reader.read_output_from(os.path.join(self.test_path, item["output_filename"]))
            self._check_nodal_stresses_with_name(actual_data, "TOTAL_STRESS_TENSOR", item)
            self._check_nodal_stresses_with_name(actual_data, "CAUCHY_STRESS_TENSOR", item)


class KratosGeoMechanicsSettlementWorkflowPyRoute(KratosGeoMechanicsSettlementWorkflow):
    """
    This test class checks the settlement workflow using the Python route.
    """
    def get_test_dir_name(self):
        return "python"


    def test_d_settlement_workflow(self):
        test_helper.run_stages(self.test_path, self.number_of_stages)

        self.check_displacements()
        self.check_nodal_stresses()



class KratosGeoMechanicsSettlementWorkflowCppRoute(KratosGeoMechanicsSettlementWorkflow):
    """
    This test class checks the settlement workflow using the C++ route.
    """
    def get_test_dir_name(self):
        return "cpp"


    def test_d_settlement_workflow(self):
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement

        status = run_geo_settlement.run_stages(self.test_path, self.project_parameters_filenames)
        self.assertEqual(status, 0)

        self.check_displacements()
        self.check_nodal_stresses()


if __name__ == '__main__':
    KratosUnittest.main()
