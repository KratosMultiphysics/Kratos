import sys
import os
import json
import math
import csv

sys.path.append(os.path.join('..', '..', '..', 'bin', 'Release'))

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class TestConsecutivePipeLines(KratosUnittest.TestCase):
    """
    Checking to see if consecutive pipelines with different configurations, work.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.is_running_under_teamcity = test_helper.is_running_under_teamcity()

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def change_head_level_polder_side(self, file_path, head_level):
        parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            for process in parameters['processes']['constraints_process_list']:
                if "Left_head" in "Left_head" in process["Parameters"]["model_part_name"]:
                    process['Parameters']['reference_coordinate'] = head_level
            parameters['processes']['constraints_process_list'][0]['Parameters']['reference_coordinate'] = head_level
        with open(parameter_file_name, 'w') as parameter_file:
            json.dump(parameters, parameter_file, indent=4)

    def model_kratos_run(self, file_path, head):
        self.change_head_level_polder_side(file_path, head)
        simulation = test_helper.run_kratos(file_path)
        pipe_active = test_helper.get_pipe_active_in_elements(simulation)
        return all(pipe_active)

    def linear_search(self, file_path, search_array):
        counter_head = 0
        while counter_head < len(search_array):
            # check if pipe elements become active
            if self.model_kratos_run(file_path, search_array[counter_head]):
                return search_array[counter_head - 1]
            counter_head = counter_head + 1
        return math.nan

    @KratosUnittest.skip("This test fails and should be investigated and solved under KRATOSGEO-36")
    def test_consecutive_pipe_lines(self):
        test_files = ["reference_geometry", "split_geometry", "split_geometry_same_pipe", "split_geometry_all_same"]
        result_dict = {}

        for test_file in test_files:
            test_name = os.path.join('test_consecutive_pipe_lines', test_file)
            file_path = test_helper.get_file_path(os.path.join('.', test_name))
            heads = [x * 1 for x in range(int(0),
                                          int(10), 1)]
            critical_head_found = self.linear_search(file_path, heads)
            result_dict[test_file] = critical_head_found

        assert math.isclose(result_dict["reference_geometry"], result_dict["split_geometry"], result_dict["split_geometry_same_pipe"], result_dict["split_geometry_all_same"])
