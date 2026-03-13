import os
import json
import math

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class TestConsecutivePipeLines(KratosUnittest.TestCase):
    """
     Class to Check consecutive pipelines with different configurations
     Note that  we are testing for convergence and FE operation/integration, No theory is available to describe 2
     materials, therefore no comparison of results.
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def change_head_level_polder_side(self, file_path, head_level):
        parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            for process in parameters['processes']['constraints_process_list']:
                if "Left_head" in process["Parameters"]["model_part_name"]:
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


    def test_reference_geometry(self):
        """
         Testing a simple Sellmeijer rule (Soil_permeability= 1.157e-12 , PIPE_D_70: 0.0001)
        """
        test_file = "reference_geometry"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(370),
                                      int(380), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 3.78, 2)


    def test_split_geometry_permeability_soil1_e10(self):
        """
         Testing a split geometry with two different permeability (Soil_1_permeability= 1.157e-10 ,
         Soil_2_permeability= 1.157e-12)
        """
        test_file = "split_geometry_permeability_soil1_e10"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(165),
                                      int(180), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 1.67, 2)


    def test_split_geometry_permeability_soil2_e10(self):
        """
         Testing a split geometry with two different permeability (Soil_1_permeability= 1.157e-12 ,
         Soil_2_permeability= 1.157e-10)
        """
        test_file = "split_geometry_permeability_soil2_e10"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(180),
                                      int(190), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 1.82, 2)


    def test_split_geometry_permeability_soil1_soil2_e10(self):
        """
         Testing a split geometry with same permeability (Soil_permeability= 1.157e-10)
        """
        test_file = "split_geometry_permeability_soil1_soil2_e10"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(90),
                                         int(100), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 0.93, 2)


    def test_split_geometry_double_lines_pipe2_added(self):
        """
         Testing a split geometry with same permeability (Soil_permeability= 1.157e-12). Pipe is separated to two parts.
         Pipe_2 with the same material is added to the geometry along the other one (PIPE_D_70: 0.0001)
        """
        test_file = "split_geometry_double_lines_pipe2_added"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(370),
                                      int(380), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 3.77, 2)


    def test_split_geometry_pipe2_D70_3e4(self):
        """
         Testing a split geometry with same permeability (Soil_permeability= 1.157e-10). Pipe is separated to two parts
         with different D70 parameter (PIPE_D_70: 0.0003)
        """
        test_file = "split_geometry_pipe2_D70_3e4"

        test_name = os.path.join('test_consecutive_pipe_lines', test_file)
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        heads = [x * 0.01 for x in range(int(240),
                                      int(250), 1)]
        critical_head_found = self.linear_search(file_path, heads)
        self.assertAlmostEqual(critical_head_found, 2.43, 2)

if __name__ == '__main__':
    KratosUnittest.main()