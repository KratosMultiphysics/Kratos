import sys
import os
import json
import math
from parameterized import parameterized

import KratosMultiphysics.KratosUnittest as KratosUnittest

sys.path.append(os.path.join('..', 'python_scripts'))
import test_helper

class log_writer:
    def __init__(self, filename):
        self.filename = filename
        self.fo = open(filename, "w")

    def write_log(self, message):
        self.fo.write(message + "\n")
        self.fo.flush()

    def close(self):
        self.fo.close()

class LatexWriterFile:
    """ Class that writes results for latex documentation """

    def __init__(self, filename=""):
        self.filename = filename
        self.value_dict_default = {"value_name": "", "test_result": 0, "kratos_results": 0}

    def write_latex_file(self, result_pair):
        with open(self.filename, "a") as output_latex_file:
            error_height = abs(result_pair['kratos_results_h'] - result_pair['test_result_h']) / \
                           (abs(result_pair['test_result_h']) + 1e-60)
            error_equivalent_software_height = abs(
                result_pair['kratos_results_h'] - result_pair['equivalent_software_h']) / \
                                               (abs(result_pair['equivalent_software_h']) + 1e-60)
            error_equivalent_software_length = abs(
                result_pair['kratos_results_l'] - result_pair['equivalent_software_l']) / \
                                               (abs(result_pair['equivalent_software_l']) + 1e-60)
            output_latex_file.write(
                f"{result_pair['value_name']} & {result_pair['test_result_h']} & "
                f"{result_pair['equivalent_software_h']} &  {round(result_pair['kratos_results_h'], 2)} & "
                f" {round(error_height * 100, 2)} &  {round(error_equivalent_software_height * 100, 2)} & "
                f" {round(result_pair['equivalent_software_l'], 2)} &  {round(result_pair['kratos_results_l'], 2)} & "
                f" {round(error_equivalent_software_length * 100, 2)} \\\\ \\hline \n")


class TestSellmeijersRuleValidation(KratosUnittest.TestCase):

    # add code that is executed before all tests once
    @classmethod
    def setUpClass(self):
        print("Setting Up Test Suite")
        self.logger = log_writer("test_log.txt")

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.latex_writer = LatexWriterFile()
        self.results = {}

    def tearDown(self):
        self.latex_writer.filename = test_helper.get_file_path('test_compare_sellmeijer_B/test_compare_sellmeijer.tex')
        self.latex_writer.write_latex_file(self.results)

    @classmethod
    def tearDownClass(self):
        self.logger.close()

    def change_material_parameters(self, file_path, kappa, d70):
        # change the values of the pipe elements
        material_file_name = os.path.join(file_path, 'MaterialParameters.json')
        with open(material_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            parameters['properties'][0]['Material']['Variables']['PERMEABILITY_XX'] = kappa
            parameters['properties'][0]['Material']['Variables']['PERMEABILITY_YY'] = kappa
            parameters['properties'][2]['Material']['Variables']['PIPE_D_70'] = d70
        with open(material_file_name, 'w') as parameter_file:
            json.dump(parameters, parameter_file, indent=4)

    def change_head_level_polder_side(self, file_path, head_level):
        parameter_file_name = os.path.join(file_path, 'ProjectParameters.json')
        with open(parameter_file_name, 'r') as parameter_file:
            parameters = json.load(parameter_file)
            if "Left" in parameters['processes']['constraints_process_list'][0]['Parameters']['model_part_name']:
                parameters['processes']['constraints_process_list'][0]['Parameters'][
                    'reference_coordinate'] = head_level
            else:
                parameters['processes']['constraints_process_list'][1]['Parameters'][
                    'reference_coordinate'] = head_level
        with open(parameter_file_name, 'w') as parameter_file:
            json.dump(parameters, parameter_file, indent=4)

    def model_kratos_run(self, file_path, head, model=None):
        # if model is not None:
        #     model.Reset()
        # open file to append to
        self.logger.write_log("HEAD :" + str(head))
        self.change_head_level_polder_side(file_path, head)
        self.logger.write_log("Head adjusted")
        simulation = test_helper.run_kratos(file_path, model)
        self.logger.write_log("Run ended")
        pipe_active = test_helper.get_pipe_active_in_elements(simulation)
        length = test_helper.get_pipe_length(simulation)
        self.logger.write_log("PIPE ACTIVE :" + str(pipe_active))
        self.logger.write_log("PIPE LENGTH :" + str(length))

        return all(pipe_active), length, None

    def find_critical_head(self, file_path, head):
        model = None
        last_head = 0.0
        last_length = 0.0

        while True:
            # check if pipe elements become active
            pipe_active, pipe_length, model = self.model_kratos_run(file_path, head, model)

            if pipe_active:
                return last_head, last_length
            else:
                last_head = head
                last_length = pipe_length
                head += 0.1


    def drange(self, start, stop, step):
        r = start
        while r < stop:
            yield r
            r += step

    def critical_head_loop(self, file_path, kappa, d70, Hc):
        self.change_material_parameters(file_path, kappa, d70)
        critical_head_found = math.nan
        length = math.nan
        critical_head_found, length = self.find_critical_head(file_path, Hc)
        return critical_head_found, length

    #@parameterized.expand(
    #    [('7.1', 1.00E-04, 1.157E-12, 3.43, 3.7, 6, 'test_compare_sellmeijer_B/HeightAquiferD10L30.gid')
    #     ])
    def test_sellmeijers_rule_height(self, name, d70, kappa, Hc, Hn, length_n, test_name_gid, isTest=False):
        self.logger.write_log("Test started: " + name + "\n")
        file_path = test_helper.get_file_path(os.path.join('./', test_name_gid))
        os.chdir(file_path)
        critical_head_found, length = self.critical_head_loop(file_path, kappa, d70, max(0.0, Hn - 1))
        self.results = {"value_name": name, "test_result_h": Hc, "kratos_results_h": critical_head_found,
                        "equivalent_software_h": Hn, "kratos_results_l": length, "equivalent_software_l": length_n}

        if critical_head_found is None and isTest:
            self.assertTrue(False, f"Critical head not found for {name}")
        else:
            print(f"Critical head found for {name} is {critical_head_found}", flush = True)
            print(f"Length found for {name} is {length}", flush = True)
            self.logger.write_log(f"Critical head found for {name} is {critical_head_found}")
            self.logger.write_log(f"Length found for {name} is {length}")

        if isTest:
            self.assertAlmostEqual(Hn, critical_head_found, 1, f"Critical head kratos: {critical_head_found}, old geo flow {Hn}")
            self.assertAlmostEqual(length_n, length, 1, f"Length kratos: {length}, old geo flow {length_n}")

if __name__ == '__main__':
    #KratosUnittest.main()
    test = TestSellmeijersRuleValidation()
    test.logger = log_writer("test_log.txt")
    test.test_sellmeijers_rule_height("7.1", 1.00E-04, 1.157E-12, 3.43, 3.7, 6, 'test_compare_sellmeijer/HeightAquiferD10L30.gid', True)
    #test.test_sellmeijers_rule_height('7.2', 1.00E-04, 1.157E-12, 6.37, 7.4, 12, 'test_compare_sellmeijer_B/HeightAquiferD10L60.gid', True)
    #test.test_sellmeijers_rule_height('7.3', 1.00E-04, 1.157E-12, 9.18, 11.2, 13.5, 'test_compare_sellmeijer_B/HeightAquiferD10L90.gid')
    #test.test_sellmeijers_rule_height('7.4', 1.00E-04, 1.157E-12, 3, 3.3, 9, 'test_compare_sellmeijer_B/HeightAquiferD20L30.gid')
    #test.test_sellmeijers_rule_height('7.5', 1.00E-04, 1.157E-12, 5.44, 6.1, 15, 'test_compare_sellmeijer_B/HeightAquiferD20L60.gid')