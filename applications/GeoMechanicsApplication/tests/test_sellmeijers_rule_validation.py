import sys
import os
import json
import math
from parameterized import parameterized

import KratosMultiphysics.KratosUnittest as KratosUnittest

sys.path.append(os.path.join('..', 'python_scripts'))
import test_helper


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

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.latex_writer = LatexWriterFile()
        self.results = {}

    def tearDown(self):
        self.latex_writer.filename = test_helper.get_file_path('test_compare_sellmeijer/test_compare_sellmeijer.tex')
        self.latex_writer.write_latex_file(self.results)

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
        self.change_head_level_polder_side(file_path, head)
        simulation = test_helper.run_kratos(file_path, model)
        pipe_active = test_helper.get_pipe_active_in_elements(simulation)
        length = test_helper.get_pipe_length(simulation)

        model = simulation.model
        return all(pipe_active), length, model

    def linear_search(self, file_path, search_array):
        counter_head = 0
        pipe_length_total = []
        model = None
        while counter_head < len(search_array):
            # check if pipe elements become active
            pipe_active, pipe_length, model = self.model_kratos_run(file_path, search_array[counter_head], model)
            pipe_length_total.append(pipe_length)
            if pipe_active:
                return search_array[counter_head - 1], pipe_length_total[counter_head - 1]
            counter_head = counter_head + 1
        return None, None

    def drange(self, start, stop, step):
        r = start
        while r < stop:
            yield r
            r += step

    def critical_head_loop(self, file_path, kappa, d70, Hc, search_type='linear'):
        self.change_material_parameters(file_path, kappa, d70)
        heads = list(self.drange(Hc - 1, Hc + 2, 0.1))
        critical_head_found = math.nan
        length = math.nan
        if search_type == 'linear':
            critical_head_found, length = self.linear_search(file_path, heads)
        return critical_head_found, length

    @parameterized.expand(
        [('7.1', 1.00E-04, 1.157E-12, 3.43, 3.7, 6, 'test_compare_sellmeijer/HeightAquiferD10L30.gid'),
         ('7.2', 1.00E-04, 1.157E-12, 6.37, 7.4, 12, 'test_compare_sellmeijer/HeightAquiferD10L60.gid'),
         ('7.3', 1.00E-04, 1.157E-12, 9.18, 11.2, 13.5, 'test_compare_sellmeijer/HeightAquiferD10L90.gid'),
         ('7.4', 1.00E-04, 1.157E-12, 3, 3.3, 9, 'test_compare_sellmeijer/HeightAquiferD20L30.gid'),
         ('7.5', 1.00E-04, 1.157E-12, 5.44, 6.1, 15, 'test_compare_sellmeijer/HeightAquiferD20L60.gid'),
         ('7.6', 1.00E-04, 1.157E-12, 7.81, 9.1, 21, 'test_compare_sellmeijer/HeightAquiferD20L90.gid'),
         ('7.7', 1.00E-04, 1.157E-12, 2.6, 3.1, 10.5, 'test_compare_sellmeijer/HeightAquiferD30L30.gid'),
         ('7.8', 1.00E-04, 1.157E-12, 5.02, 5.7, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60.gid'),
         ('7.9', 1.00E-04, 1.157E-12, 7.13, 8.1, 24, 'test_compare_sellmeijer/HeightAquiferD30L90.gid'),
         ('7.10', 3.00E-04, 1.157E-12, 10.29, 11.3, 7.5, 'test_compare_sellmeijer/HeightAquiferD10L30.gid'),
         ('7.11', 3.00E-04, 1.157E-12, 19.1, 21.9, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60.gid'),
         ('7.12', 3.00E-04, 1.157E-12, 27.54, 32.8, 15, 'test_compare_sellmeijer/HeightAquiferD10L90.gid'),
         ('7.13', 3.00E-04, 1.157E-12, 9.01, 9.9, 10.5, 'test_compare_sellmeijer/HeightAquiferD20L30.gid'),
         ('7.14', 3.00E-04, 1.157E-12, 16.33, 18.2, 16.5, 'test_compare_sellmeijer/HeightAquiferD20L60.gid'),
         ('7.15', 3.00E-04, 1.157E-12, 23.42, 27, 21, 'test_compare_sellmeijer/HeightAquiferD20L90.gid'),
         ('7.16', 3.00E-04, 1.157E-12, 7.8, 9.4, 12, 'test_compare_sellmeijer/HeightAquiferD30L30.gid'),
         ('7.17', 3.00E-04, 1.157E-12, 15.05, 16.8, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60.gid'),
         ('7.18', 3.00E-04, 1.157E-12, 21.4, 24.1, 28.5, 'test_compare_sellmeijer/HeightAquiferD30L90.gid'),
         ('7.19', 1.00E-04, 1.157E-10, 0.74, 0.8, 6, 'test_compare_sellmeijer/HeightAquiferD10L30.gid'),
         ('7.20', 1.00E-04, 1.157E-10, 1.37, 1.6, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60.gid'),
         ('7.21', 1.00E-04, 1.157E-10, 1.98, 2.4, 12, 'test_compare_sellmeijer/HeightAquiferD10L90.gid'),
         ('7.22', 1.00E-04, 1.157E-10, 0.65, 0.7, 9, 'test_compare_sellmeijer/HeightAquiferD20L30.gid'),
         ('7.23', 1.00E-04, 1.157E-10, 1.17, 1.3, 13.5, 'test_compare_sellmeijer/HeightAquiferD20L60.gid'),
         ('7.24', 1.00E-04, 1.157E-10, 1.68, 1.9, 15, 'test_compare_sellmeijer/HeightAquiferD20L90.gid'),
         ('7.25', 1.00E-04, 1.157E-10, 0.56, 0.6, 6, 'test_compare_sellmeijer/HeightAquiferD30L30.gid'),
         ('7.26', 1.00E-04, 1.157E-10, 1.08, 1.2, 16.5, 'test_compare_sellmeijer/HeightAquiferD30L60.gid'),
         ('7.27', 1.00E-04, 1.157E-10, 1.54, 1.7, 18, 'test_compare_sellmeijer/HeightAquiferD30L90.gid'),
         ('7.28', 3.00E-04, 1.157E-10, 2.22, 2.4, 6, 'test_compare_sellmeijer/HeightAquiferD10L30.gid'),
         ('7.29', 3.00E-04, 1.157E-10, 4.12, 4.8, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60.gid'),
         ('7.30', 3.00E-04, 1.157E-10, 5.93, 7.2, 12, 'test_compare_sellmeijer/HeightAquiferD10L90.gid'),
         ('7.31', 3.00E-04, 1.157E-10, 1.94, 2.1, 9, 'test_compare_sellmeijer/HeightAquiferD20L30.gid'),
         ('7.32', 3.00E-04, 1.157E-10, 3.52, 4, 13.5, 'test_compare_sellmeijer/HeightAquiferD20L60.gid'),
         ('7.33', 3.00E-04, 1.157E-10, 5.05, 5.9, 21, 'test_compare_sellmeijer/HeightAquiferD20L90.gid'),
         ('7.34', 3.00E-04, 1.157E-10, 1.68, 2, 9, 'test_compare_sellmeijer/HeightAquiferD30L30.gid'),
         ('7.35', 3.00E-04, 1.157E-10, 3.24, 3.7, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60.gid'),
         ('7.36', 3.00E-04, 1.157E-10, 4.61, 5.3, 27, 'test_compare_sellmeijer/HeightAquiferD30L90.gid')])
    def test_sellmeijers_rule_height(self, name, d70, kappa, Hc, Hn, length_n, test_name_gid):
        file_path = test_helper.get_file_path(os.path.join('./', test_name_gid))
        os.chdir(file_path)
        critical_head_found, length = self.critical_head_loop(file_path, kappa, d70, Hn, 'linear')
        self.results = {"value_name": name, "test_result_h": Hc, "kratos_results_h": critical_head_found,
                        "equivalent_software_h": Hn, "kratos_results_l": length, "equivalent_software_l": length_n}
        self.assertAlmostEqual(Hn, critical_head_found, 1,
                               f"Critical head kratos: {critical_head_found}, old geo flow {Hn}")

    @parameterized.expand(
        [('7.1', 1.00E-04, 1.157E-12, 3.43, 3.7, 6, 'test_compare_sellmeijer/HeightAquiferD10L30line'),
         ('7.2', 1.00E-04, 1.157E-12, 6.37, 7.4, 12, 'test_compare_sellmeijer/HeightAquiferD10L60line'),
         ('7.3', 1.00E-04, 1.157E-12, 9.18, 11.2, 13.5, 'test_compare_sellmeijer/HeightAquiferD10L90line'),
         ('7.4', 1.00E-04, 1.157E-12, 3, 3.3, 9, 'test_compare_sellmeijer/HeightAquiferD20L30line'),
         ('7.5', 1.00E-04, 1.157E-12, 5.44, 6.1, 15, 'test_compare_sellmeijer/HeightAquiferD20L60line'),
         ('7.6', 1.00E-04, 1.157E-12, 7.81, 9.1, 21, 'test_compare_sellmeijer/HeightAquiferD20L90line'),
         ('7.7', 1.00E-04, 1.157E-12, 2.6, 3.1, 10.5, 'test_compare_sellmeijer/HeightAquiferD30L30line'),
         ('7.8', 1.00E-04, 1.157E-12, 5.02, 5.7, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60line'),
         ('7.9', 1.00E-04, 1.157E-12, 7.13, 8.1, 24, 'test_compare_sellmeijer/HeightAquiferD30L90line'),
         ('7.10', 3.00E-04, 1.157E-12, 10.29, 11.3, 7.5, 'test_compare_sellmeijer/HeightAquiferD10L30line'),
         ('7.11', 3.00E-04, 1.157E-12, 19.1, 21.9, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60line'),
         ('7.12', 3.00E-04, 1.157E-12, 27.54, 32.8, 15, 'test_compare_sellmeijer/HeightAquiferD10L90line'),
         ('7.13', 3.00E-04, 1.157E-12, 9.01, 9.9, 10.5, 'test_compare_sellmeijer/HeightAquiferD20L30line'),
         ('7.14', 3.00E-04, 1.157E-12, 16.33, 18.2, 16.5, 'test_compare_sellmeijer/HeightAquiferD20L60line'),
         ('7.15', 3.00E-04, 1.157E-12, 23.42, 27, 21, 'test_compare_sellmeijer/HeightAquiferD20L90line'),
         ('7.16', 3.00E-04, 1.157E-12, 7.8, 9.4, 12, 'test_compare_sellmeijer/HeightAquiferD30L30line'),
         ('7.17', 3.00E-04, 1.157E-12, 15.05, 16.8, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60line'),
         ('7.18', 3.00E-04, 1.157E-12, 21.4, 24.1, 28.5, 'test_compare_sellmeijer/HeightAquiferD30L90line'),
         ('7.19', 1.00E-04, 1.157E-10, 0.74, 0.8, 6, 'test_compare_sellmeijer/HeightAquiferD10L30line'),
         ('7.20', 1.00E-04, 1.157E-10, 1.37, 1.6, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60line'),
         ('7.21', 1.00E-04, 1.157E-10, 1.98, 2.4, 12, 'test_compare_sellmeijer/HeightAquiferD10L90line'),
         ('7.22', 1.00E-04, 1.157E-10, 0.65, 0.7, 9, 'test_compare_sellmeijer/HeightAquiferD20L30line'),
         ('7.23', 1.00E-04, 1.157E-10, 1.17, 1.3, 13.5, 'test_compare_sellmeijer/HeightAquiferD20L60line'),
         ('7.24', 1.00E-04, 1.157E-10, 1.68, 1.9, 15, 'test_compare_sellmeijer/HeightAquiferD20L90line'),
         ('7.25', 1.00E-04, 1.157E-10, 0.56, 0.6, 6, 'test_compare_sellmeijer/HeightAquiferD30L30line'),
         ('7.26', 1.00E-04, 1.157E-10, 1.08, 1.2, 16.5, 'test_compare_sellmeijer/HeightAquiferD30L60line'),
         ('7.27', 1.00E-04, 1.157E-10, 1.54, 1.7, 18, 'test_compare_sellmeijer/HeightAquiferD30L90line'),
         ('7.28', 3.00E-04, 1.157E-10, 2.22, 2.4, 6, 'test_compare_sellmeijer/HeightAquiferD10L30line'),
         ('7.29', 3.00E-04, 1.157E-10, 4.12, 4.8, 10.5, 'test_compare_sellmeijer/HeightAquiferD10L60line'),
         ('7.30', 3.00E-04, 1.157E-10, 5.93, 7.2, 12, 'test_compare_sellmeijer/HeightAquiferD10L90line'),
         ('7.31', 3.00E-04, 1.157E-10, 1.94, 2.1, 9, 'test_compare_sellmeijer/HeightAquiferD20L30line'),
         ('7.32', 3.00E-04, 1.157E-10, 3.52, 4, 13.5, 'test_compare_sellmeijer/HeightAquiferD20L60line'),
         ('7.33', 3.00E-04, 1.157E-10, 5.05, 5.9, 21, 'test_compare_sellmeijer/HeightAquiferD20L90line'),
         ('7.34', 3.00E-04, 1.157E-10, 1.68, 2, 9, 'test_compare_sellmeijer/HeightAquiferD30L30line'),
         ('7.35', 3.00E-04, 1.157E-10, 3.24, 3.7, 22.5, 'test_compare_sellmeijer/HeightAquiferD30L60line'),
         ('7.36', 3.00E-04, 1.157E-10, 4.61, 5.3, 27, 'test_compare_sellmeijer/HeightAquiferD30L90line')])
    def test_sellmeijers_rule_height_line(self, name, d70, kappa, Hc, Hn, length_n, test_name_line):
        file_path = test_helper.get_file_path(os.path.join('./', test_name_line))
        os.chdir(file_path)
        critical_head_found, length = self.critical_head_loop(file_path, kappa, d70, Hn, 'linear')
        self.results = {"value_name": name, "test_result_h": Hc, "kratos_results_h": critical_head_found,
                        "equivalent_software_h": Hn, "kratos_results_l": length, "equivalent_software_l": length_n}
        self.assertAlmostEqual(Hn, critical_head_found, 1,
                               f"Critical head kratos: {critical_head_found}, old geo flow {Hn}")
