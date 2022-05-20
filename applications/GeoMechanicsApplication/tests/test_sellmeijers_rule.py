import sys
import os
import csv
import json
import math
sys.path.append(os.path.join('D:/kratos'))
import KratosMultiphysics.KratosUnittest as KratosUnittest

sys.path.append(os.path.join('..', 'python_scripts'))
import test_helper


class LatexWriterFile:
    """
    Class that writes results for latex documentation
    """

    def __init__(self, filename=""):
        self.filename = filename
        self.value_dict_default = {"value_name": "", "test_result": 0, "kratos_results": 0}

    def write_latex_file(self, result_list):
        with open(self.filename, "w+") as output_latex_file:
            for result_pair in result_list:
                if result_pair['kratos_results'] is None:
                    error = 1
                    error_equivalent_software = 1
                    output_latex_file.write(
                        f"{str(result_pair['value_name'])} & {result_pair['test_result']} & "
                        f"{result_pair['equivalent_software']} & "
                        f" NaN & {round(error, 2)} &"
                        f" {round(error_equivalent_software, 2)} \\\\ \hline \n")
                else:
                    error = abs(result_pair['kratos_results'] - result_pair['test_result']) / \
                            (abs(result_pair['test_result']) + 1e-60)
                    error_equivalent_software = abs(result_pair['kratos_results'] - result_pair['equivalent_software']) / \
                                                (abs(result_pair['equivalent_software']) + 1e-60)
                    output_latex_file.write(
                        f"{result_pair['value_name']} & {result_pair['test_result']} & "
                        f"{result_pair['equivalent_software']} & "
                        f" {round(result_pair['kratos_results'], 2)} & {round(error, 2)} &"
                        f" {round(error_equivalent_software, 2)} \\\\ \hline \n")


class TestSellmeijersRule(KratosUnittest.TestCase):

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        self.latex_writer = LatexWriterFile()
        self.test_lists = self.csv_file_reader()
        self.gid_files = {10: {30: "test_compare_sellmeijer/HeightAquiferD10L30.gid",
                               60: "test_compare_sellmeijer/HeightAquiferD10L60.gid",
                               90: "test_compare_sellmeijer/HeightAquiferD10L90.gid"},
                          20: {30: "test_compare_sellmeijer/HeightAquiferD20L30.gid",
                               60: "test_compare_sellmeijer/HeightAquiferD20L60.gid",
                               90: "test_compare_sellmeijer/HeightAquiferD20L90.gid"},
                          30: {30: "test_compare_sellmeijer/HeightAquiferD30L30.gid",
                               60: "test_compare_sellmeijer/HeightAquiferD30L60.gid",
                               90: "test_compare_sellmeijer/HeightAquiferD30L90.gid"}}
        self.is_running_under_teamcity = test_helper.is_running_under_teamcity()

    def tearDown(self):
        pass


    def csv_file_reader(self):
        with open(test_helper.get_file_path(os.path.join('.', 'test_compare_sellmeijer/tests.csv')), 'r') as file:
            results = {"name": [], "L": [], "D": [], "d70": [], "kappa": [], "Hc": [], "Hn": [], "Hc_kratos": []}
            reader = csv.reader(file, delimiter=';')
            for counter, row in enumerate(reader):
                if counter != 0:
                    results["name"].append(row[0])
                    results["L"].append(float(row[1]))
                    results["D"].append(float(row[2]))
                    results["d70"].append(float(row[3]))
                    results["kappa"].append(float(row[4]))
                    results["Hc"].append(float(row[5]))
                    results["Hn"].append(float(row[6]))
        return results

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
                parameters['processes']['constraints_process_list'][0]['Parameters']['reference_coordinate'] = head_level
            else:
                parameters['processes']['constraints_process_list'][1]['Parameters']['reference_coordinate'] = head_level
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
        return None

    def critical_head_loop(self, file_path, counter, search_type='linear'):
        self.change_material_parameters(file_path, self.test_lists["kappa"][counter], self.test_lists["d70"][counter])
        heads = [x * 0.1 for x in range(int(self.test_lists["Hc"][counter] * 10 - 40), int(self.test_lists["Hc"][counter] * 10 + 90), 1)]
        critical_head_found = math.nan
        if search_type == 'linear':
            critical_head_found = self.linear_search(file_path, heads)
        self.test_lists["Hc_kratos"].append(critical_head_found)

    def test_sellmeijers_rule_height(self):
        for counter, test_name in enumerate(self.test_lists["name"]):
            test_name_gid = self.gid_files[self.test_lists["D"][counter]][self.test_lists["L"][counter]]
            file_path = test_helper.get_file_path(os.path.join('./', test_name_gid))
            os.chdir(file_path)
            self.critical_head_loop(file_path, counter, 'linear')
        all_results = []
        for counter, test_n in enumerate(self.test_lists['name']):
            index_test = self.test_lists['name'].index(test_n)
            temp_results = {"value_name": test_n,
                            "test_result": self.test_lists['Hc'][index_test],
                            "equivalent_software": self.test_lists['Hn'][index_test],
                            "kratos_results": self.test_lists['Hc_kratos'][counter]}
            all_results.append(temp_results)
        if self.is_running_under_teamcity:
            self.latex_writer.filename = test_helper.get_file_path('test_compare_sellmeijer/test_compare_sellmeijer.tex')
            self.latex_writer.write_latex_file(all_results)
