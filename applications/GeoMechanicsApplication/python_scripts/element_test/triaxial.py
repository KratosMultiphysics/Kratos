import sys
import os
import json
import numpy as np

import KratosMultiphysics as Kratos
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis import GeoMechanicsAnalysis

tests_folder_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'tests'))
if os.path.exists(tests_folder_path):
    sys.path.append(tests_folder_path)
    import test_helper
else:
    raise ImportError(f"Tests folder not found at: {tests_folder_path}")


class TriaxialTest:
    def __init__(self, json_file_path):
        self.json_file_path = json_file_path
        self.data = self._read_json()

    def _read_json(self):
        with open(self.json_file_path, 'r') as file:
            return json.load(file)

    def _write_json(self):
        with open(self.json_file_path, 'w') as file:
            json.dump(self.data, file, indent=4)

    def modify_umat_parameters(self, new_value):
        self.data['properties'][0]['Material']['Variables']['UMAT_PARAMETERS'][0] = new_value
        self._write_json()

    def read_umat_parameters(self):
        try:
            umat_parameters = self.data['properties'][0]['Material']['Variables']['UMAT_PARAMETERS']
            cohesion = umat_parameters[2]
            friction_angle = umat_parameters[3]
            return cohesion, friction_angle
        except KeyError:
            return 0.0, 30.0  #TODO: Remove the reasonable defaults for Mohr-Coulomb failure envelope

class TriaxialTestRunner:
    def __init__(self, output_file_paths, work_dir):
        self.output_file_paths = output_file_paths
        self.work_dir = work_dir

    def run(self):
        parameters = self._load_stage_parameters()
        self._execute_analysis_stages(parameters)

        stress, mean_stress, von_mises, _, strain = self._collect_results()
        tensors = self._extract_stress_tensors(stress)
        yy_strain, vol_strain = self._compute_strains(strain)
        von_mises_values = self._compute_scalar_stresses(von_mises)
        mean_stress_values = self._compute_scalar_stresses(mean_stress)

        return tensors, yy_strain, vol_strain, von_mises_values, mean_stress_values

    def _load_stage_parameters(self):
        parameter_files = [os.path.join(self.work_dir, f'ProjectParameters.json')]
        parameters = []
        for f in parameter_files:
            with open(f, 'r') as file:
                parameters.append(Kratos.Parameters(file.read()))
        return parameters


    def _execute_analysis_stages(self, parameters):
        model = Kratos.Model()
        stages = [GeoMechanicsAnalysis(model, p) for p in parameters]
        original_cwd = os.getcwd()
        try:
            os.chdir(self.work_dir)
            for stage in stages:
                stage.Run()
        finally:
            os.chdir(original_cwd)

    def _collect_results(self):
        stress, mean_stress, von_mises, displacement, strain = [], [], [], [], []

        for path in self.output_file_paths:
            output = test_helper.GiDOutputFileReader().read_output_from(path)
            for result_name, items in output["results"].items():
                for item in items:
                    self._categorize_result(result_name, item, stress, mean_stress, von_mises, displacement, strain)

        return stress, mean_stress, von_mises, displacement, strain

    def _categorize_result(self, result_name, item, stress, mean_stress, von_mises, displacement, strain):
        values = item["values"]
        if result_name == "CAUCHY_STRESS_TENSOR":
            stress.append(item)
        elif result_name == "MEAN_EFFECTIVE_STRESS" and self._is_tri3(values):
            mean_stress.append(item)
        elif result_name == "VON_MISES_STRESS" and self._is_tri3(values):
            von_mises.append(item)
        elif result_name == "DISPLACEMENT":
            displacement.append(item)
        elif result_name == "ENGINEERING_STRAIN_TENSOR":
            strain.append(item)

    def _is_tri3(self, values):
        return isinstance(values, list) and all("value" in v and isinstance(v["value"], list) and len(v["value"]) == 3 for v in values)

    def _extract_stress_tensors(self, stress_results):
        reshaped = {}
        for result in stress_results:
            time_step = result["time"]
            values = result["values"]
            if not values:
                continue
            sublist = values[0]["value"][0]
            tensor = np.array([
                [sublist[0], sublist[3], sublist[5]],
                [sublist[3], sublist[1], sublist[4]],
                [sublist[5], sublist[4], sublist[2]],
            ])
            reshaped[time_step] = [tensor]
        return reshaped

    def _compute_strains(self, strain_results):
        yy, vol = [], []
        for result in strain_results:
            values = result["values"]
            if not values:
                continue
            eps_xx, eps_yy, eps_zz = values[0]["value"][0][:3]
            vol.append(eps_xx + eps_yy + eps_zz)
            yy.append(eps_yy)
        return yy, vol

    def _compute_scalar_stresses(self, results):
        return [r["values"][0]["value"][1] for r in results if r["values"]]
