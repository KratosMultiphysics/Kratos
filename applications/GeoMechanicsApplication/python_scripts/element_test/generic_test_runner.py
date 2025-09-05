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


class GenericTestRunner:
    def __init__(self, output_file_paths, work_dir):
        self.output_file_paths = output_file_paths
        self.work_dir = work_dir

    def run(self):
        parameters = self._load_stage_parameters()
        self._execute_analysis_stages(parameters)

        stress, mean_stress, von_mises, _, strain = self._collect_results()
        tensors = self._extract_stress_tensors(stress)
        shear_stress_xy = self._extract_shear_stress_xy(stress)
        yy_strain, vol_strain, shear_strain_xy = self._compute_strains(strain)
        von_mises_values = self._compute_scalar_stresses(von_mises)
        mean_stress_values = self._compute_scalar_stresses(mean_stress)

        return tensors, yy_strain, vol_strain, von_mises_values, mean_stress_values, shear_stress_xy, shear_strain_xy

    def _load_stage_parameters(self):
        parameter_files = [os.path.join(self.work_dir, 'ProjectParameters.json')]
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
        elif result_name == "MEAN_EFFECTIVE_STRESS" and self._is_tri3_element_gp(values):
            mean_stress.append(item)
        elif result_name == "VON_MISES_STRESS" and self._is_tri3_element_gp(values):
            von_mises.append(item)
        elif result_name == "DISPLACEMENT":
            displacement.append(item)
        elif result_name == "ENGINEERING_STRAIN_TENSOR":
            strain.append(item)

    def _is_tri3_element_gp(self, values):
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

    def _extract_shear_stress_xy(self, stress_results):
        shear_stress_xy = []
        for result in stress_results:
            values = result["values"]
            if not values:
                continue
            stress_components = values[0]["value"][0]
            shear_xy = stress_components[3]
            shear_stress_xy.append(shear_xy)
        return shear_stress_xy

    def _compute_strains(self, strain_results):
        yy, vol, shear_xy = [], [], []
        for result in strain_results:
            values = result["values"]
            if not values:
                continue
            eps_xx, eps_yy, eps_zz, eps_xy = values[0]["value"][0][:4]
            vol.append(eps_xx + eps_yy + eps_zz)
            yy.append(eps_yy)
            shear_xy.append(eps_xy)
        return yy, vol, shear_xy

    def _compute_scalar_stresses(self, results):
        return [r["values"][0]["value"][1] for r in results if r["values"]]
