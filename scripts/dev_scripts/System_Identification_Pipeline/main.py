import os
import shutil
from glob import glob
import math

import numpy as np

# from memory_profiler import profile


import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from Measurement_file_creation.MeasurementFileStructure import *
from Measurement_file_creation.MeasurementFromSimulationDataGenerator import MeasurementFromSimulationDataGenerator
from Measurement_file_creation.PerModelPartMaterialChanger import PerModelPartMaterialChanger

from Measurement_file_creation.EffectiveSensorPositionsCalculator import EffectiveSensorPositionsCalculator

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis


class SimulationSetup():

    def __init__(self) -> None:
        self.simulation_folders = []
        self.sensor_position_configurations_list = []
        self.material_parameters_to_change = {}
        self.measurement_load_cases_container = LoadCasesContainer()
        self.sensor_list = []

    def measurement_creation(self) -> None:

        data_generator = MeasurementFromSimulationDataGenerator()

        self.measurement_load_cases_container = data_generator.write_measurement_data_file(
            mdpa_model_file_name="model_file",
            primal_parameter_files=["primal_parameters.json"],
            # primal_parameter_files=["primal_parameters.json", "primal_parameters1.json", "primal_parameters2.json", "primal_parameters3.json"],
            sensor_positions=self.sensor_list,
            output_file_name="measurement_data.json",
            material_changer=self.material_changer)

    def get_additional_folder_name_from_measurement_data(self) -> str:
        n_load_cases = len(self.measurement_load_cases_container.load_cases)
        n_sensors = len(self.measurement_load_cases_container.load_cases[0].sensors_infos)
        data_name = f"_lc{n_load_cases}_s{n_sensors}"
        return data_name

    def cascading_copy_simulation_files(self, working_folder: str) -> None:

        working_folder = "./"+working_folder
        sub_folder_list = working_folder.split(sep="/")
        cumulative_folder_path = ""
        for i, folder in enumerate(sub_folder_list[0:-1]):
            if i == 0:
                cumulative_folder_path += f"{folder}/"
            else:
                cumulative_folder_path += f"/{folder}"
            for file in glob(f"{cumulative_folder_path}/*_parameters.json")+glob(f"{cumulative_folder_path}/*_parameters*[0-9].json")+glob(f"{cumulative_folder_path}/StructuralMaterials.json")+glob(f"{cumulative_folder_path}/StructuralMaterials*[0-9].json"):
                folders_in_file_path = str(file).split("/")
                shutil.copy(str(file), f"{working_folder}/"+folders_in_file_path[-1])

    def pre_operations(self, folder_suffix: str = "") -> None:
        self.measurement_creation()
        data_name = self.get_additional_folder_name_from_measurement_data()
        try:
            shutil.rmtree(f"./Optimization_Results")
        except:
            pass
        try:
            shutil.rmtree(f"./other_results")
        except:
            pass
        try:
            shutil.rmtree(f"./primal_vtk_output")
        except:
            pass
        try:
            shutil.rmtree(f"./Optimization_Results{data_name}_{folder_suffix}")
        except:
            pass
        try:
            shutil.rmtree(f"./other_results{data_name}_{folder_suffix}")
        except:
            pass
        try:
            shutil.rmtree(f"./primal_vtk_output{data_name}_{folder_suffix}")
        except:
            pass

    def post_operations(self, folder_suffix: str = "") -> None:
        data_name = self.get_additional_folder_name_from_measurement_data()
        data_name += f"_{folder_suffix}"

        # remove primal analysis files
        for file in glob('./Structure*.h5')+glob('./all_nodes_elements_model_part*.h5')+glob('./results_primal*.h5'):
            os.remove(file)

        for file in glob('./lc*.vtu')+glob('./measurement*.json')+glob('./*.csv')+glob('./*.time')+glob('./*.html'):
            os.rename(str(file), f"./Optimization_Results"+str(file)[1:])

        try:
            os.rename("./Optimization_Results", f"./Optimization_Results{data_name}")
        except:
            pass
        try:
            os.rename("./primal_vtk_output", f"./primal_vtk_output{data_name}")
        except:
            pass

    # @profile
    def run_simulations(self) -> None:
        with open("optimization_parameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())

        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)

        # try:
        analysis.Run()
        # except:
        #     print(f"Error in analysis: {working_folder}")

    def run(self) -> None:

        self.simulation_folders = [
            # "benchmark1/case1a/0_coar",
            "benchmark1/case1b/0_coar",
            # "benchmark1/case1c/0_coar",
            # "benchmark1/case1d/0_coar",
            # "benchmark1/case1a/1_medium",
            # "benchmark1/case1b/1_medium",
            # "benchmark1/case1c/1_medium", #
            # "benchmark1/case1d/1_medium",
            # "benchmark1/case1a/2_fine", #
            # "benchmark1/case1b/2_fine", #
            # "benchmark1/case1c/2_fine", #
            # "benchmark1/case1d/2_fine",

            # "benchmark2/case2a/0_coar",
            # "benchmark2/case2b/0_coar",
            # "benchmark2/case2c/0_coar",
            # "benchmark2/case2a/1_medium",
            # "benchmark2/case2b/1_medium",
            # "benchmark2/case2c/1_medium",

            # "benchmark3/case3a/0_coar",
            # "benchmark3/case3b/0_coar",
            # "benchmark3/case3c/0_coar",
            # "benchmark3/case3a/1_medium",
            # "benchmark3/case3b/1_medium",
            # "benchmark3/case3c/1_medium",
        ]

        xy_grid = []

        # # For BM 1
        # for x in np.linspace(0.0, 5.0, 20):
        #     for y in np.linspace(0.0, 0.5, 4):
        #         xy_grid.append([x, y, 0,"X"])

        # # For BM 2
        # for x in np.linspace(0,60,20):  # (10,50,5): #(5, 55, 10):
        #     for y in np.linspace(0,30,10):  # (5,25,3): #(5, 25, 5):
        #         # Exclude circle
        #         off_set_x = x - 30
        #         off_set_y = y - 15
        #         r = math.sqrt(off_set_x**2 + off_set_y**2)
        #         if r <= 5.0:
        #             continue

        #         xy_grid.append([x, y, 0,"X"])

        # For BM 3
        # for x in np.linspace(0, 14, 10):
        #     for y in np.linspace(0, 14, 10):
        #         xy_grid.append([x, y, 0,"X"])

        self.sensor_position_configurations_list = [

            # Suneth coar auto sensor
            # [
            #     [2.30683, 0.418347, 0.0, "Y"],
            #     [2.62784, 0.0, 0.0, "Y"],
            #     [2.86172, 0.5, 0.0, "Y"],
            #     [1.83304, 0.0864105, 0.0, "Y"],
            #     [3.24836, 0.167399, 0.0, "Y"],
            #     [1.53146, 0.5, 0.0, "Y"],
            #     [3.65, 0.5, 0.0, "Y"],
            #     [3.75, 0.0, 0.0, "Y"],
            #     [1.28352, 0.0, 0.0, "Y"],
            #     [1.02724, 0.452161, 0.0, "Y"],
            #     [4.17229, 0.329207, 0.0, "Y"],
            #     [4.55001, 0.0, 0.0, "Y"],
            #     [4.7, 0.5, 0.0, "Y"],
            #     [0.783517, 0.0, 0.0, "Y"],
            #     [0.419103, 0.423801, 0.0, "Y"],
            #     [0.0, 0.140249, 0.0, "Y"],
            # ],


            # Suneth medium auto sensor
            # [
            #     [2.34969, 0.412756, 0.0, "Y"],
            #     [2.65051, 0.0, 0.0, "Y"],
            #     [2.85072, 0.5, 0.0, "Y"],
            #     [2.00101, 0.0, 0.0, "Y"],
            #     [3.15001, 0.0866035, 0.0, "Y"],
            #     [1.84454, 0.5, 0.0, "Y"],
            #     [3.46523, 0.5, 0.0, "Y"],
            #     [1.49977, 0.0882805, 0.0, "Y"],
            #     [3.65, 0.0, 0.0, "Y"],
            #     [1.2, 0.5, 0.0, "Y"],
            #     [3.97249, 0.5, 0.0, "Y"],
            #     [4.15001, 0.0, 0.0, "Y"],
            #     [1.001, 0.0, 0.0, "Y"],
            #     [4.49597, 0.388373, 0.0, "Y"],
            #     [0.700041, 0.431113, 0.0, "Y"],
            #     [0.401447, 0.0, 0.0, "Y"],
            #     [5.0, 0.452691, 0.0, "Y"],
            #     [0.194539, 0.5, 0.0, "Y"],
            # ],

            # xy_grid,

            # [[0.5, 0.5, 0.0,"X"],
            #  [0.75, 0.5, 0.0,"X"],
            #  [1.0, 0.5, 0.0,"X"],
            #  [1.25, 0.5, 0.0,"X"],
            #  [1.5, 0.5, 0.0,"X"],
            #  [1.75, 0.5, 0.0,"X"],
            #  [2.0, 0.5, 0.0,"X"],
            #  [2.25, 0.5, 0.0,"X"],
            #  [2.5, 0.5, 0.0,"X"],
            #  [2.75, 0.5, 0.0,"X"],
            #  [3.0, 0.5, 0.0,"X"],
            #  [3.25, 0.5, 0.0,"X"],
            #  [3.5, 0.5, 0.0,"X"],
            #  [3.75, 0.5, 0.0,"X"],
            #  [4.0, 0.5, 0.0,"X"],
            #  [4.25, 0.5, 0.0,"X"],
            #  [4.5, 0.5, 0.0,"X"],
            #  [0.5, 0.0, 0.0,"X"],
            #  [0.75, 0.0, 0.0,"X"],
            #  [1.0, 0.0, 0.0,"X"],
            #  [1.25, 0.0, 0.0,"X"],
            #  [1.5, 0.0, 0.0,"X"],
            #  [1.75, 0.0, 0.0,"X"],
            #  [2.0, 0.0, 0.0,"X"],
            #  [2.25, 0.0, 0.0,"X"],
            #  [2.5, 0.0, 0.0,"X"],
            #  [2.75, 0.0, 0.0,"X"],
            #  [3.0, 0.0, 0.0,"X"],
            #  [3.25, 0.0, 0.0,"X"],
            #  [3.5, 0.0, 0.0,"X"],
            #  [3.75, 0.0, 0.0,"X"],
            #  [4.0, 0.0, 0.0,"X"],
            #  [4.25, 0.0, 0.0,"X"],
            #  [4.5, 0.0, 0.0,"X"]],

            # [[0.5, 0.5, 0.0,"X"],
            #  [0.75, 0.5, 0.0,"X"],
            #  [1.0, 0.5, 0.0,"X"],
            #  [1.25, 0.5, 0.0,"X"],
            #  [1.5, 0.5, 0.0,"X"],
            #  [1.75, 0.5, 0.0,"X"],
            #  [2.0, 0.5, 0.0,"X"],
            #  [2.25, 0.5, 0.0,"X"],
            #  [2.5, 0.5, 0.0,"X"],
            #  [2.75, 0.5, 0.0,"X"],
            #  [3.0, 0.5, 0.0,"X"],
            #  [3.25, 0.5, 0.0,"X"],
            #  [3.5, 0.5, 0.0,"X"],
            #  [3.75, 0.5, 0.0,"X"],
            #  [4.0, 0.5, 0.0,"X"],
            #  [4.25, 0.5, 0.0,"X"],
            #  [4.5, 0.5, 0.0,"X"]],

            # [[0.5, 0.0, 0.0,"X"],
            #  [1.0, 0.0, 0.0,"X"],
            #  [1.5, 0.0, 0.0,"X"],
            #  [2.0, 0.0, 0.0,"X"],
            #  [2.5, 0.0, 0.0,"X"],
            #  [3.0, 0.0, 0.0,"X"],
            #  [3.5, 0.0, 0.0,"X"],
            #  [4.0, 0.0, 0.0,"X"],
            #  [4.5, 0.0, 0.0,"X"]],

            # [[0.5, 0.0, 0.0,"X"],
            #  [1.0, 0.0, 0.0,"X"],
            #  [1.5, 0.0, 0.0,"X"],
            #  [2.0, 0.0, 0.0,"X"],
            #  [2.5, 0.0, 0.0,"X"],
            #  [3.0, 0.0, 0.0,"X"],
            #  [3.5, 0.0, 0.0,"X"],
            #  [4.0, 0.0, 0.0,"X"],
            #  [4.5, 0.0, 0.0,"X"],
            #  [0.5, 0.5, 0.0,"X"],
            #  [1.0, 0.5, 0.0,"X"],
            #  [1.5, 0.5, 0.0,"X"],
            #  [2.0, 0.5, 0.0,"X"],
            #  [2.5, 0.5, 0.0,"X"],
            #  [3.0, 0.5, 0.0,"X"],
            #  [3.5, 0.5, 0.0,"X"],
            #  [4.0, 0.5, 0.0,"X"],
            #  [4.5, 0.5, 0.0,"X"]],

            # [[0.5, 0.5, 0.0,"X"],
            #  [1.5, 0.5, 0.0,"X"],
            #  [2.5, 0.5, 0.0,"X"],
            #  [3.5, 0.5, 0.0,"X"],
            #  [4.5, 0.5, 0.0,"X"]],

            ##################

            # suneth coar auto sensor
            # [
            #     [60.0,  30.0,       0.0, "X"],
            #     [60.0,   0.0,       0.0, "Y"],
            #     [60.0,  15.0001,       0.0, "Y"],
            #     [60.0,   7.50003,       0.0, "Y"],
            #     [60.0,  22.5001,       0.0, "Y"],
            #     [54.748,   4.56132,       0.0, "Y"],
            #     [54.0231,  27.3096,       0.0, "Y"],
            #     [53.2274,  20.0081,       0.0, "Y"],
            #     [53.4303,  11.2092,       0.0, "Y"],
            #     [48.2662,  15.5791,       0.0, "Y"],
            #     [46.9118,   9.44347,       0.0, "Y"],
            #     [47.2241,  22.0901,       0.0, "Y"],
            #     [33.379,  25.6485,       0.0, "Y"],
            #     [30.0,  10.0,       0.0, "Y"],
            #     [26.805,  18.8458,       0.0, "Y"],
            #     [24.8895,  28.3791,       0.0, "Y"],
            #     [46.7033,   1.89605,       0.0, "X"],
            #     [33.8458,  18.195,       0.0, "X"],
            #     [39.6996,  13.9511,       0.0, "X"],
            #     [37.5,   0.0,       0.0, "X"],
            #     [27.5,   0.0,       0.0, "X"],
            #     [36.0218,   6.78236,       0.0, "X"],
            #     [12.5,  30.0,       0.0, "X"],
            #     [5.00002,   0.0,       0.0, "X"],
            # ],

            # suneth medium auto sensor
            # [
            #     [60.0,  30.0,       0.0, "Y"],
            #     [60.0,   2.50001,       0.0, "Y"],
            #     [60.0,  23.7501,       0.0, "Y"],
            #     [60.0,   8.75004,       0.0, "Y"],
            #     [60.0,  15.0001,       0.0, "Y"],
            #     [54.6844,  26.8897,       0.0, "Y"],
            #     [54.3573,   5.41175,       0.0, "Y"],
            #     [53.9305,  20.687,       0.0, "Y"],
            #     [51.25,   0.0,       0.0, "Y"],
            #     [54.2268,  11.9216,       0.0, "Y"],
            #     [47.2608,  15.8356,       0.0, "Y"],
            #     [47.4013,   9.20476,       0.0, "Y"],
            #     [48.7499,  30.0,       0.0, "Y"],
            #     [45.0,   1.25001,       0.0, "Y"],
            #     [43.0007,  22.7291,       0.0, "Y"],
            #     [30.018,  21.4715,       0.0, "Y"],
            #     [27.4999,  30.0,       0.0, "Y"],
            #     [27.4407,   6.44933,       0.0, "Y"],
            #     [26.25,   0.0,       0.0, "Y"],
            #     [20.0,  26.25,       0.0, "Y"],
            #     [18.7089,   4.15801,       0.0, "Y"],
            #     [19.3087,  18.6327,       0.0, "X"],
            #     [34.9254,  14.1383,       0.0, "X"],
            #     [36.2501,   0.0,       0.0, "X"],
            #     [35.0,  30.0,       0.0, "X"],
            #     [10.2961,  15.1852,       0.0, "X"],
            #     [11.25,   0.0,       0.0, "X"],
            #     [7.6532,  25.5873,       0.0, "X"],
            #     [2.23349,  28.6948,       0.0, "X"],
            # ],

            # xy_grid,

            # [[50, 25, 0,"X"],
            #  [50, 5, 0,"X"],
            #  [40, 25, 0,"X"],
            #  [40, 5, 0,"X"],
            #  [30, 25, 0,"X"],
            #  [30, 5, 0,"X"],
            #  [20, 25, 0,"X"],
            #  [20, 5, 0,"X"],
            #  [10, 25, 0,"X"],
            #  [10, 5, 0,"X"],
            #  [10, 15, 0,"X"],
            #  [20, 15, 0,"X"],
            #  [40, 15, 0,"X"],
            #  [50, 15, 0,"X"]],

            ##################

            # Suneth coar auto sensor
            # [
            #     [14.0, 14.0,       0.0, "Y"],
            #     [14.0, 1.0,       0.0, "Y"],
            #     [14.0, 11.0,       0.0, "Y"],
            #     [14.0, 4.00002,       0.0, "Y"],
            #     [14.0, 7.50003,       0.0, "Y"],
            #     [11.3685, 12.4534,       0.0, "Y"],
            #     [11.0, 9.46629e-11,       0.0, "Y"],
            #     [10.9825, 3.07101,       0.0, "Y"],
            #     [11.124, 9.30573,       0.0, "Y"],
            #     [10.3254, 6.3515,       0.0, "Y"],
            #     [8.49997, 14.0,       0.0, "Y"],
            #     [8.10208, 8.66092,       0.0, "Y"],
            #     [7.55054, 2.74574,       0.0, "Y"],
            #     [6.00001, 0.0,       0.0, "Y"],
            #     [6.19644, 11.8532,       0.0, "Y"],
            #     [3.12466, 11.3437,       0.0, "Y"],
            #     [6.75677, 5.71431,       0.0, "Y"],
            #     [1.50001, 0.0,       0.0, "X"],
            #     [3.72271, 2.08911,       0.0, "X"],
            #     [2.5, 5.0,       0.0, "X"],
            # ],

            # Suneth medium auto sensor
            # [
            #     [14.0, 0.0,       0.0, "X"],
            #     [14.0, 14.0,       0.0, "X"],
            #     [11.0, 14.0,       0.0, "X"],
            #     [14.0, 11.0,       0.0, "Y"],
            #     [14.0, 3.50001,       0.0, "Y"],
            #     [14.0, 7.75003,       0.0, "Y"],
            #     [11.5064, 1.73662,       0.0, "Y"],
            #     [11.7209, 5.7206,       0.0, "Y"],
            #     [10.9665, 10.7562,       0.0, "Y"],
            #     [9.00001, 1.57772e-10,       0.0, "Y"],
            #     [9.28823, 3.76963,       0.0, "Y"],
            #     [9.63662, 8.03884,       0.0, "Y"],
            #     [7.11965, 6.33218,       0.0, "Y"],
            #     [8.26957, 12.4788,       0.0, "Y"],
            #     [6.67378, 2.18829,       0.0, "Y"],
            #     [5.3282, 13.2268,       0.0, "Y"],
            #     [3.00427, 11.2863,       0.0, "Y"],
            #     [4.72512, 8.41281,       0.0, "X"],
            #     [1.50001, 0.0,       0.0, "Y"],
            #     [4.50002, 0.0,       0.0, "X"],
            #     [3.17417, 5.01326,       0.0, "X"],
            #     [0.86566, 7.00028,       0.0, "X"],
            #     [0.902568, 2.97293,       0.0, "X"],
            # ],

            # xy_grid,

            # [
            #  [2, 0, 0,"X"],
            #  [4.5, 0, 0,"X"],
            #  [7, 0, 0,"X"],
            #  [9.5, 0, 0,"X"],
            #  [12, 0, 0,"X"],

            #  [2, 2, 0,"X"],
            #  [4.5, 2, 0,"X"],
            #  [7, 2, 0,"X"],
            #  [9.5, 2, 0,"X"],
            #  [12, 2, 0,"X"],

            #  [2, 4.5, 0,"X"],
            #  [4.5, 4.5, 0,"X"],
            #  [7, 4.5, 0,"X"],
            #  [9.5, 4.5, 0,"X"],
            #  [12, 4.5, 0,"X"],

            #  [2, 7, 0,"X"],
            #  [4.5, 7, 0,"X"],
            #  [7, 7, 0,"X"],
            #  [9.5, 7, 0,"X"],
            #  [12, 7, 0,"X"],

            #  [2, 9.5, 0,"X"],
            #  [4.5, 9.5, 0,"X"],
            #  [7, 9.5, 0,"X"],
            #  [9.5, 9.5, 0,"X"],
            #  [12, 9.5, 0,"X"],

            #  [2, 12, 0,"X"],
            #  [4.5, 12, 0,"X"],
            #  [7, 12, 0,"X"],
            #  [9.5, 12, 0,"X"],
            #  [12, 12, 0,"X"],

            #  [2, 14, 0,"X"],
            #  [4.5, 14, 0,"X"],
            #  [7, 14, 0,"X"],
            #  [9.5, 14, 0,"X"],
            #  [12, 14, 0,"X"]]

        ]

        material_parameters_to_change = {
            "ElemMat1": {"YOUNG_MODULUS": 30000000000.0},   # Basis
            "ElemMat2": {"YOUNG_MODULUS": 10000000000.0},   # BM1 | BM2 | BM3(Stahl)
            # "ElemMat3": {"YOUNG_MODULUS": 210000000000.0},   # x   | BM2 | BM3(Stahl)
            # "ElemMat4": {"YOUNG_MODULUS": 210000000000.0},   # x   | BM2 | BM3(Stahl)
            # "ElemMat5": {"YOUNG_MODULUS": 210000000000.0},   # x   | BM2 | BM3(Stahl)
            # "ElemMat6": {"YOUNG_MODULUS": 10000000000.0},   # x   | BM2 | BM3
            # "ElemMat7": {"YOUNG_MODULUS": 10000000000.0},   # x   | x   | BM3
            # "ElemMat8": {"YOUNG_MODULUS": 10000000000.0},   # x   | x   | BM3
            # "ElemMat9": {"YOUNG_MODULUS": 10000000000.0},   # x   | x   | BM3
        }

        self.material_changer = PerModelPartMaterialChanger(material_parameters_to_change)

        ############################################################################

        oldPath = os.getcwd()

        for working_folder in self.simulation_folders:

            simulation_description = "smoothing3_sensX_autoSensor_test3"

            os.chdir(oldPath)

            self.cascading_copy_simulation_files(working_folder)

            os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

            sensor_placer = EffectiveSensorPositionsCalculator()
            auto_found_sensor_positions = sensor_placer.calculate_sensor_positions(model_file_name="model_file",
                                                                                    primal_settings_file_name="primal_parameters.json",
                                                                                    adjoint_settings_file_name="adjoint_parameters.json",
                                                                                    potential_positions_model_part="all_nodes_elements_model_part",
                                                                                    num_of_sensors_to_select=9,
                                                                                    sensor_direction=[1.0,0.0,0.0],
                                                                                    redundancy_factor=1.0)
            for pos in auto_found_sensor_positions: pos.append("X")
            self.sensor_position_configurations_list.append(auto_found_sensor_positions)

            for sensor_config in self.sensor_position_configurations_list:

                ############################################################################

                self.sensor_list.clear()
                for position in sensor_config:
                    if position[3] == "X":
                        self.sensor_list.append(
                            SensorDataContainer(position_of_mesh_node=[position[0], position[1], position[2]], type_of_sensor="DISPLACEMENT", measurement_direction_normal=[1, 0, 0])
                        )
                    elif position[3] == "Y":
                        self.sensor_list.append(
                            SensorDataContainer(position_of_mesh_node=[position[0], position[1], position[2]], type_of_sensor="DISPLACEMENT", measurement_direction_normal=[0, 1, 0])
                        )

            ############################################################################

                self.pre_operations(folder_suffix=simulation_description)

                self.run_simulations()

                self.post_operations(folder_suffix=simulation_description)

            ############################################################################


sim_setup = SimulationSetup()
sim_setup.run()
