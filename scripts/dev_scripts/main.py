import os
import shutil
from glob import glob
import math

import numpy as np


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
            # "benchmark1/case1a/1_medium",
            # "benchmark1/case1a/2_fine", #
            # "benchmark1/case1b/0_coar",
            # "benchmark1/case1b/1_medium",
            # "benchmark1/case1b/2_fine", #
            "benchmark1/case1c/0_coar",
            # "benchmark1/case1c/1_medium", #
            # "benchmark1/case1c/2_fine", #
            # "benchmark1/case1d/0_coar",
            # "benchmark1/case1d/1_medium", #
            # "benchmark1/case1d/2_fine",
            # "benchmark2/case2a/0_coar",
            # "benchmark2/case2b/0_coar",
            # "benchmark2/case2c/0_coar",
            # "benchmark3/case3a/0_coar",
            # "benchmark3/case3b/0_coar",
            # "benchmark3/case3c/0_coar",
        ]

        xy_grid = []

        # # # For BM 1
        # for x in np.linspace(0.5, 4.5, 9):
        #     for y in np.linspace(0.0, 0.5, 3):
        #         xy_grid.append([x, y, 0])

        # # For BM 2
        # for x in np.linspace(10,50,5):  # (10,50,5): #(5, 55, 10):
        #     for y in np.linspace(5,25,3):  # (5,25,3): #(5, 25, 5):
        #         # Exclude circle
        #         off_set_x = x - 30
        #         off_set_y = y - 15
        #         r = math.sqrt(off_set_x**2 + off_set_y**2)
        #         if r <= 5.0:
        #             continue

        #         xy_grid.append([x, y, 0])

        # For BM 3
        # for x in np.linspace(0.5,13.5,5):
        #     for y in np.linspace(0.5,13.5,5):
        #         xy_grid.append([x, y, 0])

        self.sensor_position_configurations_list = [

            # xy_grid,

            # [[0.5, 0.5, 0.0],
            #  [0.75, 0.5, 0.0],
            #  [1.0, 0.5, 0.0],
            #  [1.25, 0.5, 0.0],
            #  [1.5, 0.5, 0.0],
            #  [1.75, 0.5, 0.0],
            #  [2.0, 0.5, 0.0],
            #  [2.25, 0.5, 0.0],
            #  [2.5, 0.5, 0.0],
            #  [2.75, 0.5, 0.0],
            #  [3.0, 0.5, 0.0],
            #  [3.25, 0.5, 0.0],
            #  [3.5, 0.5, 0.0],
            #  [3.75, 0.5, 0.0],
            #  [4.0, 0.5, 0.0],
            #  [4.25, 0.5, 0.0],
            #  [4.5, 0.5, 0.0],
            #  [0.5, 0.0, 0.0],
            #  [0.75, 0.0, 0.0],
            #  [1.0, 0.0, 0.0],
            #  [1.25, 0.0, 0.0],
            #  [1.5, 0.0, 0.0],
            #  [1.75, 0.0, 0.0],
            #  [2.0, 0.0, 0.0],
            #  [2.25, 0.0, 0.0],
            #  [2.5, 0.0, 0.0],
            #  [2.75, 0.0, 0.0],
            #  [3.0, 0.0, 0.0],
            #  [3.25, 0.0, 0.0],
            #  [3.5, 0.0, 0.0],
            #  [3.75, 0.0, 0.0],
            #  [4.0, 0.0, 0.0],
            #  [4.25, 0.0, 0.0],
            #  [4.5, 0.0, 0.0]],

            # [[0.5, 0.5, 0.0],
            #  [0.75, 0.5, 0.0],
            #  [1.0, 0.5, 0.0],
            #  [1.25, 0.5, 0.0],
            #  [1.5, 0.5, 0.0],
            #  [1.75, 0.5, 0.0],
            #  [2.0, 0.5, 0.0],
            #  [2.25, 0.5, 0.0],
            #  [2.5, 0.5, 0.0],
            #  [2.75, 0.5, 0.0],
            #  [3.0, 0.5, 0.0],
            #  [3.25, 0.5, 0.0],
            #  [3.5, 0.5, 0.0],
            #  [3.75, 0.5, 0.0],
            #  [4.0, 0.5, 0.0],
            #  [4.25, 0.5, 0.0],
            #  [4.5, 0.5, 0.0]],

            [[0.5, 0.5, 0.0],
             [1.0, 0.5, 0.0],
             [1.5, 0.5, 0.0],
             [2.0, 0.5, 0.0],
             [2.5, 0.5, 0.0],
             [3.0, 0.5, 0.0],
             [3.5, 0.5, 0.0],
             [4.0, 0.5, 0.0],
             [4.5, 0.5, 0.0]],

            # [[0.5, 0.5, 0.0],
            #  [1.5, 0.5, 0.0],
            #  [2.5, 0.5, 0.0],
            #  [3.5, 0.5, 0.0],
            #  [4.5, 0.5, 0.0]],

            ##################
            # xy_grid,

            # [[50, 25, 0],
            #  [50, 5, 0],
            #  [40, 25, 0],
            #  [40, 5, 0],
            #  [30, 25, 0],
            #  [30, 5, 0],
            #  [20, 25, 0],
            #  [20, 5, 0],
            #  [10, 25, 0],
            #  [10, 5, 0],
            #  [10, 15, 0],
            #  [20, 15, 0],
            #  [40, 15, 0],
            #  [50, 15, 0],],

            ##################

            # xy_grid,

            # [[12, 3, 0]
            #  ]

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

            simulation_description = "smoothing3_msY"

            os.chdir(oldPath)

            self.cascading_copy_simulation_files(working_folder)

            os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

            # sensor_placer = EffectiveSensorPositionsCalculator()
            # auto_found_sensor_positions = sensor_placer.calculate_sensor_positions(model_file_name="model_file",
            #                                                                         primal_settings_file_name="primal_parameters.json",
            #                                                                         adjoint_settings_file_name="adjoint_parameters.json",
            #                                                                         potential_positions_model_part="all_nodes_elements_model_part",
            #                                                                         num_of_sensors_to_select=5,
            #                                                                         redundancy_factor=0.2)
            # self.sensor_position_configurations_list.append(auto_found_sensor_positions)

            for sensor_config in self.sensor_position_configurations_list:

                ############################################################################

                self.sensor_list.clear()
                for position in sensor_config:
                    # self.sensor_list.append(
                    #     SensorDataContainer(position_of_mesh_node=position, type_of_sensor="DISPLACEMENT", measurement_direction_normal=[1, 0, 0])
                    # )
                    self.sensor_list.append(
                        SensorDataContainer(position_of_mesh_node=position, type_of_sensor="DISPLACEMENT", measurement_direction_normal=[0, 1, 0])
                    )

            ############################################################################

                self.pre_operations(folder_suffix=simulation_description)

                self.run_simulations()

                self.post_operations(folder_suffix=simulation_description)

            ############################################################################


sim_setup = SimulationSetup()
sim_setup.run()
