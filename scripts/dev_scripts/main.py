import os
import shutil
from glob import glob

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
        material_changer = PerModelPartMaterialChanger(self.material_parameters_to_change)

        data_generator = MeasurementFromSimulationDataGenerator()

        self.measurement_load_cases_container = data_generator.write_measurement_data_file(
            mdpa_model_file_name="model_file",
            primal_parameter_files=["primal_parameters.json"],
            sensor_positions=self.sensor_list,
            output_file_name="measurement_data.json",
            material_changer=material_changer)

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
            for file in glob(f"{cumulative_folder_path}/*_parameters.json")+glob(f"{cumulative_folder_path}/StructuralMaterials.json"):
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
        for file in glob('./Structure*.h5')+glob('./all_nodes_elements_model_part*.h5'):
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
        #     print(f"Error in analysis: {working_folder}. Continue running other examples")

    def run(self) -> None:

        self.simulation_folders = [
            "benchmark1/case1a/0_coar"
        ]

        self.sensor_position_configurations_list = [
            [[0.5, 0.5, 0.0],
             [0.75, 0.5, 0.0],
             [1.0, 0.5, 0.0],
             [1.25, 0.5, 0.0],
             [1.5, 0.5, 0.0],
             [1.75, 0.5, 0.0],
             [2.0, 0.5, 0.0],
             [2.25, 0.5, 0.0],
             [2.5, 0.5, 0.0],
             [2.75, 0.5, 0.0],
             [3.0, 0.5, 0.0],
             [3.25, 0.5, 0.0],
             [3.5, 0.5, 0.0],
             [3.75, 0.5, 0.0],
             [4.0, 0.5, 0.0],
             [4.25, 0.5, 0.0],
             [4.5, 0.5, 0.0]],
        ]

        self.material_parameters_to_change = {
            "ElemMat1": {"YOUNG_MODULUS": 30000000000.0},  # BM1
            "ElemMat2": {"YOUNG_MODULUS": 10000000000.0},  # BM1
        }

        ############################################################################

        oldPath = os.getcwd()

        for working_folder in self.simulation_folders:

            simulation_description = "no_smoothing"

            os.chdir(oldPath)

            self.cascading_copy_simulation_files(working_folder)

            os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

            # sensor_placer = EffectiveSensorPositionsCalculator()
            # auto_found_sensor_positions = sensor_placer.calculate_sensor_positions(model_file_name="model_file",
            #                                                                         primal_settings_file_name="primal_parameters.json",
            #                                                                         adjoint_settings_file_name="adjoint_parameters.json",
            #                                                                         potential_positions_model_part="all_nodes_elements_model_part",
            #                                                                         num_of_sensors_to_select=3,
            #                                                                         redundancy_factor=0.2)
            # self.sensor_position_configurations_list.append(auto_found_sensor_positions)

            for sensor_config in self.sensor_position_configurations_list:

                ############################################################################

                for position in sensor_config:
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
