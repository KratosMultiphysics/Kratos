import os
import shutil
from glob import glob

import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from Measurement_file_creation.MeasurementFileStructure import *
from Measurement_file_creation.MeasurementFromSimulationDataGenerator import MeasurementFromSimulationDataGenerator
from Measurement_file_creation.PerModelPartMaterialChanger import PerModelPartMaterialChanger

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis


def measurement_creation() -> None:
    material_changer = PerModelPartMaterialChanger(material_parameters_to_change)

    data_generator = MeasurementFromSimulationDataGenerator(sensor_and_load_information=raw_measurement_file, material_changer=material_changer)

    with open("primal_parameters.json", 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    data_generator.write_measurement_data_file(simulation_parameters=parameters, mdpa_model_file_name="model_file", output_file_name="measurement_data.json")


def get_additional_folder_name_from_measurement_data() -> str:
    n_load_cases = len(raw_measurement_file.load_cases)
    n_sensors = len(raw_measurement_file.load_cases[0].sensors_infos)
    data_name = f"_lc{n_load_cases}_s{n_sensors}"
    return data_name


def cascading_copy_simulation_files(working_folder: str) -> None:
    # os.chdir(oldPath)
    # os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

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


def pre_operations(folder_suffix: str = "") -> None:
    measurement_creation()
    data_name = get_additional_folder_name_from_measurement_data()
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


def post_operations(folder_suffix: str = "") -> None:
    data_name = get_additional_folder_name_from_measurement_data()
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


simulation_folders = [
    # "benchmark1/case1a/0_coar",
    # "benchmark1/case1a/1_medium",
    "benchmark1/case1a/2_fine",
    # "benchmark1/case1b/0_coar",
    # "benchmark1/case1b/1_medium",
    # "benchmark1/case1b/2_fine", #
    # "benchmark1/case1c/0_coar",
    # "benchmark1/case1c/1_medium", #
    # "benchmark1/case1c/2_fine", #
    # "benchmark1/case1d/0_coar",
    # "benchmark1/case1d/1_medium", #
    # "benchmark1/case1d/2_fine",
    # "benchmark2/case2a/0_coar",
    # "benchmark2/case2b/0_coar",
    # "benchmark2/case2c/0_coar",
    # "benchmark3/case3a/0_coar",
]

xy_grid = []
for x in np.linspace(0.5, 4.5, 9):
    for y in np.linspace(0.0, 0.5, 3):
        xy_grid.append([x, y, 0])

print(xy_grid)


sensor_position_configurations_list = [

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

    # [[0.5, 0.5, 0.0],
    #  [1.0, 0.5, 0.0],
    #  [1.5, 0.5, 0.0],
    #  [2.0, 0.5, 0.0],
    #  [2.5, 0.5, 0.0],
    #  [3.0, 0.5, 0.0],
    #  [3.5, 0.5, 0.0],
    #  [4.0, 0.5, 0.0],
    #  [4.5, 0.5, 0.0]],

    [[0.5, 0.5, 0.0],
     [1.5, 0.5, 0.0],
     [2.5, 0.5, 0.0],
     [3.5, 0.5, 0.0],
     [4.5, 0.5, 0.0]],

    ##################

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

    # [[12, 3, 0]]

]

oldPath = os.getcwd()

for working_folder in simulation_folders:

    simulation_description = "no_smoothing"

    os.chdir(oldPath)

    cascading_copy_simulation_files(working_folder)

    os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

    for sensor_config in sensor_position_configurations_list:

        ############################################################################

        sensor_list = []

        for position in sensor_config:
            sensor_list.append(
                SensorDataContainer(position_of_mesh_node=position, type_of_sensor="DISPLACEMENT", measurement_direction_normal=[0, 1, 0])
            )

        l = LoadDataContainer()

        raw_measurement_file = LoadCasesContainer(load_cases=[
            PerLoadCaseMeasurementDataContainer(load_info=l, sensors_infos=sensor_list)
        ])

        material_parameters_to_change = {
            "ElemMat1": {"YOUNG_MODULUS": 30000000000.0},  # BM1
            "ElemMat2": {"YOUNG_MODULUS": 10000000000.0},  # BM1
            # "ElemMat3": {"YOUNG_MODULUS": 10000000000.0},  # BM3 BM2
            # "ElemMat4": {"YOUNG_MODULUS": 10000000000.0},  # BM3 BM2
            # "ElemMat5": {"YOUNG_MODULUS": 10000000000.0},  # BM3 BM2
            # "ElemMat6": {"YOUNG_MODULUS": 10000000000.0},  # BM3 BM2
            # "ElemMat7": {"YOUNG_MODULUS": 10000000000.0},  # BM3
            # "ElemMat8": {"YOUNG_MODULUS": 10000000000.0},  # BM3
            # "ElemMat9": {"YOUNG_MODULUS": 10000000000.0},  # BM3
        }

############################################################################

        pre_operations(folder_suffix=simulation_description)

        with open("optimization_parameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())

        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)

        try:
            analysis.Run()
        except:
            print(f"Error in analysis: {working_folder}")

        post_operations(folder_suffix=simulation_description)

############################################################################
