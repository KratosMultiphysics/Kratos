import os

from Measurement_file_creation.MeasurementFromSimulationDataGenerator import MeasurementFromSimulationDataGenerator

from Measurement_file_creation.MeasurementFileStructure import *

import KratosMultiphysics as Kratos

s1 = sensor_infos(position_of_mesh_node=[1.0, 0.5, 0.0])
s2 = sensor_infos(position_of_mesh_node=[2.0, 0.5, 0.0])
l = load_info(strength_in_N=-10000000.0, position_of_mesh_vertex=[2.5, 0.5, 0.0])
s_list = sensors_infos([s1, s2])

raw_measurement_file = MeasurementFile([(l, s_list)])

data_generator = MeasurementFromSimulationDataGenerator(sensor_and_load_information=raw_measurement_file)

working_folder = "measurement_residual_test"
os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

with open("primal_parameters_import_mdpa.json", 'r') as parameter_file:
    parameters = Kratos.Parameters(parameter_file.read())

data_generator.write_measurement_data_file(simulation_parameters=parameters, file_name="test_measurement_file.json")
