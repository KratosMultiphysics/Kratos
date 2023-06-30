import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

from Measurement_file_creation.MeasurementDataGenerator import MeasurementDataGenerator
from Measurement_file_creation.MeasurementFileStructure import *
import json


class MeasurementFromSimulationDataGenerator(MeasurementDataGenerator):

    def __init__(self, sensor_and_load_information: MeasurementFile):
        self.sensor_and_load_information = sensor_and_load_information

    def write_measurement_data_file(self, simulation_parameters: Kratos.Parameters, file_name: str):

        for load_info, sensors_info in self.sensor_and_load_information.load_cases:

            if load_info.type_of_load == "PointLoad":
                load_parameters = self.get_point_load_parameters(load_info)
            else:
                RuntimeError("The automatic measurement file creation from simulation data currently only supports 'PointLoad's")

            model = Kratos.Model()
            simulation_parameters["processes"]["loads_process_list"] = Kratos.Parameters(load_parameters)

            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Run()

        #     for sensor in self.measurement_data["load_cases"][0]["sensors_infos"]:

        #     point = Kratos.Point(sensor["position_of_mesh_node"])
        #     search_configuration = Kratos.Configuration.Initial
        #     search_tolerance = 1e-6

        #     point_locator = Kratos.BruteForcePointLocator(self.model_part)

        #     found_node_id = point_locator.FindNode(
        #         point, search_configuration, search_tolerance
        #     )

        #     sensor["mesh_node_id"] = found_node_id
        #     # Add node to sensor_positions_model_part
        #     self.sensor_positions_model_part.AddNode(self.model_part.GetNode(found_node_id))

        # # Update measurement data file with newly added infos
        # with open(self.measurement_data_file, "w") as outfile:
        #     json.dump(self.measurement_data, outfile)

    def get_point_load_parameters(self, load_info: load_info) -> Kratos.Parameters:

        load_json = """
                [{
                  "python_module": "assign_vector_by_direction_to_condition_process",
                  "kratos_module": "KratosMultiphysics",
                  "help": "This process sets a vector variable value over a condition",
                  "check": "DirectorVectorNonZero direction",
                  "process_name": "AssignModulusAndDirectionToConditionsProcess",
                  "Parameters": {
                    "mesh_id": 0,
                    "model_part_name": "Structure.point_load_positions_model_part",
                    "variable_name": "POINT_LOAD",
                    "modulus": 1.0,
                    "direction": [0.0, 0.0, 0.0],
                    "interval": [0.0, "End"]
                  }
                }]
                    """

        load_parameters = Kratos.Parameters(load_json)
        load_parameters[0]["Parameters"]["modulus"] = Kratos.Parameters(str(load_info.strength_in_N))
        load_parameters[0]["Parameters"]["direction"] = Kratos.Parameters(str(load_info.direction_normal))

        return load_parameters
