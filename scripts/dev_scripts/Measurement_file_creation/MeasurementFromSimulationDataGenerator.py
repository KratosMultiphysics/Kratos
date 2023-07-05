import KratosMultiphysics as Kratos
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

from Measurement_file_creation.MeasurementDataGenerator import MeasurementDataGenerator
from Measurement_file_creation.MeasurementFileStructure import *
import json
import dataclasses as dc


class MeasurementFromSimulationDataGenerator(MeasurementDataGenerator):

    def __init__(self, sensor_and_load_information: LoadCasesContainer):
        self.sensor_and_load_information = sensor_and_load_information

    def write_measurement_data_file(self, simulation_parameters: Kratos.Parameters, mdpa_model_file_name:str ,output_file_name: str):

        for load_case in self.sensor_and_load_information.load_cases:
            sensors_infos:List[SensorDataContainer] = load_case.sensors_infos
            load_info:LoadDataContainer = load_case.load_info

            model = Kratos.Model()
            # Adjust import settings for the model
            new_settings = """{"input_type": %s,
            "input_filename":%s }"""%('"mdpa"','"'+mdpa_model_file_name+'"')
            simulation_parameters["solver_settings"]["model_import_settings"] = Kratos.Parameters(new_settings)

            # Create dummy analysis to be able to access model part and to find node
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Initialize()
            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            # Adjust load data in parameters and dataclass
            found_load_node_id = self.get_node_id_at_position(load_info.position_of_mesh_vertex,model_part)
            load_info.mesh_node_id = found_load_node_id
            if load_info.type_of_load == "PointLoad":
                load_parameters = self.get_point_load_parameters(load_info)
            else:
                raise RuntimeError("The automatic measurement file creation from simulation data currently only supports 'PointLoad's")

            model = Kratos.Model()
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Initialize()
            simulation.RunSolutionLoop()
            simulation.Finalize()
            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            # Transfer data to sensor data in dataclass
            for sensor in sensors_infos:

                found_sensor_node_id = self.get_node_id_at_position(sensor.position_of_mesh_node,model_part)

                node:Kratos.Node = model_part.GetNode(found_sensor_node_id)
                sensor.mesh_node_id = found_sensor_node_id
                try:
                    var = Kratos.KratosGlobals.GetVariable(sensor.type_of_sensor)
                except:
                    raise RuntimeError(f"Please provide a sensor type that matches a Kratos global variable (like: DISPLACEMENT). {sensor.type_of_sensor} was not found")
                simulated_displacement = node.GetSolutionStepValue(var)
                measurement_normal = sensor.measurement_direction_normal
                sensor.measured_value = simulated_displacement[0]*measurement_normal[0]+simulated_displacement[1]*measurement_normal[1]+simulated_displacement[2]*measurement_normal[2]


            # Update measurement data file with newly added infos
            if not output_file_name.endswith(".json"):
                output_file_name += ".json"

            with open(output_file_name, "w") as outfile:
                json.dump(dc.asdict(self.sensor_and_load_information), outfile)

    def get_point_load_parameters(self, load_info: LoadDataContainer) -> Kratos.Parameters:

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

    def get_node_id_at_position(self,position:List[float],model_part:Kratos.ModelPart) -> int:
        point = Kratos.Point(position)

        # Search for points in mesh
        search_configuration = Kratos.Configuration.Initial
        search_tolerance = 1e-6
        point_locator = Kratos.BruteForcePointLocator(model_part)
        found_node_id = point_locator.FindNode(
            point, search_configuration, search_tolerance
        )

        return found_node_id
