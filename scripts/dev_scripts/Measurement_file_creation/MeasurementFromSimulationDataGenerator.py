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

            model = Kratos.Model()
            # Adjust import settings for the model
            new_settings = """{"input_type": %s,
            "input_filename":%s }"""%('"mdpa"','"'+mdpa_model_file_name+'"')
            simulation_parameters["solver_settings"]["model_import_settings"] = Kratos.Parameters(new_settings)

            # Create dummy analysis to be able to access model part and to find node
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Initialize()
            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            model = Kratos.Model()
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Run()
            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            # Transfer data to sensor data in dataclass
            for sensor in sensors_infos:
                self.update_sensor_infos(sensor=sensor,model_part=model_part)

            # Update measurement data file with newly added infos
            if not output_file_name.endswith(".json"):
                output_file_name += ".json"

            with open(output_file_name, "w") as outfile:
                json.dump(dc.asdict(self.sensor_and_load_information), outfile)

    def update_sensor_infos(self,sensor_data:SensorDataContainer,model_part:Kratos.ModelPart):
            found_sensor_node_id = self.get_node_id_at_position(sensor_data.position_of_mesh_node,model_part)

            node:Kratos.Node = model_part.GetNode(found_sensor_node_id)
            sensor_data.mesh_node_id = found_sensor_node_id
            try:
                var = Kratos.KratosGlobals.GetVariable(sensor_data.type_of_sensor)
            except:
                raise RuntimeError(f"Please provide a sensor type that matches a Kratos global variable (like: DISPLACEMENT). {sensor_data.type_of_sensor} was not found")
            simulated_displacement = node.GetSolutionStepValue(var)
            measurement_normal = sensor_data.measurement_direction_normal
            sensor_data.measured_value = simulated_displacement[0]*measurement_normal[0]+simulated_displacement[1]*measurement_normal[1]+simulated_displacement[2]*measurement_normal[2]

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
