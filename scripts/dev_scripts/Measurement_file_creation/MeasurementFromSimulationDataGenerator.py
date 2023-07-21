import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

from Measurement_file_creation.MeasurementDataGenerator import MeasurementDataGenerator
from Measurement_file_creation.MaterialChanger import MaterialChanger
from Measurement_file_creation.MeasurementFileStructure import *
import json
import dataclasses as dc
import math


class MeasurementFromSimulationDataGenerator(MeasurementDataGenerator):

    def __init__(self, sensor_and_load_information: LoadCasesContainer, material_changer: MaterialChanger = None):
        self.sensor_and_load_information = sensor_and_load_information
        self.material_changer = material_changer

    def write_measurement_data_file(self, simulation_parameters: Kratos.Parameters, mdpa_model_file_name: str, output_file_name: str):

        print("MeasurementFromSimulationDataGenerator:: Start writing process")

        for i, load_case in enumerate(self.sensor_and_load_information.load_cases):
            sensors_infos: List[SensorDataContainer] = load_case.sensors_infos

            # model = Kratos.Model()
            # Adjust import settings for the model
            new_settings = """{"input_type": %s,
            "input_filename":%s }""" % ('"mdpa"', '"'+mdpa_model_file_name+'"')
            simulation_parameters["solver_settings"]["model_import_settings"] = Kratos.Parameters(new_settings)

            # Create dummy analysis to be able to access model part and to find node
            # simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            # simulation.Initialize()
            # model_part = model.GetModelPart(model.GetModelPartNames()[0])

            model = Kratos.Model()
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Initialize()
            print("MeasurementFromSimulationDataGenerator:: Finished StructuralMechanicsAnalysis initialization")
            if not (self.material_changer == None):
                self.material_changer.adjust_material_of_model_part(model.GetModelPart("Structure"))
            print("MeasurementFromSimulationDataGenerator:: Finished adjusting the materials")
            simulation.RunSolutionLoop()
            print("MeasurementFromSimulationDataGenerator:: Finished StructuralMechanicsAnalysis solution loop")
            simulation.Finalize()
            print("MeasurementFromSimulationDataGenerator:: Finished StructuralMechanicsAnalysis finalize")
            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            print("MeasurementFromSimulationDataGenerator:: Start extraction of sensor data")
            # Transfer data to sensor data in dataclass
            model_part.CreateSubModelPart("sensors")
            for sensor in sensors_infos:
                self.update_sensor_infos_and_add_to_sensor_model_part(sensor_data=sensor, model_part=model_part, sensor_model_part=model_part.GetSubModelPart("sensors"))

            # Update measurement data file with newly added infos
            if not output_file_name.endswith(".json"):
                output_file_name += ".json"

            print("MeasurementFromSimulationDataGenerator:: Start writing measurement json")
            with open(output_file_name, "w") as outfile:
                json.dump(dc.asdict(self.sensor_and_load_information), outfile)
            print("MeasurementFromSimulationDataGenerator:: Finished writing measurement json")

            self.output_measurement_model(model, f"lc{i}")

    def update_sensor_infos_and_add_to_sensor_model_part(self, sensor_data: SensorDataContainer, model_part: Kratos.ModelPart, sensor_model_part: Kratos.ModelPart):
        found_sensor_node_id = self.get_node_id_closest_to_position(sensor_data.position_of_mesh_node, model_part)

        if found_sensor_node_id == -1:
            raise ValueError(f"Sensor location {sensor_data.position_of_mesh_node} does not match a node location of the provided mesh.")

        node: Kratos.Node = model_part.GetNode(found_sensor_node_id)
        sensor_data.mesh_node_id = found_sensor_node_id
        sensor_data.position_of_mesh_node = [node.X - node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0],
                                             node.Y - node.GetSolutionStepValue(Kratos.DISPLACEMENT)[1],
                                             node.Z - node.GetSolutionStepValue(Kratos.DISPLACEMENT)[2]]
        sensor_model_part.AddNode(node)

        try:
            var = Kratos.KratosGlobals.GetVariable(sensor_data.type_of_sensor)
        except:
            raise RuntimeError(f"Please provide a sensor type that matches a Kratos global variable (like: DISPLACEMENT). {sensor_data.type_of_sensor} was not found")

        simulated_displacement = node.GetSolutionStepValue(var)
        measurement_normal = sensor_data.measurement_direction_normal
        sensor_data.measured_value = simulated_displacement[0]*measurement_normal[0]+simulated_displacement[1]*measurement_normal[1]+simulated_displacement[2]*measurement_normal[2]

    def get_node_id_closest_to_position(self, position: List[float], model_part: Kratos.ModelPart) -> int:
        point = Kratos.Point(position)

        # Search for points in mesh
        search_configuration = Kratos.Configuration.Initial
        search_tolerance = 1e-6
        point_locator = Kratos.BruteForcePointLocator(model_part)
        found_element_id = point_locator.FindElement(
            point, Kratos.Vector(), search_configuration, search_tolerance
        )

        # Then the point lies outside the mesh or on a node -> so we try to receive the node id directly
        if found_element_id == -1:
            found_node_id = point_locator.FindNode(
                point, search_configuration, search_tolerance
            )
            return found_node_id

        closest_node_id = None
        last_closest_distance = 10e20
        current_distance = 0
        dx = 0
        dy = 0
        dz = 0
        for node in model_part.GetElement(found_element_id).GetNodes():
            dx = position[0] - node.X-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0]
            dy = position[1] - node.Y-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[1]
            dz = position[2] - node.Z-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[2]

            current_distance = math.sqrt(dx*dx + dy*dy + dz*dz)

            if last_closest_distance > current_distance:
                last_closest_distance = current_distance
                closest_node_id = node.Id

        return closest_node_id

    def output_measurement_model(self, model: Kratos.Model, file_name_prefix: str) -> None:
        for model_part_name in model.GetModelPartNames():
            model_part = model.GetModelPart(model_part_name)
            vtu_output_process = Kratos.VtuOutput(model_part, True, Kratos.VtuOutput.BINARY, 9)

            # Assign youngs modulus to element values
            for element in model_part.Elements:
                element.SetValue(Kratos.YOUNG_MODULUS, element.Properties[Kratos.YOUNG_MODULUS])
            if len(model_part.Elements) > 0:
                vtu_output_process.AddNonHistoricalVariable(Kratos.YOUNG_MODULUS, Kratos.VtuOutput.ELEMENTS)

            # Assign displacement to node values
            for node in model_part.Nodes:
                node.SetValue(Kratos.DISPLACEMENT, node.GetSolutionStepValue(Kratos.DISPLACEMENT))
            if len(model_part.Nodes) > 0:
                vtu_output_process.AddNonHistoricalVariable(Kratos.DISPLACEMENT, Kratos.VtuOutput.NODES)

            vtu_output_process.PrintOutput(f"{file_name_prefix}_{model_part_name}")
