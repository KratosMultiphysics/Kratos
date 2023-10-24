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

    """_summary_
    Implementation of the MeasurementDataGenerator interface. Class is capable of running Kratos simulations based on provided parameter files and extracting simulations results at the positions of sensors.
    The main function to use is 'write_measurement_data_file()'. It will write a measurement.json file based on the simulation parameters and the sensor_positions provided.
    """

    def __init__(self):
        pass

    def write_measurement_data_file(self, mdpa_model_file_name: str, primal_parameter_files: "list[str]", sensor_positions: "list[SensorDataContainer]", output_file_name: str, material_changer: MaterialChanger = None) -> LoadCasesContainer:
        self.material_changer = material_changer
        self.load_cases_list = []

        print("MeasurementFromSimulationDataGenerator:: Start writing process")

        # Loop over all provided primal parameter files
        for i, parameter_file in enumerate(primal_parameter_files):

            with open(parameter_file, 'r') as parameter_file:
                simulation_parameters = Kratos.Parameters(parameter_file.read())

            # Adjust import settings for the model so that it gets loaded fresh every time.
            new_settings = """{"input_type": %s,
            "input_filename":%s }""" % ('"mdpa"', '"'+mdpa_model_file_name+'"')
            simulation_parameters["solver_settings"]["model_import_settings"] = Kratos.Parameters(new_settings)

            # Model and simulation initialization
            model = Kratos.Model()
            simulation = StructuralMechanicsAnalysis(model, simulation_parameters)
            simulation.Initialize()

            # If specified - let the material changer do its thing
            if not (self.material_changer == None):
                self.material_changer.adjust_material_of_model_part(model.GetModelPart("Structure"))

            # Rest of the simulation procedure
            simulation.RunSolutionLoop()
            simulation.Finalize()

            model_part = model.GetModelPart(model.GetModelPartNames()[0])

            # Transfer data to sensor data in dataclass
            model_part.CreateSubModelPart("sensors")
            current_sensor_infos = []
            for sensor in sensor_positions:
                self.update_sensor_infos_and_add_to_sensor_model_part(sensor_data=sensor, model_part=model_part, sensor_model_part=model_part.GetSubModelPart("sensors"))
                current_sensor_infos.append(sensor.copy())

            # TODO adjust
            # load_infos = self.get_load_info_from_parameters(simulation_parameters=simulation_parameters)
            load_infos = LoadDataContainer(type_of_load=parameter_file.name[:-5])

            load_case_info = PerLoadCaseMeasurementDataContainer(load_info=load_infos, sensors_infos=current_sensor_infos)
            self.load_cases_list.append(load_case_info)

            self.output_measurement_model(model, f"lc{i}")

        self.sensor_and_load_information = LoadCasesContainer(load_cases=self.load_cases_list)

        # Update measurement data file with newly added infos
        with open(output_file_name, "w") as outfile:
            json.dump(dc.asdict(self.sensor_and_load_information), outfile)
        print("MeasurementFromSimulationDataGenerator:: Finished writing measurement json")

        return self.sensor_and_load_information

    def get_load_info_from_parameters(self, simulation_parameters: Kratos.Parameters) -> LoadDataContainer:
        raise ValueError("Not implemented")
        return LoadDataContainer()

    def update_sensor_infos_and_add_to_sensor_model_part(self, sensor_data: SensorDataContainer, model_part: Kratos.ModelPart, sensor_model_part: Kratos.ModelPart):

        """_summary_
        Fills a provided SensorDataContainer with the remaining information and adds the node a sensor is placed on to a provided model part.
        """

        # Find the node that is closest to the sensors position
        found_sensor_node_id = self.get_node_id_closest_to_position(sensor_data.position_of_mesh_node, model_part)

        if found_sensor_node_id == -1:
            raise ValueError(f"Sensor location {sensor_data.position_of_mesh_node} does not match a node location of the provided mesh.")

        # Start - Update information
        node: Kratos.Node = model_part.GetNode(found_sensor_node_id)
        sensor_data.mesh_node_id = found_sensor_node_id
        # Correct the nodes position by the displacement it was moved in the simulation
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
        # Get the measured value in the direction of the specified measurement normal
        sensor_data.measured_value = simulated_displacement[0]*measurement_normal[0]+simulated_displacement[1]*measurement_normal[1]+simulated_displacement[2]*measurement_normal[2]
        # End - Update information

    def get_node_id_closest_to_position(self, position: List[float], model_part: Kratos.ModelPart) -> int:

        """_summary_
        Finds and returns the id of a node in the provided model part to which a provided position is closest to.

        Returns:
            int: Id of a node in the provided model part
        """

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

            if found_node_id == -1:
                raise RuntimeError(f"Provided position: {position} does not lie inside the provided model part ({model_part.Name}) or on a node.")

            return found_node_id

        # Sensor lies inside an element -> find the closest node
        closest_node_id = None
        last_closest_distance = 10e20
        current_distance = 0
        dx = 0
        dy = 0
        dz = 0
        # Loop over all nodes
        for node in model_part.GetElement(found_element_id).GetNodes():
            # Correct node positions by the displacement it was moved in the simulation
            dx = position[0] - node.X-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0]
            dy = position[1] - node.Y-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[1]
            dz = position[2] - node.Z-node.GetSolutionStepValue(Kratos.DISPLACEMENT)[2]

            current_distance = math.sqrt(dx*dx + dy*dy + dz*dz)

            if last_closest_distance > current_distance:
                last_closest_distance = current_distance
                closest_node_id = node.Id

        return closest_node_id

    def output_measurement_model(self, model: Kratos.Model, file_name_prefix: str) -> None:

        """_summary_
        Outputs the model which was used in the generation of the artificial measurement values in vtu format.
        Currently only the Kratos.YOUNG_MODULUS and Kratos.DISPLACEMENT values are saved
        """

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
