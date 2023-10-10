import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import numpy as np
import plotly.express as px


class EffectiveSensorPositionsCalculator():

    """_summary_
    Class to wrap the calculation of 'ideal' sensor positions.
    Primary function to use is 'calculate_sensor_positions'. Returns the sensor positions as a 2d list of floats [[x,y,z], ...]
    """

    def __init__(self):
        pass

    def run_primal_analysis(self, primal_settings_file_name: str, model_file_name: str) -> None:

        # Running the primal is necessary for running an adjoint analysis

        with open(primal_settings_file_name, 'r') as parameter_file:
            simulation_parameters = Kratos.Parameters(parameter_file.read())

        new_settings = """{"input_type": %s,
            "input_filename":%s }""" % ('"mdpa"', '"'+model_file_name+'"')
        simulation_parameters["solver_settings"]["model_import_settings"] = Kratos.Parameters(new_settings)

        model = Kratos.Model()
        simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(model, simulation_parameters)
        simulation.Initialize()
        simulation.RunSolutionLoop()
        simulation.Finalize()

    def change_single_node_model_part_to_node_id(self, model_file_name: str, model_part: Kratos.ModelPart, node_id: int) -> None:

        # The "single_node" model part is a dummy model part to help with calculating the gradient with respect to one nodal displacement
        # This method helps with selecting the next node for the gradient calculation

        # Create "single_node"-model part if needed and remove all nodes from it
        if "Structure.single_node" not in model_part.GetModel().GetModelPartNames():
            single_node_model_part = model_part.CreateSubModelPart("single_node")
        else:
            single_node_model_part = model_part.GetSubModelPart("single_node")
        for node in single_node_model_part.GetNodes():
            if node is not None:
                single_node_model_part.RemoveNode(node)
        # Add new node to model part
        single_node_model_part.AddNodes([node_id])
        # Write the model to the original file
        Kratos.ModelPartIO(str(model_file_name), Kratos.ModelPartIO.WRITE | Kratos.ModelPartIO.MESH_ONLY).WriteModelPart(model_part)

    def get_node_ids_in_model_part(self, model_part: Kratos.ModelPart, model_part_name: str) -> "list[int]":
        ids = []
        for node in model_part.GetSubModelPart(model_part_name).GetNodes():
            ids.append(node.Id)
        return ids

    def normalize_sensitivities_to_element_area(self, expression: Kratos.Expression) -> None:
        sensitivities = expression.Evaluate()

        for i, element in enumerate(expression.GetModelPart().GetElements()):
            nodes = element.GetNodes()
            area = Kratos.Triangle3D3(nodes[0], nodes[1], nodes[2]).Area()
            sensitivities[i] /= area

        Kratos.Expression.CArrayExpressionIO.Read(expression, sensitivities)

    def calculate_displacement_to_young_modulus_sensitivity_matrix(self, model_part: Kratos.ModelPart, model_file_name: str, primal_settings_file_name: str, adjoint_settings_file_name: str, potential_positions_model_part: str) -> "(np.ndarray, list(int))":

        # Calculates the jacobian of d_u/d_E for all node-element combinations.
        # The "potential_positions_model_part" parameter specifies to which nodes/u's a gradient should be calculated

        self.run_primal_analysis(primal_settings_file_name=primal_settings_file_name, model_file_name=model_file_name)

        with open(adjoint_settings_file_name, 'r') as parameter_file:
            adjoint_parameters = Kratos.Parameters(parameter_file.read())

        adjoint_parameters["solver_settings"]["sensitivity_settings"]["element_data_value_sensitivity_variables"].SetStringArray(["YOUNG_MODULUS"])

        adjoint_parameters["solver_settings"]["response_function_settings"]["response_type"].SetString("adjoint_nodal_displacement")

        adjoint_parameters["solver_settings"]["response_function_settings"]["response_part_name"].SetString("single_node")

        node_ids_of_potential_sensor_positions = self.get_node_ids_in_model_part(model_part=model_part, model_part_name=potential_positions_model_part)
        sensitivities_per_node = np.zeros((len(model_part.GetElements()), len(node_ids_of_potential_sensor_positions)))

        for i, node_id in enumerate(node_ids_of_potential_sensor_positions):

            self.change_single_node_model_part_to_node_id(model_file_name=model_file_name, model_part=model_part, node_id=node_id)

            adjoint_model = Kratos.Model()
            adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(adjoint_model, adjoint_parameters)
            adjoint_model_part: Kratos.ModelPart = adjoint_analysis._GetSolver().GetComputingModelPart()

            adjoint_analysis.Run()

            # now fill the collective expressions
            adjoint_young_modulus_sensitivity = Kratos.Expression.ElementExpression(adjoint_model_part)
            Kratos.Expression.VariableExpressionIO.Read(adjoint_young_modulus_sensitivity, KratosOA.YOUNG_MODULUS_SENSITIVITY)
            self.normalize_sensitivities_to_element_area(adjoint_young_modulus_sensitivity)
            sensitivities = adjoint_young_modulus_sensitivity.Evaluate()

            sensitivities_per_node[:, i] = sensitivities

        return sensitivities_per_node, node_ids_of_potential_sensor_positions

    def get_positions_from_node_ids(self, model_part: Kratos.ModelPart, node_ids: "list[int]"):

        node_positions = []
        for ID in node_ids:
            node = model_part.GetNode(ID)
            node_positions.append([node.X, node.Y, node.Z])

        return node_positions

    def output_vector_as_vtu(self, file_name: str, vector: np.ndarray, model_part: Kratos.ModelPart) -> None:

        # Outputs a vector with size of elements in the model part as a vtu file

        expression = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(expression, vector)

        Kratos.Expression.VariableExpressionIO.Write(expression, KratosOA.YOUNG_MODULUS_SENSITIVITY)

        vtu_output_process = Kratos.VtuOutput(model_part, True, Kratos.VtuOutput.BINARY, 9)
        vtu_output_process.AddNonHistoricalVariable(KratosOA.YOUNG_MODULUS_SENSITIVITY, Kratos.VtuOutput.ELEMENTS)
        vtu_output_process.PrintOutput(file_name)

    def calculate_sensor_positions(self, model_file_name: str, primal_settings_file_name: str, adjoint_settings_file_name: str, potential_positions_model_part: str, num_of_sensors_to_select: int, redundancy_factor: float = 0.0) -> "list(list(float))":

        # Read model from file
        model = Kratos.Model()
        model_part = model.CreateModelPart("Structure")
        Kratos.ModelPartIO(model_file_name, Kratos.ModelPartIO.READ | Kratos.ModelPartIO.MESH_ONLY).ReadModelPart(model_part)

        sensitivities_per_node, node_ids_of_potential_sensor_positions = self.calculate_displacement_to_young_modulus_sensitivity_matrix(
            model_part=model_part,
            model_file_name=model_file_name,
            primal_settings_file_name=primal_settings_file_name,
            adjoint_settings_file_name=adjoint_settings_file_name,
            potential_positions_model_part=potential_positions_model_part)

        # Map sensitivities to the adjacent information space
        def mapping_to_information_space_function(x): return 1.0 if x > 1e-13 else -1.0
        vectorized_mapping = np.vectorize(mapping_to_information_space_function)
        sensitivities_per_node = vectorized_mapping(sensitivities_per_node)

        # Helper functions
        def cap_to_zero_function(x): return 0.0 if x < 0.0 else x
        vectorized_cap_to_zero_function = np.vectorize(cap_to_zero_function)
        def make_zero_to_one_function(x): return 1.0 if x == 0.0 else 0.0
        vectorized_make_zero_to_one_function = np.vectorize(make_zero_to_one_function)
        def make_zero_or_one(x): return 1.0 if x > 1e-13 else 0.0
        vectorized_make_zero_or_one = np.vectorize(make_zero_or_one)

        selected_sensor_node_ids = []
        observation_map = np.zeros((sensitivities_per_node.shape[0]))

        redundancy_weight = 0.75  # 0 only new. elements count <-> 1 equally weight of new. and red. elements count <-> inf only red. elements count
        best_picks_batch = 10  # Batch size that is selected based on the score ranking

        for n in range(num_of_sensors_to_select):
            # fig = px.imshow(sensitivities_per_node)
            # fig.show()
            new_red_unobs_map = (sensitivities_per_node.T - observation_map).T  # Is <0 for all unobserved elements, ==0 for all redundant element and >0 for all newly observed elements
            # fig = px.imshow(new_red_unobs_map)
            # fig.show()

            new_map = vectorized_cap_to_zero_function(new_red_unobs_map)
            # fig = px.imshow(new_map)
            # fig.show()
            red_map = vectorized_make_zero_to_one_function(new_red_unobs_map)
            # fig = px.imshow(red_map)
            # fig.show()

            n_new_visible_elements_score = np.sum(new_map, axis=0)+1
            n_red_visible_elements_score = np.sum(red_map, axis=0)+1
            fig = px.imshow([n_new_visible_elements_score, n_red_visible_elements_score])
            fig.show()

            total_score = np.abs((n_new_visible_elements_score/n_red_visible_elements_score)-redundancy_weight)
            fig = px.imshow([total_score])
            fig.show()

            sorted_indices = np.argsort(total_score)
            best_picks = sorted_indices[-(best_picks_batch+1):-1]

            if n == 0:
                max_or_min = 0.5
                target_score = max_or_min * total_score.max() + (1-max_or_min) * total_score.min()
                index_closest_to_target = -1
                last_closest_ds = 10e16
                for m, score in enumerate(total_score):
                    delta_score = abs(target_score - score)
                    if delta_score < last_closest_ds:
                        last_closest_ds = delta_score
                        index_closest_to_target = m

                best_picks = [index_closest_to_target]

            max_visible_elements = 0
            best_id = -1
            for pick in best_picks:
                if max_visible_elements < n_new_visible_elements_score[pick]:
                    max_visible_elements = n_new_visible_elements_score[pick]
                    best_id = pick

            print("Selected sensor: ",best_id, n_new_visible_elements_score[best_id])
            selected_sensor_node_ids.append(best_id)
            observation_map += new_map[:, best_id]

            # observation_score = np.sum(sensitivities_per_node, axis=0)
            # index_of_max = np.argmax(observation_score)
            # # put to separate vector for output
            # observation_map += vectorized_make_zero_or_one(sensitivities_per_node[:, index_of_max])
            # selected_sensor_node_ids.append(node_ids_of_potential_sensor_positions[index_of_max])
            # # remove information gathered by the sensor
            # sensitivities_per_node = (sensitivities_per_node.T - sensitivities_per_node[:, index_of_max] * (1-redundancy_factor)).T
            # sensitivities_per_node = vectorized_cap_to_zero_function(sensitivities_per_node)
            # print(f"EffectiveSensorPositionsCalculator:: Selected sensor adds an information score of {observation_score[index_of_max]}")

        self.output_vector_as_vtu(f"observability_map_n{num_of_sensors_to_select}", vector=observation_map, model_part=model_part)

        selected_sensor_positions = self.get_positions_from_node_ids(model_part=model_part, node_ids=selected_sensor_node_ids)

        print(f"EffectiveSensorPositionsCalculator:: Selected sensor positions: {selected_sensor_positions}")

        return selected_sensor_positions
