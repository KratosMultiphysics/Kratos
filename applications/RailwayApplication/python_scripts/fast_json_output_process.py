from KratosMultiphysics.json_output_process import JsonOutputProcess
import KratosMultiphysics

import orjson


def fast_read_external_json(file_name: str):
    """
    Reads a JSON file using orjson for faster performance.

    Args:
        - file_name (str): The name of the input JSON file.
    """
    with open(file_name, 'rb') as outfile:
        data = orjson.loads(outfile.read())

    return data

def fast_write_external_json(file_name: str, data: dict):
    """
    Writes data to a JSON file using orjson for faster performance.

    Args:
        - file_name (str): The name of the output JSON file.
        - data (dict): The data to write to the JSON file.
    """
    with open(file_name, 'wb') as outfile:
        outfile.write(orjson.dumps(data, option=orjson.OPT_INDENT_2))

class FastJsonOutputProcess(JsonOutputProcess):
    """
    This class is used in order to create a json file containing
    the solution a given model part with a certain frequency

    This class is a fast version of the JsonOutputProcess, which uses orjson
    for reading and writing the json file. This is significantly faster than the
    standard JsonOutputProcess, especially for large models or frequent output.

    - Inheritance:
        :class:`KratosMultiphysics.JsonOutputProcess`

    """

    def __init__(self, model: KratosMultiphysics.Model, params: KratosMultiphysics.Parameters):
        """
        Constructor of the FastJsonOutputProcess class.

        Args:
            - model (KratosMultiphysics.Model): The model containing the sub_model_part to output.
            - params (KratosMultiphysics.Parameters): The parameters for the process, including:
                - "sub_model_part_name": Name of the sub_model_part to output.
                - "output_file_name": Name of the output json file.
                - "output_variables": List of variables to output.
                - "gauss_points_output_variables": List of gauss point variables to output.
                - "frequency": Frequency of output in seconds.
                - "historical_value": Whether to output historical values or not.
                - "resultant_solution": Whether to compute a resultant solution or not.
        """

        super().__init__(model, params)

    def ExecuteFinalizeSolutionStep(self):
        """
        Finalize the solution step by writing the output to a json file.
        """

        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:

            data = fast_read_external_json(self.output_file_name)

            self.time_counter = 0.0
            data["TIME"].append(time)

            # Nodal values
            if self.output_variables:
                count = 0
                for node in self.sub_model_part.Nodes:
                    compute = self._JsonOutputProcess__check_flag(node)

                    if compute:
                        node_identifier = "NODE_" + self._JsonOutputProcess__get_node_identifier(node)

                        for variable in self.output_variables:
                            variable_name = variable.Name()
                            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                            if self.historical_value:
                                value = node.GetSolutionStepValue(variable, 0)
                            else:
                                value = node.GetValue(variable)

                            if variable_type == "Double":
                                if not self.resultant_solution:
                                    data[node_identifier][variable_name].append(value)
                                else:
                                    if count == 0:
                                        data["RESULTANT"][variable_name].append(value)
                                    else:
                                        data["RESULTANT"][variable_name][-1] += value
                            elif variable_type == "Array":
                                if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                                    if not self.resultant_solution:
                                        data[node_identifier][variable_name + "_X"].append(value[0])
                                        data[node_identifier][variable_name + "_Y"].append(value[1])
                                        data[node_identifier][variable_name + "_Z"].append(value[2])
                                    else:
                                        if count == 0:
                                            data["RESULTANT"][variable_name + "_X"].append(value[0])
                                            data["RESULTANT"][variable_name + "_Y"].append(value[1])
                                            data["RESULTANT"][variable_name + "_Z"].append(value[2])
                                        else:
                                            data["RESULTANT"][variable_name + "_X"][-1] += value[0]
                                            data["RESULTANT"][variable_name + "_Y"][-1] += value[1]
                                            data["RESULTANT"][variable_name + "_Z"][-1] += value[2]
                                else:
                                    if not self.resultant_solution:
                                        list = self._JsonOutputProcess__kratos_vector_to_python_list(value)
                                        data[node_identifier][variable_name].append(list)
                                    else:
                                        aux = 0.0
                                        for index in range(len(value)):
                                            aux += value[index]
                                        if count == 0:
                                            data["RESULTANT"][variable_name].append(aux)
                                        else:
                                            data["RESULTANT"][variable_name][-1] += aux
                            elif variable_type == "Vector":
                                if not self.resultant_solution:
                                    data[node_identifier][variable_name].append(self._JsonOutputProcess__kratos_vector_to_python_list(value))
                                else:
                                    if count == 0:
                                        data["RESULTANT"][variable_name].append(self._JsonOutputProcess__kratos_vector_to_python_list(value))
                                    else:
                                        data["RESULTANT"][variable_name][-1] += self._JsonOutputProcess__kratos_vector_to_python_list(value)

                            # TODO: Add pending classes
                        count += 1

            count = 0
            # Gauss points values
            if self.gauss_points_output_variables:
                for elem in self.sub_model_part.Elements:
                    compute = self._JsonOutputProcess__check_flag(elem)

                    if compute:
                        for variable in self.gauss_points_output_variables:
                            variable_name = variable.Name()
                            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                            value = elem.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)

                            gauss_point_number = len(value)

                            if variable_type == "Double":
                                if not self.resultant_solution:
                                    for gp in range(gauss_point_number):
                                        data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(value[gp])
                                else:
                                    if count == 0:
                                        for gp in range(gauss_point_number):
                                            data["RESULTANT"][variable_name][str(gp)].append(value[gp])
                                    else:
                                        for gp in range(gauss_point_number):
                                            data["RESULTANT"][variable_name][str(gp)][-1] += value[gp]
                            elif variable_type == "Array":
                                if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                                    if not self.resultant_solution:
                                        for gp in range(gauss_point_number):
                                            data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)].append(
                                                value[gp][0])
                                            data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"][str(gp)].append(
                                                value[gp][1])
                                            data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)].append(
                                                value[gp][2])
                                    else:
                                        if count == 0:
                                            for gp in range(gauss_point_number):
                                                data["RESULTANT"][variable_name + "_X"][str(gp)].append(value[gp][0])
                                                data["RESULTANT"][variable_name + "_Y"][str(gp)].append(value[gp][1])
                                                data["RESULTANT"][variable_name + "_Z"][str(gp)].append(value[gp][2])
                                        else:
                                            for gp in range(gauss_point_number):
                                                data["RESULTANT"][variable_name + "_X"][str(gp)][-1] += value[gp][0]
                                                data["RESULTANT"][variable_name + "_Y"][str(gp)][-1] += value[gp][1]
                                                data["RESULTANT"][variable_name + "_Z"][str(gp)][-1] += value[gp][2]
                                else:
                                    if not self.resultant_solution:
                                        list = self._JsonOutputProcess__kratos_vector_to_python_list(value)
                                        for gp in range(gauss_point_number):
                                            data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                                    else:
                                        if count == 0:
                                            for gp in range(gauss_point_number):
                                                aux = 0.0
                                                for index in range(len(value[gp])):
                                                    aux += value[index]
                                                data["RESULTANT"][variable_name][str(gp)].append(aux)
                                        else:
                                            for gp in range(gauss_point_number):
                                                aux = 0.0
                                                for index in range(len(value[gp])):
                                                    aux += value[index]
                                                data["RESULTANT"][variable_name][str(gp)][-1] += aux
                            elif variable_type == "Vector":
                                if not self.resultant_solution:
                                    for gp in range(gauss_point_number):
                                        list = self._JsonOutputProcess__kratos_vector_to_python_list(value[gp])
                                        data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                                else:
                                    if count == 0:
                                        for gp in range(gauss_point_number):
                                            list = self._JsonOutputProcess__kratos_vector_to_python_list(value[gp])
                                            data["RESULTANT"][variable_name][str(gp)][-1] += list

                                # TODO: Add pending classes
                    count += 1

            fast_write_external_json(self.output_file_name, data)


def Factory(settings: KratosMultiphysics.Parameters, Model: KratosMultiphysics.Model):
    """
    Factory method to create an instance of the FastJsonOutputProcess.

    Args:
        - settings (KratosMultiphysics.Parameters): The parameters for the process
        - Model (KratosMultiphysics.Model): The model containing the sub_model_part to output.
    """
    return FastJsonOutputProcess(Model, settings["Parameters"])