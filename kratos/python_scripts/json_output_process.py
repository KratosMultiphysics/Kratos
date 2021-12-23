# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.json_utilities import read_external_json, write_external_json

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return JsonOutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class JsonOutputProcess(KratosMultiphysics.Process):
    """This class is used in order to create a json file containing
    the solution a given model part with a certain frequency

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model used to construct the process.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                          : "This process generates a json file containing the solution of a list of variables from a given submodelpart",
            "output_variables"              : [],
            "gauss_points_output_variables" : [],
            "output_file_name"              : "",
            "model_part_name"               : "",
            "sub_model_part_name"           : "",
            "check_for_flag"                : "",
            "time_frequency"                : 1.00,
            "historical_value"              : true,
            "resultant_solution"            : false,
            "use_node_coordinates"          : false
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.params = params

        self.time_counter = 0.0

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We get the submodelpart
        model_part_name = self.params["model_part_name"].GetString()
        sub_model_part_name = self.params["sub_model_part_name"].GetString()
        if sub_model_part_name != "":
            self.sub_model_part = self.model.GetModelPart(model_part_name).GetSubModelPart(sub_model_part_name)
        else:
            self.sub_model_part = self.model.GetModelPart(model_part_name)

        if self.sub_model_part.GetCommunicator().TotalProcesses() > 1: # mpi-execution
            raise Exception("This process cannot be used for writing output in MPI!")

        # If we consider any flag
        flag_name = self.params["check_for_flag"].GetString()
        if flag_name != "":
            self.flag = KratosMultiphysics.KratosGlobals.GetFlag(flag_name)
        else:
            self.flag = None

        self.output_file_name = self.params["output_file_name"].GetString()
        self.output_variables = self.__generate_variable_list_from_input(self.params["output_variables"])
        self.gauss_points_output_variables = self.__generate_variable_list_from_input(self.params["gauss_points_output_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.resultant_solution = self.params["resultant_solution"].GetBool()
        self.historical_value = self.params["historical_value"].GetBool()
        self.use_node_coordinates = self.params["use_node_coordinates"].GetBool()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        This step generates the structure of the dictionary

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data = {}
        data["TIME"] = []
        count = 0
        for node in self.sub_model_part.Nodes:
            compute = self.__check_flag(node)

            if compute:
                node_identifier = "NODE_" + self.__get_node_identifier(node)

                if not self.resultant_solution:
                    data[node_identifier] = {}
                else:
                    if count == 0:
                        data["RESULTANT"] = {}

                for i in range(self.params["output_variables"].size()):
                    out = self.params["output_variables"][i]
                    variable_name = out.GetString()
                    variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                    variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                    if self.historical_value:
                        value = node.GetSolutionStepValue(variable, 0)
                    else:
                        value = node.GetValue(variable)

                    if variable_type == "Double":
                        if not self.resultant_solution:
                            data[node_identifier][variable_name] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = []
                    elif variable_type == "Array":
                        if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                            if not self.resultant_solution:
                                data[node_identifier][variable_name + "_X"] = []
                                data[node_identifier][variable_name + "_Y"] = []
                                data[node_identifier][variable_name + "_Z"] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name + "_X"] = []
                                    data["RESULTANT"][variable_name + "_Y"] = []
                                    data["RESULTANT"][variable_name + "_Z"] = []
                        else:
                            if not self.resultant_solution:
                                data[node_identifier][variable_name] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name] = []
                    elif variable_type == "Vector":
                        if not self.resultant_solution:
                            data[node_identifier][variable_name] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = []
                    # TODO: Add pending classes
                count += 1

        count = 0
        # Gauss points values
        for elem in self.sub_model_part.Elements:
            compute = self.__check_flag(elem)

            if compute:
                if not self.resultant_solution:
                    data["ELEMENT_" + str(elem.Id)] = {}
                else:
                    data["RESULTANT"] = {}

                for i in range(self.params["gauss_points_output_variables"].size()):
                    out = self.params["gauss_points_output_variables"][i]
                    variable_name = out.GetString()
                    variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                    variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                    value = elem.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)

                    gauss_point_number = len(value)

                    if variable_type == "Double":
                        if not self.resultant_solution:
                            data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][variable_name][str(gp)] = []
                    elif variable_type == "Array":
                        if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                            if not self.resultant_solution:
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_X"] = {}
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"] = {}
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"] = {}
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)] = []
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"][str(gp)] = []
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name + "_X"] = {}
                                    data["RESULTANT"][variable_name + "_Y"] = {}
                                    data["RESULTANT"][variable_name + "_Z"] = {}
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][variable_name + "_X"][str(gp)] = []
                        else:
                            if not self.resultant_solution:
                                data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name] = {}
                    elif variable_type == "Vector":
                        if not self.resultant_solution:
                            data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][variable_name][str(gp)] = []
                    # TODO: Add pending classes
                count += 1

        write_external_json(self.output_file_name, data)

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data = read_external_json(self.output_file_name)

        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            data["TIME"].append(time)
            count = 0

            # Nodal values
            for node in self.sub_model_part.Nodes:
                compute = self.__check_flag(node)

                if compute:
                    node_identifier = "NODE_" + self.__get_node_identifier(node)

                    for i in range(self.params["output_variables"].size()):
                        out = self.params["output_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
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
                                    list = self.__kratos_vector_to__python_list(value)
                                    data[node_identifier][variable_name ].append(list)
                                else:
                                    aux = 0.0
                                    for index in range(len(value)):
                                        aux += value[index]
                                    if count == 0:
                                        data["RESULTANT"][variable_name ].append(aux)
                                    else:
                                        data["RESULTANT"][variable_name ][-1] += aux
                        elif variable_type == "Vector":
                            if not self.resultant_solution:
                                data[node_identifier][variable_name].append(value)
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name].append(value)
                                else:
                                    data["RESULTANT"][variable_name][-1] += value

                        # TODO: Add pending classes
                    count += 1

            count = 0
            # Gauss points values
            for elem in self.sub_model_part.Elements:
                compute = self.__check_flag(elem)

                if compute:
                    for i in range(self.params["gauss_points_output_variables"].size()):
                        out = self.params["gauss_points_output_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
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
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)].append(value[gp][0])
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"][str(gp)].append(value[gp][1])
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)].append(value[gp][2])
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
                                    list = self.__kratos_vector_to__python_list(value)
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
                                    list = self.__kratos_vector_to__python_list(value[gp])
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                            else:
                                if count == 0:
                                    for gp in range(gauss_point_number):
                                        list = self.__kratos_vector_to__python_list(value[gp])
                                        data["RESULTANT"][variable_name][str(gp)][-1] += list

                            # TODO: Add pending classes
                count += 1

        write_external_json(self.output_file_name, data)

    def __kratos_vector_to__python_list(self, value):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        list = []
        for index in range(len(value)):
            list.append(value[index])
        return list

    def __generate_variable_list_from_input(self, param):
        """ Parse a list of variables from input.

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [KratosMultiphysics.KratosGlobals.GetVariable(param[i].GetString()) for i in range( 0, param.size())]

    def __check_flag(self, component):
        """ Checks the flag over a component

        Keyword arguments:
        self -- It signifies an instance of a class.
        component -- The Kratos node or element to check
        """

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True

    def __get_node_identifier(self, node):
        """ returns the identifier/key for saving nodal results in the json
        this can be either the node Id or its coordinates
        The coordinates can be used to check the nodal results in MPI

        Keyword arguments:
        self -- It signifies an instance of a class.
        node -- The Kratos node to get the identifier for
        """
        if self.use_node_coordinates:
            return 'X_{0:.{digits}f}_Y_{1:.{digits}f}_Z_{2:.{digits}f}'.format(node.X0, node.Y0, node.Z0, digits=6)
        else:
            return str(node.Id)
