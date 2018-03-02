from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
from json_utilities import *
import json
KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return JsonOutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "python_process"

class JsonOutputProcess(KratosMultiphysics.Process):
    """This class is used in order to create a json file containing 
    the solution a given model part with a certain frequency

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model_part -- the model part used to construct the process.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model_part, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model_part -- the model part used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "output_variables"              : [],
            "gauss_points_output_variables" : [],
            "output_file_name"              : "",
            "model_part_name"               : "",
            "sub_model_part_name"           : "",
            "time_frequency"                : 1.00,
            "historical_value"              : true,
            "resultant_solution"            : false
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model_part

        self.params = params

        self.output_file_name = ""
        self.output_variables = []
        self.gauss_points_output_variables = []
        self.frequency = 0.0
        self.time_counter = 0.0
        self.resultant_solution = False

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        self.output_file_name = self.params["output_file_name"].GetString()
        if (len(self.params["sub_model_part_name"].GetString()) > 0):
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["sub_model_part_name"].GetString())
        else:
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]

        self.output_variables = self.__generate_variable_list_from_input(self.params["output_variables"])
        self.gauss_points_output_variables = self.__generate_variable_list_from_input(self.params["gauss_points_output_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.resultant_solution = self.params["resultant_solution"].GetBool()
        self.historical_value = self.params["historical_value"].GetBool()

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
            if (self.resultant_solution == False):
                data["NODE_" + str(node.Id)] = {}
            else:
                data["RESULTANT"] = {}

            for i in range(self.params["output_variables"].size()):
                out = self.params["output_variables"][i]
                variable_name = out.GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                variable_type = self.__check_variable_type(variable_name)
                if (self.historical_value == True):
                    value = node.GetSolutionStepValue(variable, 0)
                else:
                    value = node.GetValue(variable)
                    
                if (variable_type == "Double" or variable_type == "Component"):
                    if (self.resultant_solution == False):
                        data["NODE_" + str(node.Id)][out.GetString()] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString()] = []
                elif variable_type == "Array":
                    if (self.__check_variable_type(variable_name + "_X") == "Component"):
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][out.GetString() + "_X"] = []
                            data["NODE_" + str(node.Id)][out.GetString() + "_Y"] = []
                            data["NODE_" + str(node.Id)][out.GetString() + "_Z"] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString() + "_X"] = []
                                data["RESULTANT"][out.GetString() + "_Y"] = []
                                data["RESULTANT"][out.GetString() + "_Z"] = []
                    else:
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][out.GetString()] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString()] = []
                elif variable_type == "Vector":
                    if (self.resultant_solution == False):
                        data["NODE_" + str(node.Id)][out.GetString()] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString()] = []

                # TODO: Add pending classes

        # Gauss points values
        for elem in self.sub_model_part.Elements:
            if (self.resultant_solution == False):
                data["ELEMENT_" + str(elem.Id)] = {}
            else:
                data["RESULTANT"] = {}

            for i in range(self.params["gauss_points_output_variables"].size()):
                out = self.params["gauss_points_output_variables"][i]
                variable_name = out.GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                variable_type = self.__check_variable_type(variable_name)

                value = elem.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)

                gauss_point_number = len(value)

                if (variable_type == "Double" or variable_type == "Component"):
                    if (self.resultant_solution == False):
                        data["ELEMENT_" + str(elem.Id)][out.GetString()] = {}
                        for gp in range(gauss_point_number):
                            data["ELEMENT_" + str(elem.Id)][out.GetString()][gp] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString()] = {}
                            for gp in range(gauss_point_number):
                                data["RESULTANT"][out.GetString()][gp] = []
                elif variable_type == "Array":
                    if (self.__check_variable_type(variable_name + "_X") == "Component"):
                        if (self.resultant_solution == False):
                            data["ELEMENT_" + str(elem.Id)][out.GetString() + "_X"] = {}
                            data["ELEMENT_" + str(elem.Id)][out.GetString() + "_Y"] = {}
                            data["ELEMENT_" + str(node.Id)][out.GetString() + "_Z"] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][out.GetString() + "_X"][gp] = []
                                data["ELEMENT_" + str(elem.Id)][out.GetString() + "_Y"][gp] = []
                                data["ELEMENT_" + str(elem.Id)][out.GetString() + "_Z"][gp] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString() + "_X"] = {}
                                data["RESULTANT"][out.GetString() + "_Y"] = {}
                                data["RESULTANT"][out.GetString() + "_Z"] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][out.GetString() + "_X"][gp] = []
                    else:
                        if (self.resultant_solution == False):
                            data["ELEMENT_" + str(elem.Id)][out.GetString()] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][out.GetString()][gp] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString()] = {}
                elif variable_type == "Vector":
                    if (self.resultant_solution == False):
                        data["ELEMENT_" + str(elem.Id)][out.GetString()] = {}
                        for gp in range(gauss_point_number):
                            data["ELEMENT_" + str(elem.Id)][out.GetString()][gp] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString()] = {}
                            for gp in range(gauss_point_number):
                                data["RESULTANT"][out.GetString()][gp] = []

                # TODO: Add pending classes

            count += 1

        write_external_json(self.output_file_name, data)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data =  read_external_json(self.output_file_name)

        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            data["TIME"].append(time)
            count = 0

            # Nodal values
            for node in self.sub_model_part.Nodes:
                for i in range(self.params["output_variables"].size()):
                    out = self.params["output_variables"][i]
                    variable_name = out.GetString()
                    variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                    variable_type = self.__check_variable_type(variable_name)
                    if (self.historical_value == True):
                        value = node.GetSolutionStepValue(variable, 0)
                    else:
                        value = node.GetValue(variable)

                    if (variable_type == "Double" or variable_type == "Component"):
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][out.GetString()].append(value)
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString()].append(value)
                            else:
                                data["RESULTANT"][out.GetString()][-1] += value
                    elif variable_type == "Array":
                        if (self.__check_variable_type(variable_name + "_X") == "Component"):
                            if (self.resultant_solution == False):
                                data["NODE_" + str(node.Id)][out.GetString() + "_X"].append(value[0])
                                data["NODE_" + str(node.Id)][out.GetString() + "_Y"].append(value[1])
                                data["NODE_" + str(node.Id)][out.GetString() + "_Z"].append(value[2])
                            else:
                                if (count == 0):
                                    data["RESULTANT"][out.GetString() + "_X"].append(value[0])
                                    data["RESULTANT"][out.GetString() + "_Y"].append(value[1])
                                    data["RESULTANT"][out.GetString() + "_Z"].append(value[2])
                                else:
                                    data["RESULTANT"][out.GetString() + "_X"][-1] += value[0]
                                    data["RESULTANT"][out.GetString() + "_Y"][-1] += value[1]
                                    data["RESULTANT"][out.GetString() + "_Z"][-1] += value[2]
                        else:
                            if (self.resultant_solution == False):
                                list = self.__kratos_vector_to__python_list(value)
                                data["NODE_" + str(node.Id)][out.GetString() ].append(list)
                            else:
                                aux = 0.0
                                for index in range(len(value)):
                                    aux += value[index]
                                if (count == 0):
                                    data["RESULTANT"][out.GetString() ].append(aux)
                                else:
                                    data["RESULTANT"][out.GetString() ][-1] += aux
                    elif variable_type == "Vector":
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][out.GetString()].append(value)
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString()][-1] += value

                    # TODO: Add pending classes

            # Gauss points values
            for elem in self.sub_model_part.Elements:
                for i in range(self.params["output_variables"].size()):
                    out = self.params["output_variables"][i]
                    variable_name = out.GetString()
                    variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                    variable_type = self.__check_variable_type(variable_name)
                    value = elem.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)

                    gauss_point_number = len(value)

                    if (variable_type == "Double" or variable_type == "Component"):
                        if (self.resultant_solution == False):
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][out.GetString()][gp].append(value)
                        else:
                            if (count == 0):
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][out.GetString()][gp].append(value)
                            else:
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][out.GetString()][gp][-1] += value
                    elif variable_type == "Array":
                        if (self.__check_variable_type(variable_name + "_X") == "Component"):
                            if (self.resultant_solution == False):
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][out.GetString() + "_X"][gp].append(value[0])
                                    data["ELEMENT_" + str(elem.Id)][out.GetString() + "_Y"][gp].append(value[1])
                                    data["ELEMENT_" + str(elem.Id)][out.GetString() + "_Z"][gp].append(value[2])
                            else:
                                if (count == 0):
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][out.GetString() + "_X"][gp].append(value[0])
                                        data["RESULTANT"][out.GetString() + "_Y"][gp].append(value[1])
                                        data["RESULTANT"][out.GetString() + "_Z"][gp].append(value[2])
                                else:
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][out.GetString() + "_X"][gp][-1] += value[0]
                                        data["RESULTANT"][out.GetString() + "_Y"][gp][-1] += value[1]
                                        data["RESULTANT"][out.GetString() + "_Z"][gp][-1] += value[2]
                        else:
                            if (self.resultant_solution == False):
                                list = self.__kratos_vector_to__python_list(value)
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][out.GetString()][gp].append(list)
                            else:
                                aux = 0.0
                                for index in range(len(value)):
                                    aux += value[index]
                                if (count == 0):
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][out.GetString()][gp].append(aux)
                                else:
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][out.GetString()][gp][-1] += aux
                    elif variable_type == "Vector":
                        if (self.resultant_solution == False):
                            for gp in range(gauss_point_number):
                                data["NODE_" + str(node.Id)][out.GetString()][gp].append(value)
                        else:
                            if (count == 0):
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][out.GetString()][gp][-1] += value
                                    
                        # TODO: Add pending classes
                count += 1

        write_external_json(self.output_file_name, data)

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

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

    def __check_variable_type(self, variable):
        """ This method checks the type of variable to be saved

        Keyword arguments:
        self -- It signifies an instance of a class.
        variable -- The name of the variable to create
        """
        if KratosMultiphysics.HasBoolVariable(variable):
            return "Bool"
        elif KratosMultiphysics.HasIntVariable(variable):
            return "Integer"
        elif KratosMultiphysics.HasUnsignedIntVariable(variable):
            return "Unsigned Integer"
        elif KratosMultiphysics.HasDoubleVariable(variable):
            return "Double"
        elif KratosMultiphysics.HasArrayVariable(variable):
            return "Array"
        elif KratosMultiphysics.HasVectorVariable(variable):
            return "Vector"
        elif KratosMultiphysics.HasMatrixVariable(variable):
            return "Matrix"
        elif KratosMultiphysics.HasStringVariable(variable):
            return "String"
        elif KratosMultiphysics.HasVariableComponent(variable):
            return "Component"
        elif KratosMultiphysics.HasFlagsVariable(variable):
            return "Flag"
        else:
            return "NONE"


