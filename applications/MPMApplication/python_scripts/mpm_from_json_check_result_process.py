# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.json_utilities import read_external_json

# Import KratosUnittest
from KratosMultiphysics.KratosUnittest import isclose as t_isclose
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Math import
from math import log10, ceil

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return MPMFromJsonCheckResultProcess(Model, settings["Parameters"])

class LegacyFromJsonCheckResultProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This class is used in order to check results using a json file
    containing the solution a given model part with a certain frequency

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model contaning the model_parts
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the model contaning the model_parts
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""{
            "help"                 : "This process checks the solution obtained from a given json file. It can be used for generating tests for a problem",
            "check_variables"      : [],
            "gauss_points_check_variables" : [],
            "input_file_name"      : "",
            "model_part_name"      : "",
            "sub_model_part_name"  : "",
            "check_for_flag"       : "",
            "historical_value"     : true,
            "tolerance"            : 1e-3,
            "relative_tolerance"   : 1e-6,
            "time_frequency"       : 1.00,
            "use_node_coordinates" : false,
            "check_only_local_entities" : false
        }""")

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)
        self.params = params
        self.model  = model

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        # We get the submodelpart
        model_part_name = self.params["model_part_name"].GetString()
        sub_model_part_name = self.params["sub_model_part_name"].GetString()
        if (sub_model_part_name != ""):
            self.model_part = self.model[model_part_name].GetSubModelPart(sub_model_part_name)
        else:
            self.model_part = self.model[model_part_name]

        # If we consider any flag
        flag_name = self.params["check_for_flag"].GetString()
        if flag_name != "":
            self.flag = KratosMultiphysics.KratosGlobals.GetFlag(flag_name)
        else:
            self.flag = None

        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.gauss_points_check_variables = self.__generate_variable_list_from_input(self.params["gauss_points_check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.historical_value = self.params["historical_value"].GetBool()
        self.data = read_external_json(self.params["input_file_name"].GetString())
        self.abs_tol = self.params["tolerance"].GetDouble()
        self.rel_tol = self.params["relative_tolerance"].GetDouble()
        self.use_node_coordinates = self.params["use_node_coordinates"].GetBool()
        self.check_only_local_entities = self.params["check_only_local_entities"].GetBool()
        self.rel_tol_digits = ComputeRelevantDigits(self.rel_tol)
        self.abs_tol_digits = ComputeRelevantDigits(self.abs_tol)

        # Initialize counter
        self.step_counter = 0
        self.time_counter = 0.0

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        time = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]

            # Nodal values
            for node in self.__get_nodes():
                compute = self.__check_flag(node)

                if compute:
                    node_identifier = "NODE_" + self.__get_node_identifier(node)

                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        if self.historical_value:
                            value = node.GetSolutionStepValue(variable, 0)
                        else:
                            value = node.GetValue(variable)

                        # Scalar variable
                        if variable_type == "Double":
                            values_json = self.data[node_identifier][variable_name]
                            value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            self.__check_values(node.Id, "Node", value, value_json, variable_name)
                        # Array variable
                        elif variable_type == "Array":
                            if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                                for component_index, component in enumerate(["_X", "_Y", "_Z"]):
                                    values_json = self.data[node_identifier][variable_name+component]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    self.__check_values(node.Id, "Node", value[component_index], value_json, variable_name+component)
                            else:
                                values_json = self.data[node_identifier][variable_name][self.step_counter]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    self.__check_values(node.Id, "Node", value[index], value_json, variable_name)
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data[node_identifier][variable_name][self.step_counter]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                self.__check_values(node.Id, "Node", value[index], value_json, variable_name)
            # Nodal values
            for elem in self.__get_elements():
                compute = self.__check_flag(elem)

                if compute is True:
                    for i in range(self.params["gauss_points_check_variables"].size()):
                        out = self.params["gauss_points_check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        value = elem.CalculateOnIntegrationPoints(variable, self.model_part.ProcessInfo)

                        gauss_point_number = len(value)

                        # Scalar variable
                        if variable_type == "Double":
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name][str(gp)]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                self.__check_values(elem.Id, "Element", value[gp], value_json, variable_name)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double"):
                                for gp in range(gauss_point_number):
                                    for component_index, component in enumerate(["_X", "_Y", "_Z"]):
                                        values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name+component][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        self.__check_values(elem.Id, "Element", value[gp][component_index], value_json, variable_name+component)
                            else:
                                for gp in range(gauss_point_number):
                                    values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][self.step_counter]
                                    for index in range(len(value[gp])):
                                        value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                        self.__check_values(elem.Id, "Element", value[gp][index], value_json, variable_name)
                        # Vector variable
                        elif variable_type == "Vector":
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][self.step_counter]
                                for index in range(len(value[gp])):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    self.__check_values(elem.Id, "Element", value[gp][index], value_json, variable_name)

                        # TODO: Add pending classes

            self.step_counter += 1

    def __linear_interpolation(self, x, x_list, y_list):
        """ This method is defined to interpolate values of a
        list using the PiecewiseLinearTable from Kratos

        Keyword arguments:
        self -- It signifies an instance of a class.
        x -- The value to interpolate
        x_list -- Values in X axis
        y_list -- Values in Y axis
        """

        tb = KratosMultiphysics.PiecewiseLinearTable()
        for i in range(len(x_list)):
            tb.AddRow(x_list[i], y_list[i])

        return tb.GetNearestValue(x)

    def __generate_variable_list_from_input(self,param):
        """ Parse a list of variables from input.

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]

    def __check_flag(self, component):
        """ Checks the flag over a component

        Keyword arguments:
        self -- It signifies an instance of a class.
        component -- The Kratos node or element to check
        """

        if self.flag is not None:
            if component.Is(self.flag) == False:
                return False

        return True

    def __check_values(self, entity_id, entity_type, value_entity, value_json, variable_name):
        """ Checks if two values are the same and issues a detailed error message
        in case they do not match up to the specified tolerance

        Keyword arguments:
        self -- It signifies an instance of a class.
        entity_id -- The Kratos node or element to check
        entity_type -- The type of the entity
        value_entity -- The value on the entity
        value_json -- The reference value from the json
        variable_name -- The name of the variable
        """
        relevant_digits = int(max(self.rel_tol_digits, self.abs_tol_digits))+1 # +1 for one more digit of output
        isclosethis = t_isclose(value_entity, value_json, rel_tol=self.rel_tol, abs_tol=self.abs_tol)
        msg  = 'Error checking {} #{} for variable {} results:\n'.format(entity_type, entity_id, variable_name)
        msg += '%.*f != %.*f; rel_tol=%.*f, abs_tol=%.*f' % (relevant_digits, value_entity, relevant_digits, value_json, self.rel_tol_digits, self.rel_tol, self.abs_tol_digits, self.abs_tol)
        self.assertTrue(isclosethis, msg=msg)

    def __get_node_identifier(self, node):
        """ returns the identifier/key for saving nodal results in the json
        this can be either the node Id or its coordinates
        The coordinates can be used to check the nodal results in MPI

        Keyword arguments:
        self -- It signifies an instance of a class.
        node -- The Kratos node to get the identifier for
        """
        if self.use_node_coordinates:
            digits = 6
            return 'X_%.*f_Y_%.*f_Z_%.*f' % (digits, node.X0, digits, node.Y0, digits, node.Z0)
        else:
            return str(node.Id)

    def __get_nodes(self):
        """ returns the nodes to be checked
        Either only local or all (local and ghost)
        This is ONLY relevant in MPI

        Keyword arguments:
        self -- It signifies an instance of a class.
        node -- The Kratos node to get the identifier for
        """
        if self.check_only_local_entities:
            return self.model_part.GetCommunicator().LocalMesh().Nodes
        else:
            return self.model_part.Nodes

    def __get_elements(self):
        """ returns the elements to be checked
        Either only local or all (local and ghost)
        This is ONLY relevant in MPI

        Keyword arguments:
        self -- It signifies an instance of a class.
        node -- The Kratos node to get the identifier for
        """
        if self.check_only_local_entities:
            return self.model_part.GetCommunicator().LocalMesh().Elements
        else:
            return self.model_part.Elements

def ComputeRelevantDigits(number):
    """ Computes the relevant digits

    Keyword arguments:
    self -- It signifies an instance of a class.
    """
    return int(ceil(abs(log10(number))))
class MPMFromJsonCheckResultProcess(LegacyFromJsonCheckResultProcess, KratosUnittest.TestCase): # TODO: This must be updated to the new C++ version

    def __init__(self, model_part, params):
        super(MPMFromJsonCheckResultProcess, self).__init__(model_part, params)

    def ExecuteFinalizeSolutionStep(self):

        tol = self.params["tolerance"].GetDouble()
        reltol = self.params["relative_tolerance"].GetDouble()
        time = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]

            # Material points values
            for mp in self.model_part.Elements:
                compute = self.__check_flag(mp)

                if (compute == True):
                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        value = mp.CalculateOnIntegrationPoints(variable,self.model_part.ProcessInfo)[0]

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Integer" or variable_type == "Component"):
                            values_json = self.data["MP_" + str(mp.Id)][variable_name]
                            value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            isclosethis = t_isclose(value, value_json, rel_tol=reltol, abs_tol=tol)
                            self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking material point " + str(mp.Id) + " " + variable_name + " results."))
                        # Array variable
                        elif variable_type == "Array":

                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double"):
                                for component_index, component in enumerate(["_X", "_Y", "_Z"]):
                                    values_json = self.data["MP_" + str(mp.Id)][variable_name +component]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = t_isclose(value[component_index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value[component_index]) + " != "+str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking material point " + str(mp.Id) + " " + variable_name + " results."))
                            else:
                                values_json = self.data["MP_"+str(mp.Id)][variable_name][step - 1]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    isclosethis = t_isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking material point " + str(mp.Id) + " " + variable_name + " results."))
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data["MP_"+str(mp.Id)][variable_name][step - 1]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                isclosethis = t_isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking material point " + str(mp.Id) + " " + variable_name + " results."))

    def __linear_interpolation(self, x, x_list, y_list):

        tb = KratosMultiphysics.PiecewiseLinearTable()
        for i in range(len(x_list)):
            tb.AddRow(x_list[i], y_list[i])

        return tb.GetNearestValue(x)

    def __generate_variable_list_from_input(self,param):

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]

    def __check_flag(self, component):

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True
