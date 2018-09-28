from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import json
import math
from json_utilities import *
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FromJsonCheckResultProcess(Model, settings["Parameters"])

class FromJsonCheckResultProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):
    """This class is used in order to check results using a json file
    containing the solution a given model part with a certain frequency

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
            "time_frequency"       : 1.00
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        params.ValidateAndAssignDefaults(default_parameters)
        self.params = params

        self.model_part = model_part

        self.iscloseavailable = hasattr(math,  "isclose")

        self.check_variables = []
        self.frequency    = 0.0
        self.time_counter = 0.0
        self.data = {}

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We get the submodelpart
        model_part_name = self.params["model_part_name"].GetString()
        sub_model_part_name = self.params["sub_model_part_name"].GetString()
        if (sub_model_part_name != ""):
            self.sub_model_part = self.model_part[model_part_name].GetSubModelPart(sub_model_part_name)
        else:
            self.sub_model_part = self.model_part[model_part_name]

        # If we consider any flag
        flag_name = self.params["check_for_flag"].GetString()
        if (flag_name != ""):
            self.flag = KratosMultiphysics.KratosGlobals.GetFlag(flag_name)
        else:
            self.flag = None

        input_file_name = self.params["input_file_name"].GetString()
        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.gauss_points_check_variables = self.__generate_variable_list_from_input(self.params["gauss_points_check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.historical_value = self.params["historical_value"].GetBool()
        self.data =  read_external_json(input_file_name)

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        This step generates the structure of the dictionary

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

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

        tol = self.params["tolerance"].GetDouble()
        reltol = self.params["relative_tolerance"].GetDouble()
        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]

            # Nodal values
            for node in self.sub_model_part.Nodes:
                compute = self.__check_flag(node)

                if (compute == True):
                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        if (self.historical_value == True):
                            value = node.GetSolutionStepValue(variable, 0)
                        else:
                            value = node.GetValue(variable)

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            values_json = self.data["NODE_" + str(node.Id)][variable_name]
                            value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            if (self.iscloseavailable == True):
                                isclosethis = math.isclose(value, value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                            else:
                                self.assertAlmostEqual(value, value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " results."), delta=tol)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                if (self.iscloseavailable == True):
                                    # X-component
                                    values_json = self.data["NODE_" + str(node.Id)][variable_name + "_X"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[0], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                                    # Y-component
                                    values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Y"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[1], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != "+str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                                    # Z-component
                                    values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Z"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[2], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value)+" != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                                else:
                                    # X-component
                                    values_json = self.data["NODE_"+str(node.Id)][variable_name  + "_X"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[0], value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " X-component results."), delta=tol)
                                    # Y-component
                                    values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Y"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[1], value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " Y-component results."), delta=tol)
                                    # Z-component
                                    values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Z"]
                                    value_json = self._linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[2], value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " Z-component results."), delta=tol)
                            else:
                                values_json = self.data["NODE_"+str(node.Id)][variable_name][step - 1]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    if (self.iscloseavailable == True):
                                        isclosethis = math.isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                                    else:
                                        self.assertAlmostEqual(value[index], value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " results."), delta=tol)
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data["NODE_"+str(node.Id)][variable_name][step - 1]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                if (self.iscloseavailable == True):
                                    isclosethis = math.isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking node " + str(node.Id) + " " + variable_name + " results."))
                                else:
                                    self.assertAlmostEqual(value[index], value_json, msg=("Error checking node " + str(node.Id) + " " + variable_name + " results."), delta=tol)
            # Nodal values
            for elem in self.sub_model_part.Elements:
                compute = self.__check_flag(elem)

                if (compute == True):
                    for i in range(self.params["gauss_points_check_variables"].size()):
                        out = self.params["gauss_points_check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        value = elem.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)

                        gauss_point_number = len(value)

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name][str(gp)]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                if (self.iscloseavailable == True):
                                    isclosethis = math.isclose(value[gp], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem " + str(elem.Id) + " " + variable_name + " results."))
                                else:
                                    self.assertAlmostEqual(value[gp], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " results."), delta=tol)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                if (self.iscloseavailable == True):
                                    for gp in range(gauss_point_number):
                                        # X-component
                                        values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        isclosethis = math.isclose(value[gp][0], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem " + str(elem.Id) + " " + variable_name + " results."))
                                        # Y-component
                                        values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Y"][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        isclosethis = math.isclose(value[gp][1], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem " + str(elem.Id) + " " + variable_name + " results."))
                                        # Z-component
                                        values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Z"][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        isclosethis = math.isclose(value[gp][2], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem "+str(elem.Id) + " " + variable_name + " results."))
                                else:
                                    for gp in range(gauss_point_number):
                                        # X-component
                                        values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name  + "_X"][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        self.assertAlmostEqual(value[gp][0], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " X-component results."), delta=tol)
                                        # Y-component
                                        values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name  + "_Y"][str(gp)]
                                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                        self.assertAlmostEqual(value[gp][1], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " Y-component results."), delta=tol)
                                        # Z-component
                                        values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)]
                                        value_json = self._linear_interpolation(time, input_time_list, values_json)
                                        self.assertAlmostEqual(value[gp][2], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " Z-component results."), delta=tol)
                            else:
                                for gp in range(gauss_point_number):
                                    values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][step - 1]
                                    for index in range(len(value[gp])):
                                        value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                        if (self.iscloseavailable == True):
                                            isclosethis = math.isclose(value[gp][index], value_json, rel_tol=reltol, abs_tol=tol)
                                            self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem " + str(elem.Id) + " " + variable_name + " results."))
                                        else:
                                            self.assertAlmostEqual(value[gp][index], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " results."), delta=tol)
                        # Vector variable
                        elif variable_type == "Vector":
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][step - 1]
                                for index in range(len(value[gp])):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    if (self.iscloseavailable == True):
                                        isclosethis = math.isclose(value[gp][index], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value[gp]) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking elem " + str(elem.Id) + " " + variable_name + " results."))
                                    else:
                                        self.assertAlmostEqual(value[gp][index], value_json, msg=("Error checking elem " + str(elem.Id) + " " + variable_name + " results."), delta=tol)

                        # TODO: Add pending classes

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

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True
