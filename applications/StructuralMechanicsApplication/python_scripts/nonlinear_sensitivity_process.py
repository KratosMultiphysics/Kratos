from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import json
import math
import numpy as np
from json_utilities import *
KratosMultiphysics.CheckForPreviousImport()

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NonlinearSensitivityProcess(Model, settings["Parameters"])

class NonlinearSensitivityProcess(KratosMultiphysics.Process):
    """This class is used to compute sensitivity information to charactarize
    to structural nonlinear behaviour.

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
            "help"                 : "This process is used to compute sensitivity information to charactarize to structural nonlinear behaviour.",
            "check_variables"      : [],
            "gauss_points_check_variables" : [],
            "input_file_name"      : "",
            "model_part_name"      : "",
            "sub_model_part_name"  : "",
            "check_for_flag"       : "",
            "historical_value"     : true,
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

        self.input_file_name = self.params["input_file_name"].GetString()
        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.gauss_points_check_variables = self.__generate_variable_list_from_input(self.params["gauss_points_check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.historical_value = self.params["historical_value"].GetBool()
        #self.data =  read_external_json(input_file_name)

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
        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.data =  read_external_json(self.input_file_name)
        #time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        #dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        input_time_list = self.data["TIME"]

        if len(input_time_list) is 3:
            print("Hurra")

            # Nodal values
            for node in self.sub_model_part.Nodes:
                compute = self.__check_flag(node)

                if (compute == True):
                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        #if (self.historical_value == True):
                        #    value = node.GetSolutionStepValue(variable, 0)
                        #else:
                        #    value = node.GetValue(variable)

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            values_json = self.data["NODE_" + str(node.Id)][variable_name]
                            #value_json = self.__linear_interpolation(time, input_time_list, values_json)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                # X-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_X"]
                                curvature_x = self.__ComputeEFCurvature(values_json, input_time_list)
                                sen_first_x, sen_second_x = self.__ComputeFirstAndSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                if variable_name == "DISPLACEMENT":
                                    node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_X, curvature_x)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_X, sen_first_x)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_X, sen_second_x)
                                # Y-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Y"]
                                curvature_y = self.__ComputeEFCurvature(values_json, input_time_list)
                                sen_first_y, sen_second_y = self.__ComputeFirstAndSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                if variable_name == "DISPLACEMENT":
                                    node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Y, curvature_y)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Y, sen_first_y)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Y, sen_second_y)
                                # Z-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Z"]
                                curvature_z = self.__ComputeEFCurvature(values_json, input_time_list)
                                sen_first_z, sen_second_z = self.__ComputeFirstAndSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                if variable_name == "DISPLACEMENT":
                                    node.SetValue(StructuralMechanicsApplication.DISPLACEMENT_NL_SENSITIVITY_Z, curvature_z)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_FIRST_ORDER_Z, sen_first_z)
                                    node.SetValue(StructuralMechanicsApplication.NL_SENSITIVITY_SECOND_ORDER_Z, sen_second_z)
                            else:
                                values_json = self.data["NODE_"+str(node.Id)][variable_name][step - 1]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data["NODE_"+str(node.Id)][variable_name][step - 1]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])

            # Elemental values
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
                                #value_json = self.__linear_interpolation(time, input_time_list, values_json)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                for gp in range(gauss_point_number):
                                    # X-component
                                    values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)]
                                    print("************************************************************************")
                                    curvature = self.__ComputeEFCurvature(values_json, input_time_list)
                                    sen_first_x, sen_second_x = self.__ComputeFirstAndSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                    #print(curvature)
                                    #print(sen_first_x)
                                    #print(sen_second_x)
                                    print("************************************************************************")
                                    #value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    # Y-component
                                    values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Y"][str(gp)]
                                    #value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    # Z-component
                                    values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Z"][str(gp)]
                                    #value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            else:
                                for gp in range(gauss_point_number):
                                    values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][step - 1]
                                    for index in range(len(value[gp])):
                                        value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                        # Vector variable
                        elif variable_type == "Vector":
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)][step - 1]
                                for index in range(len(value[gp])):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])


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

    def __ComputeEFCurvature(self, response_value_array, load_factor_array):
        polynom = np.polyfit(load_factor_array, response_value_array, 2)
        curvature = 2 * polynom[0]
        return curvature

    def __ComputeFirstAndSecondOrderNLSensitivityFactors(self, response_value_array, load_factor_array):
        lambda_0 = load_factor_array[0]
        lambda_1 = load_factor_array[1]
        lambda_2 = load_factor_array[2]
        f_1 = lambda_1 / lambda_0
        delta_10 = lambda_1 - lambda_0
        delta_20 = lambda_2 - lambda_0

        response_0 = response_value_array[0]
        response_1 = response_value_array[1]
        response_2 = response_value_array[2]

        sensitivity_first_order = 0.0
        sensitivity_second_order = 0.0

        if abs(response_0) > 1e-8:
            sensitivity_first_order = response_1 / ( response_0 * f_1 )
            slope_10 = (response_1 - response_0) / delta_10
            slope_20 = (response_2 - response_1) / delta_20
            sensitivity_second_order = slope_20 / slope_10

        return sensitivity_first_order, sensitivity_second_order
