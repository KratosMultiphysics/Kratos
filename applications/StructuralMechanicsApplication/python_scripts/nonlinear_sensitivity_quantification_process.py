from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import sys
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics import json_utilities
import json
import math
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NonlinearSensitivityQuantificationProcess(model, settings["Parameters"])


def zero_vector(size):
    v = KratosMultiphysics.Vector(size)
    for i in range(size):
        v[i] = 0.0
    return v


def zero_matrix(row_size, column_size):
    m = KratosMultiphysics.Matrix(row_size, column_size)
    for row_index in range(row_size):
        for column_index in range(column_size):
            m[row_index, column_index] = 0.0
    return m


def ComputeEFCurvature(response_value_array, load_factor_array, take_absolute_value):
    if take_absolute_value == True:
        for i in range(len(response_value_array)):
            response_value_array[i] = abs(response_value_array[i])

    polynom = np.polyfit(load_factor_array, response_value_array, 2)
    curvature = 2 * polynom[0]
    return curvature


def ComputeFirstOrderNLSensitivityFactors(response_value_array, load_factor_array):
    lambda_0 = load_factor_array[0]
    lambda_1 = load_factor_array[1]
    f_1 = lambda_1 / lambda_0
    if f_1 < 1e-10:
        raise Exception("Pseudo-time steps are not valid!")

    response_0 = response_value_array[0]
    response_1 = response_value_array[1]

    sensitivity_first_order = 0.0

    if abs(response_0) > 1e-10:
        sensitivity_first_order = response_1 / ( response_0 * f_1 )

    return sensitivity_first_order


def ComputeSecondOrderNLSensitivityFactors(response_value_array, load_factor_array):
    lambda_0 = load_factor_array[0]
    lambda_1 = load_factor_array[1]
    lambda_2 = load_factor_array[2]
    delta_10 = lambda_1 - lambda_0
    delta_20 = lambda_2 - lambda_0
    if (delta_10 < 1e-10) or (delta_20 < 1e-10):
        raise Exception("Pseudo-time steps are not valid!")

    response_0 = response_value_array[0]
    response_1 = response_value_array[1]
    response_2 = response_value_array[2]

    sensitivity_second_order = 0.0


    slope_10 = (response_1 - response_0) / delta_10
    slope_20 = (response_2 - response_0) / delta_20

    if abs(slope_10) > 1e-10:
        sensitivity_second_order = slope_20 / slope_10

    return sensitivity_second_order


def AssembleVectorValuesIntoMatrix(given_vector, row_size, column_size):
    vector_size = given_vector.Size()
    if not vector_size == (row_size*column_size):
        raise Exception("Vector size does not fit!")
    m = KratosMultiphysics.Matrix(row_size, column_size)
    index = 0
    for row_index in range(row_size):
        for column_index in range(column_size):
            m[row_index, column_index] = given_vector[index]
            index += 1
    return m

def AssembleResults(component, variable, result_list):
    output_variable_1 = KratosMultiphysics.KratosGlobals.GetVariable(variable.Name() + "_NL_SENSITIVITY")
    output_variable_2 = KratosMultiphysics.KratosGlobals.GetVariable(variable.Name() + "_NL_SENSITIVITY_FIRST_ORDER")
    output_variable_3 = KratosMultiphysics.KratosGlobals.GetVariable(variable.Name() + "_NL_SENSITIVITY_SECOND_ORDER")
    variable_list = ((output_variable_1, output_variable_2, output_variable_3))
    for var_i, res_i in zip(variable_list, result_list):
        component.SetValue(var_i, res_i)

class NonlinearSensitivityQuantificationProcess(KratosMultiphysics.Process):
    """
        This process is used to quantifiy the non-linear behaviour of given quantities of interest.
    """

    def __init__(self, model, parameter):
        KratosMultiphysics.Process.__init__(self)
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model_part -- the model part used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                          : "This process is used to quantifiy the non-linear behaviour of given quantities of interest",
            "nodal_variables"               : [],
            "gauss_point_variables"         : [],
            "file_name"                     : "",
            "model_part_name"               : "",
            "sensitivity_model_part_name"   : "",
            "check_for_flag"                : "",
            "traced_time_steps"             : [],
            "historical_value"              : true,
            "resultant_solution"            : false,
            "curvature_absolute_value"      : true
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        parameter.ValidateAndAssignDefaults(default_parameters)

        model_part_name = parameter["model_part_name"].GetString()
        sensitivity_model_part_name = parameter["sensitivity_model_part_name"].GetString()

        self.model_part = model[model_part_name]
        # Get sensitivity model part
        if (sensitivity_model_part_name != ""):
            self.sensitivity_model_part = model[model_part_name].GetSubModelPart(sensitivity_model_part_name)
        else:
            self.sensitivity_model_part = model[model_part_name]

        # If we consider any flag
        flag_name = parameter["check_for_flag"].GetString()
        if (flag_name != ""):
            self.flag = KratosMultiphysics.KratosGlobals.GetFlag(flag_name)
        else:
            self.flag = None

        self.file_name = parameter["file_name"].GetString()
        self.nodal_variables = self.__GenerateVariableListFromInput(parameter["nodal_variables"])
        self.gauss_points_variables = self.__GenerateVariableListFromInput(parameter["gauss_point_variables"])
        self.historical_value = parameter["historical_value"].GetBool()
        self.absolute_value = parameter["curvature_absolute_value"].GetBool() # defines if the curvature computation is done with absolute values
        self.resultant_solution = parameter["resultant_solution"].GetBool()
        self.traced_time_steps = [ parameter["traced_time_steps"][i].GetDouble() for i in range( 0, parameter["traced_time_steps"].size() ) ]
        self.sensitivities_computed = False
        if len(self.traced_time_steps) is not 3:
            raise Exception("You have to provide three time steps!")


    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        for variable in self.nodal_variables:
            if variable.Name() == "DISPLACEMENT":
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(DISPLACEMENT_NL_SENSITIVITY, zero_vector(3), self.model_part.Nodes)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(DISPLACEMENT_NL_SENSITIVITY_FIRST_ORDER, zero_vector(3), self.model_part.Nodes)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(DISPLACEMENT_NL_SENSITIVITY_SECOND_ORDER, zero_vector(3), self.model_part.Nodes)
        for variable in self.gauss_points_variables:
            if variable.Name() == "FORCE":
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(FORCE_NL_SENSITIVITY, zero_vector(3), self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(FORCE_NL_SENSITIVITY_FIRST_ORDER, zero_vector(3), self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(FORCE_NL_SENSITIVITY_SECOND_ORDER, zero_vector(3), self.model_part.Elements)
            elif variable.Name() == "SHELL_MOMENT_GLOBAL":
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(SHELL_MOMENT_GLOBAL_NL_SENSITIVITY, zero_matrix(3, 3), self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(SHELL_MOMENT_GLOBAL_NL_SENSITIVITY_FIRST_ORDER, zero_matrix(3, 3), self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(SHELL_MOMENT_GLOBAL_NL_SENSITIVITY_SECOND_ORDER, zero_matrix(3, 3), self.model_part.Elements)
            elif variable.Name() == "VON_MISES_STRESS":
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(VON_MISES_STRESS_NL_SENSITIVITY, 0.0, self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(VON_MISES_STRESS_NL_SENSITIVITY_FIRST_ORDER, 0.0, self.model_part.Elements)
                KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(VON_MISES_STRESS_NL_SENSITIVITY_SECOND_ORDER, 0.0, self.model_part.Elements)


    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        This step generates the structure of the dictionary

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data = {}
        data["TIME"] = []
        count = 0
        for node in self.sensitivity_model_part.Nodes:
            compute = self.__CheckFlag(node)

            if (compute == True):
                if (self.resultant_solution == False):
                    data["NODE_" + str(node.Id)] = {}
                else:
                    if (count == 0):
                        data["RESULTANT"] = {}

                for variable in self.nodal_variables:
                    variable_name = variable.Name()
                    variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                    if (self.historical_value == True):
                        value = node.GetSolutionStepValue(variable, 0)
                    else:
                        value = node.GetValue(variable)

                    if (variable_type == "Double" or variable_type == "Component"):
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][variable_name] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][variable_name] = []
                    elif variable_type == "Array":
                        if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                            if (self.resultant_solution == False):
                                data["NODE_" + str(node.Id)][variable_name + "_X"] = []
                                data["NODE_" + str(node.Id)][variable_name + "_Y"] = []
                                data["NODE_" + str(node.Id)][variable_name + "_Z"] = []
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name + "_X"] = []
                                    data["RESULTANT"][variable_name + "_Y"] = []
                                    data["RESULTANT"][variable_name + "_Z"] = []
                        else:
                            if (self.resultant_solution == False):
                                data["NODE_" + str(node.Id)][variable_name] = []
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name] = []
                    elif variable_type == "Vector":
                        if (self.resultant_solution == False):
                            data["NODE_" + str(node.Id)][variable_name] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][variable_name] = []
                    # TODO: Add pending classes
                count += 1

        count = 0
        # Gauss points values
        for elem in self.sensitivity_model_part.Elements:
            compute = self.__CheckFlag(elem)

            if (compute == True):
                if (self.resultant_solution == False):
                    data["ELEMENT_" + str(elem.Id)] = {}
                else:
                    data["RESULTANT"] = {}

                for variable in self.gauss_points_variables:
                    variable_name = variable.Name()
                    variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                    value = elem.CalculateOnIntegrationPoints(variable, self.sensitivity_model_part.ProcessInfo)

                    gauss_point_number = len(value)

                    if (variable_type == "Double" or variable_type == "Component"):
                        if (self.resultant_solution == False):
                            data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][variable_name][str(gp)] = []
                    elif variable_type == "Array":
                        if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                            if (self.resultant_solution == False):
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_X"] = {}
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"] = {}
                                data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"] = {}
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)] = []
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"][str(gp)] = []
                                    data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)] = []
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name + "_X"] = {}
                                    data["RESULTANT"][variable_name + "_Y"] = {}
                                    data["RESULTANT"][variable_name + "_Z"] = {}
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][variable_name + "_X"][str(gp)] = []
                        else:
                            if (self.resultant_solution == False):
                                data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name] = {}
                    elif variable_type == "Vector":
                        if (self.resultant_solution == False):
                            data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][variable_name][str(gp)] = []
                    elif variable_type == "Matrix":
                        if (self.resultant_solution == False):
                            data["ELEMENT_" + str(elem.Id)][variable_name] = {}
                            for gp in range(gauss_point_number):
                                data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)] = []
                        else:
                            if (count == 0):
                                data["RESULTANT"][variable_name] = {}
                                for gp in range(gauss_point_number):
                                    data["RESULTANT"][variable_name][str(gp)] = []
                    # TODO: Add pending classes
                count += 1

        json_utilities.write_external_json(self.file_name, data)


    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data =  json_utilities.read_external_json(self.file_name)

        time = self.sensitivity_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        #dt = self.sensitivity_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        #self.time_counter += dt
        if round(time,1) in self.traced_time_steps:
            #self.time_counter = 0.0
            data["TIME"].append(time)
            count = 0

            # Nodal values
            for node in self.sensitivity_model_part.Nodes:
                compute = self.__CheckFlag(node)

                if (compute == True):
                    for variable in self.nodal_variables:
                        variable_name = variable.Name()
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        if (self.historical_value == True):
                            value = node.GetSolutionStepValue(variable, 0)
                        else:
                            value = node.GetValue(variable)

                        if (variable_type == "Double" or variable_type == "Component"):
                            if (self.resultant_solution == False):
                                data["NODE_" + str(node.Id)][variable_name].append(value)
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name].append(value)
                                else:
                                    data["RESULTANT"][variable_name][-1] += value
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                if (self.resultant_solution == False):
                                    data["NODE_" + str(node.Id)][variable_name + "_X"].append(value[0])
                                    data["NODE_" + str(node.Id)][variable_name + "_Y"].append(value[1])
                                    data["NODE_" + str(node.Id)][variable_name + "_Z"].append(value[2])
                                else:
                                    if (count == 0):
                                        data["RESULTANT"][variable_name + "_X"].append(value[0])
                                        data["RESULTANT"][variable_name + "_Y"].append(value[1])
                                        data["RESULTANT"][variable_name + "_Z"].append(value[2])
                                    else:
                                        data["RESULTANT"][variable_name + "_X"][-1] += value[0]
                                        data["RESULTANT"][variable_name + "_Y"][-1] += value[1]
                                        data["RESULTANT"][variable_name + "_Z"][-1] += value[2]
                            else:
                                if (self.resultant_solution == False):
                                    list = self.__kratos_vector_to__python_list(value)
                                    data["NODE_" + str(node.Id)][variable_name ].append(list)
                                else:
                                    aux = 0.0
                                    for index in range(len(value)):
                                        aux += value[index]
                                    if (count == 0):
                                        data["RESULTANT"][variable_name ].append(aux)
                                    else:
                                        data["RESULTANT"][variable_name ][-1] += aux
                        elif variable_type == "Vector":
                            if (self.resultant_solution == False):
                                data["NODE_" + str(node.Id)][variable_name].append(value)
                            else:
                                if (count == 0):
                                    data["RESULTANT"][variable_name].append(value)
                                else:
                                    data["RESULTANT"][variable_name][-1] += value

                        # TODO: Add pending classes
                    count += 1

            count = 0
            # Gauss points values
            for elem in self.sensitivity_model_part.Elements:
                compute = self.__CheckFlag(elem)

                if (compute == True):
                    for variable in self.gauss_points_variables:
                        variable_name = variable.Name()
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        value = elem.CalculateOnIntegrationPoints(variable, self.sensitivity_model_part.ProcessInfo)

                        gauss_point_number = len(value)

                        if (variable_type == "Double" or variable_type == "Component"):
                            if (self.resultant_solution == False):
                                for gp in range(gauss_point_number):
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(value[gp])
                            else:
                                if (count == 0):
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][variable_name][str(gp)].append(value[gp])
                                else:
                                    for gp in range(gauss_point_number):
                                        data["RESULTANT"][variable_name][str(gp)][-1] += value[gp]
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                if (self.resultant_solution == False):
                                    for gp in range(gauss_point_number):
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)].append(value[gp][0])
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_Y"][str(gp)].append(value[gp][1])
                                        data["ELEMENT_" + str(elem.Id)][variable_name + "_Z"][str(gp)].append(value[gp][2])
                                else:
                                    if (count == 0):
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
                                if (self.resultant_solution == False):
                                    list = self.__kratos_vector_to__python_list(value)
                                    for gp in range(gauss_point_number):
                                        data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                                else:
                                    if (count == 0):
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
                            if (self.resultant_solution == False):
                                for gp in range(gauss_point_number):
                                    list = self.__kratos_vector_to__python_list(value[gp])
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                            else:
                                if (count == 0):
                                    for gp in range(gauss_point_number):
                                        list = self.__kratos_vector_to__python_list(value[gp])
                                        data["RESULTANT"][variable_name][str(gp)][-1] += list
                        elif variable_type == "Matrix":
                            if (self.resultant_solution == False):
                                for gp in range(gauss_point_number):
                                    list = self.__kratos_matrix_to__python_list(value[gp])
                                    data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)].append(list)
                count += 1

        json_utilities.write_external_json(self.file_name, data)


    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.data =  json_utilities.read_external_json(self.file_name)
        input_time_list = self.data["TIME"]
        time = self.sensitivity_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        if len(input_time_list) is 3 and self.sensitivities_computed is False and round(time,1) == self.traced_time_steps[2]:
            self.sensitivities_computed = True
            # Nodal values
            for node in self.sensitivity_model_part.Nodes:
                compute = self.__CheckFlag(node)
                curvature_array = KratosMultiphysics.Array3()
                sen_first_array = KratosMultiphysics.Array3()
                sen_second_array = KratosMultiphysics.Array3()
                if (compute == True):
                    for variable in self.nodal_variables:
                        variable_name = variable.Name()
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        result_list = []
                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            values_json = self.data["NODE_" + str(node.Id)][variable_name]
                            result_list.append(ComputeEFCurvature(values_json, input_time_list, self.absolute_value))
                            result_list.append(ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list))
                            result_list.append(ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list))

                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                # X-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_X"]
                                curvature_array[0] = ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                sen_first_array[0] = ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                sen_second_array[0] = ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                # Y-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Y"]
                                curvature_array[1] = ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                sen_first_array[1] = ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                sen_second_array[1] = ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                # Z-component
                                values_json = self.data["NODE_" + str(node.Id)][variable_name + "_Z"]
                                curvature_array[2] = ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                sen_first_array[2] = ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                sen_second_array[2] = ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)

                                result_list.extend([curvature_array, sen_first_array, sen_second_array])

                        AssembleResults(node, variable, result_list)

            # Elemental values
            for elem in self.sensitivity_model_part.Elements:
                compute = self.__CheckFlag(elem)
                curvature_array = KratosMultiphysics.Array3()
                sen_first_array = KratosMultiphysics.Array3()
                sen_second_array = KratosMultiphysics.Array3()
                if (compute == True):
                    for variable in self.gauss_points_variables:
                        variable_name = variable.Name()
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                        value = elem.CalculateOnIntegrationPoints(variable, self.sensitivity_model_part.ProcessInfo)

                        gauss_point_number = len(value)
                        result_list = []
                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name][str(gp)]
                                curvature = ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                sen_first = ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                sen_second = ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                #TODO test it not used so far!
                                curvature /= gauss_point_number
                                sen_first /= gauss_point_number
                                sen_second /= gauss_point_number
                                result_list.extend([curvature, sen_first, sen_second])

                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                for gp in range(gauss_point_number):
                                    # X-component
                                    values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name + "_X"][str(gp)]
                                    curvature_array[0] += ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                    sen_first_array[0] += ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                    sen_second_array[0] += ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                    # Y-component
                                    values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Y"][str(gp)]
                                    curvature_array[1] += ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                    sen_first_array[1] += ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                    sen_second_array[1] += ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                                    # Z-component
                                    values_json = self.data["ELEMENT_"+str(elem.Id)][variable_name + "_Z"][str(gp)]
                                    curvature_array[2] += ComputeEFCurvature(values_json, input_time_list, self.absolute_value)
                                    sen_first_array[2] += ComputeFirstOrderNLSensitivityFactors(values_json, input_time_list)
                                    sen_second_array[2] += ComputeSecondOrderNLSensitivityFactors(values_json, input_time_list)
                            curvature_array /= gauss_point_number
                            sen_first_array /= gauss_point_number
                            sen_second_array /= gauss_point_number
                            result_list.extend([curvature_array, sen_first_array, sen_second_array])

                        # Matrix variable
                        elif variable_type == "Matrix":
                            row_size = value[0].Size1()
                            column_size = value[0].Size2()
                            result_size = row_size*column_size
                            curvature_vector = zero_vector(result_size)
                            sen_first_vector = zero_vector(result_size)
                            sen_second_vector = zero_vector(result_size)
                            for gp in range(gauss_point_number):
                                values_json = self.data["ELEMENT_" + str(elem.Id)][variable_name][str(gp)]
                                if len(values_json) is 3:
                                    index = 0
                                    for value1, value2, value3 in zip(values_json[0], values_json[1], values_json[2]):
                                        values_array = [value1, value2, value3]
                                        curvature_vector[index] += ComputeEFCurvature(values_array, input_time_list, self.absolute_value)
                                        sen_first_vector[index] += ComputeFirstOrderNLSensitivityFactors(values_array, input_time_list)
                                        sen_second_vector[index] += ComputeSecondOrderNLSensitivityFactors(values_array, input_time_list)
                                        index += 1
                                else:
                                    raise Exception("Number of response values does not fit!")
                            curvature_vector /= gauss_point_number
                            sen_first_vector /= gauss_point_number
                            sen_second_vector /= gauss_point_number
                            curvature_matrix = AssembleVectorValuesIntoMatrix(curvature_vector, row_size, column_size)
                            sen_first_matrix = AssembleVectorValuesIntoMatrix(sen_first_vector, row_size, column_size)
                            sen_second_matrix = AssembleVectorValuesIntoMatrix(sen_second_vector, row_size, column_size)
                            result_list.extend([curvature_matrix, sen_first_matrix, sen_second_matrix])

                        AssembleResults(elem, variable, result_list)



    def ExecuteFinalize(self):
        if self.sensitivities_computed is False:
            raise Exception("Process was not able to compute non-linear sensitivities!")


    def __GenerateVariableListFromInput(self, parameter):
        """ Parse a list of variables from input.

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        # At least verify that the input is a string
        if not parameter.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( parameter[i].GetString() ) for i in range( 0, parameter.size() ) ]


    def __CheckFlag(self, component):
        """ Checks the flag over a component

        Keyword arguments:
        self -- It signifies an instance of a class.
        component -- The Kratos node or element to check
        """

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True






