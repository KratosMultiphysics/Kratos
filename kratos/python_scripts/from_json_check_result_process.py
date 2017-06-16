from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import json
#import numpy as np # This cannot be here, manually interpolated
from json_utilities import *
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def linear_interpolation(x, x_list, y_list):
    ind_inf = 0
    ind_sup = -1
    x_inf = x_list[ind_inf]
    x_sup = x_list[ind_sup]
    for i in range(len(x_list)):
        if x_list[i] <= x:
            ind_inf = i
            x_inf = x_list[ind_inf]
        if x_list[-(1 + i)] >= x:
            ind_sup = -(1 + i)
            x_sup = x_list[ind_sup]

    if (x_sup-x_inf == 0):
        y = y_list[ind_inf]
    else:
        prop_sup = (x-x_inf)/(x_sup-x_inf)
        prop_inf = 1.0 - prop_sup
        y = y_list[ind_inf] * prop_inf + y_list[ind_sup] * prop_sup

    return y

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FromJsonCheckResultProcess(Model, settings["Parameters"])

class FromJsonCheckResultProcess(KratosMultiphysics.Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : [],
            "input_file_name"      : "",
            "model_part_name"      : "",
            "sub_model_part_name"  : "",
            "tolerance"            : 1e-3,
            "time_frequency"       : 1.00
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model_part

        self.params = params

        self.check_variables = []
        self.frequency    = 0.0
        self.time_counter = 0.0
        self.data = {}

    def ExecuteInitialize(self):
        input_file_name = self.params["input_file_name"].GetString()
        if (len(self.params["sub_model_part_name"].GetString()) > 0):
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["sub_model_part_name"].GetString())
        else:
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.data =  read_external_json(input_file_name)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        tol = self.params["tolerance"].GetDouble()
        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]
            for node in self.sub_model_part.Nodes:
                for i in range(self.params["check_variables"].size()):
                    out = self.params["check_variables"][i]
                    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )
                    val = node.GetSolutionStepValue(variable, 0)
                    if isinstance(val,float):
                        value_scalar = True
                    else:
                        value_scalar = False

                    value = node.GetSolutionStepValue(variable, 0)
                    # Scalar variable
                    if value_scalar:
                        values_json = self.data["NODE_"+str(node.Id)][out.GetString() ]
                        value_json = linear_interpolation(time, input_time_list, values_json)
                        self.assertAlmostEqual(value, value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" results."), delta=tol)
                    # Vector variable
                    else:
                        # X-component
                        values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_X"]
                        value_json = linear_interpolation(time, input_time_list, values_json)
                        self.assertAlmostEqual(value[0], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" X-component results."), delta=tol)
                        # Y-component
                        values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Y"]
                        value_json = linear_interpolation(time, input_time_list, values_json)
                        self.assertAlmostEqual(value[1], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" Y-component results."), delta=tol)
                        # Z-component
                        values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Z"]
                        value_json = linear_interpolation(time, input_time_list, values_json)
                        self.assertAlmostEqual(value[2], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" Z-component results."), delta=tol)

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def __generate_variable_list_from_input(self,param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
      return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]
