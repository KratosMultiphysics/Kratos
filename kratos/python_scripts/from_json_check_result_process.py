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

    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : [],
            "input_file_name"      : "",
            "model_part_name"      : "",
            "sub_model_part_name"  : "",
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
        input_file_name = self.params["input_file_name"].GetString()
        if (len(self.params["sub_model_part_name"].GetString()) > 0):
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["sub_model_part_name"].GetString())
        else:
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.historical_value = self.params["historical_value"].GetBool()
        self.data =  read_external_json(input_file_name)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        tol = self.params["tolerance"].GetDouble()
        reltol = self.params["relative_tolerance"].GetDouble()
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

                    if (self.historical_value == True):
                        value = node.GetSolutionStepValue(variable, 0)
                    else:
                        value = node.GetValue(variable)
                    # Scalar variable
                    if isinstance(value,float):
                        values_json = self.data["NODE_"+str(node.Id)][out.GetString() ]
                        value_json = self.__linear_interpolation(time, input_time_list, values_json)
                        if (self.iscloseavailable == True):
                            isclosethis = math.isclose(value, value_json, rel_tol=reltol, abs_tol=tol)
                            self.assertTrue(isclosethis, msg=(str(value)+" != "+str(value_json)+", rel_tol = "+str(reltol)+", abs_tol = "+str(tol)+
                                                              " : Error checking node "+str(node.Id)+" "+out.GetString()+" results."))
                        else:
                            self.assertAlmostEqual(value, value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" results."), delta=tol)
                    # Vector variable
                    else:
                        if (KratosMultiphysics.KratosGlobals.HasVariable( out.GetString() + "_X" )): # We will asume to be components
                            if (self.iscloseavailable == True):
                                # X-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_X"]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                isclosethis = math.isclose(value[0], value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value)+" != "+str(value_json)+", rel_tol = "+str(reltol)+", abs_tol = "+str(tol)+
                                                              " : Error checking node "+str(node.Id)+" "+out.GetString()+" results."))
                                # Y-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Y"]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                isclosethis = math.isclose(value[1], value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value)+" != "+str(value_json)+", rel_tol = "+str(reltol)+", abs_tol = "+str(tol)+
                                                              " : Error checking node "+str(node.Id)+" "+out.GetString()+" results."))
                                # Z-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Z"]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                isclosethis = math.isclose(value[2], value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value)+" != "+str(value_json)+", rel_tol = "+str(reltol)+", abs_tol = "+str(tol)+
                                                              " : Error checking node "+str(node.Id)+" "+out.GetString()+" results."))
                            else:
                                # X-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_X"]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                self.assertAlmostEqual(value[0], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" X-component results."), delta=tol)
                                # Y-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Y"]
                                value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                self.assertAlmostEqual(value[1], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" Y-component results."), delta=tol)
                                # Z-component
                                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Z"]
                                value_json = _self._linear_interpolation(time, input_time_list, values_json)
                                self.assertAlmostEqual(value[2], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" Z-component results."), delta=tol)
                        else:
                            values_json = self.data["NODE_"+str(node.Id)][out.GetString() ]
                            for index in range(len(value)):
                                value_json = self.__linear_interpolation(time, input_time_list, values_json[index])
                                if (self.iscloseavailable == True):
                                    isclosethis = math.isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value)+" != "+str(value_json)+", rel_tol = "+str(reltol)+", abs_tol = "+str(tol)+
                                                              " : Error checking node "+str(node.Id)+" "+out.GetString()+" results."))
                                else:
                                    self.assertAlmostEqual(value[index], value_json, msg=("Error checking node "+str(node.Id)+" "+out.GetString()+" results."), delta=tol)

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def __linear_interpolation(self, x, x_list, y_list):
        tb = KratosMultiphysics.PiecewiseLinearTable()
        for i in range(len(x_list)):
            tb.AddRow(x_list[i], y_list[i])

        return tb.GetNearestValue(x)

    def __generate_variable_list_from_input(self,param):
      '''Parse a list of variables from input.'''
      # At least verify that the input is a string
      if not param.IsArray():
          raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

      # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
      return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]