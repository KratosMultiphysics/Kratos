from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
import json
import numpy as np
from json_utilities import *
CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return FromJsonCheckResultProcess(Model, settings["Parameters"])

class FromJsonCheckResultProcess(Process, KratosUnittest.TestCase):
  
    def __init__(self,model_part,params):

        self.model_part = model_part

        self.params = params

        self.check_variables = []
        self.frequency    = 0.0
        self.time_counter = 0.0
        self.data = {}
        
    def ExecuteInitialize(self):
        input_file_name = self.params["input_file_name"].GetString()
        self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
        self.check_variables = self.__generate_variable_list_from_input(self.params["check_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.data =  read_external_json(input_file_name)
        
    def ExecuteBeforeSolutionLoop(self):
        pass
    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        time = self.sub_model_part.ProcessInfo.GetValue(TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
          self.time_counter = 0.0
          input_time_list = self.data["TIME"]
          for node in self.sub_model_part.Nodes:
            for i in range(self.params["check_variables"].size()):
              out = self.params["check_variables"][i]
              variable = KratosGlobals.GetVariable( out.GetString() )
              val = node.GetSolutionStepValue(variable, 0)
              if isinstance(val,float):
                value_scalar = True
              else:
                value_scalar = False

              value = node.GetSolutionStepValue(variable, 0)
              if value_scalar:
                values_json = self.data["NODE_"+str(node.Id)][out.GetString() ]
                value_json = np.interp(time, input_time_list, values_json)
                self.assertAlmostEqual(value, value_json)
              else: # It is a vector
                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_X"]
                value_json = np.interp(time, input_time_list, values_json)
                self.assertAlmostEqual(value[0], value_json)
                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Y"]
                value_json = np.interp(time, input_time_list, values_json)
                self.assertAlmostEqual(value[1], value_json)
                values_json = self.data["NODE_"+str(node.Id)][out.GetString()  + "_Z"]
                value_json = np.interp(time, input_time_list, values_json)
                self.assertAlmostEqual(value[2], value_json)
              
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
      return [ KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]

