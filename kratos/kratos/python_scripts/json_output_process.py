from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
from KratosMultiphysics import *
from json_utilities import *
import json
CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return JsonOutputProcess(Model, settings["Parameters"])

class JsonOutputProcess(Process):
  
    def __init__(self,model_part,params):

        self.model_part = model_part

        self.params = params

        self.output_file_name = ""
        self.output_variables = []
        self.frequency = 0.0
        self.time_counter = 0.0
        
    def ExecuteInitialize(self):
        self.output_file_name = self.params["output_file_name"].GetString()
        self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
        self.output_variables = self.__generate_variable_list_from_input(self.params["output_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        
    def ExecuteBeforeSolutionLoop(self):
      data = {}
      data["TIME"] = []
      for node in self.sub_model_part.Nodes:
        data["NODE_"+str(node.Id)] = {}

        for i in range(self.params["output_variables"].size()):
          out = self.params["output_variables"][i]
          variable = KratosGlobals.GetVariable( out.GetString() )
          val = node.GetSolutionStepValue(variable, 0)
          if isinstance(val,float):
              data["NODE_"+str(node.Id)][out.GetString() ] = []
          else: # It is a vector
              data["NODE_"+str(node.Id)][out.GetString()  + "_X"] = []
              data["NODE_"+str(node.Id)][out.GetString()  + "_Y"] = []
              data["NODE_"+str(node.Id)][out.GetString()  + "_Z"] = []
        write_external_json(self.output_file_name, data)
            
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        data =  read_external_json(self.output_file_name)
        
        time = self.sub_model_part.ProcessInfo.GetValue(TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
          self.time_counter = 0.0
          data["TIME"].append(time)
          for node in self.sub_model_part.Nodes:
            for i in range(self.params["output_variables"].size()):
              out = self.params["output_variables"][i]
              variable = KratosGlobals.GetVariable( out.GetString() )
              val = node.GetSolutionStepValue(variable, 0)
              if isinstance(val,float):
                value_scalar = True
              else:
                value_scalar = False

              value = node.GetSolutionStepValue(variable, 0)
              if value_scalar:
                data["NODE_"+str(node.Id)][out.GetString() ].append(value)
              else: # It is a vector
                data["NODE_"+str(node.Id)][out.GetString()  + "_X"].append(value[0])
                data["NODE_"+str(node.Id)][out.GetString()  + "_Y"].append(value[1])
                data["NODE_"+str(node.Id)][out.GetString()  + "_Z"].append(value[2])
              
        write_external_json(self.output_file_name, data)

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


