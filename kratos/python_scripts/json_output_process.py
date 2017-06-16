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

class JsonOutputProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "output_variables"      : [],
            "output_file_name"      : "",
            "model_part_name"      : "",
            "sub_model_part_name"  : "",
            "time_frequency"       : 1.00,
            "resultant_solution"   : false
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        
        self.model_part = model_part

        self.params = params

        self.output_file_name = ""
        self.output_variables = []
        self.frequency = 0.0
        self.time_counter = 0.0
        self.resultant_solution = False
        
    def ExecuteInitialize(self):
        self.output_file_name = self.params["output_file_name"].GetString()
        if (len(self.params["sub_model_part_name"].GetString()) > 0):
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["sub_model_part_name"].GetString())
        else:
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
            
        self.output_variables = self.__generate_variable_list_from_input(self.params["output_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.resultant_solution = self.params["resultant_solution"].GetBool()
        
    def ExecuteBeforeSolutionLoop(self):
        data = {}
        data["TIME"] = []
        count = 0
        for node in self.sub_model_part.Nodes:
            if (self.resultant_solution == False):
                data["NODE_"+str(node.Id)] = {}
            else:
                data["RESULTANT"] = {}
                    
            for i in range(self.params["output_variables"].size()):
                out = self.params["output_variables"][i]
                variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )
                val = node.GetSolutionStepValue(variable, 0)
                if isinstance(val,float):
                    if (self.resultant_solution == False):
                        data["NODE_"+str(node.Id)][out.GetString() ] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString() ] = []
                else: # It is a vector
                    if (self.resultant_solution == False):
                        data["NODE_"+str(node.Id)][out.GetString()  + "_X"] = []
                        data["NODE_"+str(node.Id)][out.GetString()  + "_Y"] = []
                        data["NODE_"+str(node.Id)][out.GetString()  + "_Z"] = []
                    else:
                        if (count == 0):
                            data["RESULTANT"][out.GetString()  + "_X"] = []
                            data["RESULTANT"][out.GetString()  + "_Y"] = []
                            data["RESULTANT"][out.GetString()  + "_Z"] = []
            count += 1
            
        write_external_json(self.output_file_name, data)
        
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        data =  read_external_json(self.output_file_name)
        
        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            data["TIME"].append(time)
            count = 0
            for node in self.sub_model_part.Nodes:
                for i in range(self.params["output_variables"].size()):
                    out = self.params["output_variables"][i]
                    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )
                    val = node.GetSolutionStepValue(variable, 0)

                    value = node.GetSolutionStepValue(variable, 0)
                    if isinstance(val,float):
                        if (self.resultant_solution == False):
                            data["NODE_"+str(node.Id)][out.GetString() ].append(value)
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString() ].append(value)
                            else:
                                data["RESULTANT"][out.GetString() ][-1] += value
                    else: # It is a vector
                        if (self.resultant_solution == False):
                            data["NODE_"+str(node.Id)][out.GetString()  + "_X"].append(value[0])
                            data["NODE_"+str(node.Id)][out.GetString()  + "_Y"].append(value[1])
                            data["NODE_"+str(node.Id)][out.GetString()  + "_Z"].append(value[2])
                        else:
                            if (count == 0):
                                data["RESULTANT"][out.GetString()  + "_X"].append(value[0])
                                data["RESULTANT"][out.GetString()  + "_Y"].append(value[1])
                                data["RESULTANT"][out.GetString()  + "_Z"].append(value[2])
                            else:
                                data["RESULTANT"][out.GetString()  + "_X"][-1] += value[0]
                                data["RESULTANT"][out.GetString()  + "_Y"][-1] += value[1]
                                data["RESULTANT"][out.GetString()  + "_Z"][-1] += value[2]
                count += 1
              
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
      return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]


