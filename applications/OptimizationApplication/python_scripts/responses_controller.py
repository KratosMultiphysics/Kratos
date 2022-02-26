# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# additional imports
# Kratos Core and Apps
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
# Importing the analysis base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Additional imports
import KratosMultiphysics.OptimizationApplication.responses.response_function_factory as response_function_factory

import time as timer

# ==============================================================================
def CreateController(reponses_settings,model,model_parts_controller,analyses_controller):
    return ResponsesController(reponses_settings,model,model_parts_controller,analyses_controller)

# ==============================================================================
class ResponsesController:
    # --------------------------------------------------------------------------
    def __init__(self,reponses_settings,model,model_parts_controller,analyses_controller):
        
        self.reponses_settings = reponses_settings
        self.analyses_controller = analyses_controller
        self.model_parts_controller = model_parts_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "name"                : "REPONSE_NAME",
            "type"                : "RESPONSE_TYPE",
            "settings"                : {
                "analysis_name"   : "",
                "evaluate_model_parts": [],
                "design_types": [],
                "design_model_parts": [],
                "gradient_settings" : {}
            }
        }""")

        for itr in range(self.reponses_settings.size()):
            for key in default_settings.keys():
                if not self.reponses_settings[itr].Has(key):
                    raise RuntimeError("ResponsesController: Required entry '{}' missing in 'response Nr.{}'!".format(key,itr+1))  
            self.reponses_settings[itr].ValidateAndAssignDefaults(default_settings)            
            self.reponses_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])  


        self.responses = {}
        for itr in range(self.reponses_settings.size()):
            response_settings = self.reponses_settings[itr]
            response_name = response_settings["name"].GetString()            
            response_type = response_settings["type"].GetString()

            # check for name
            if  response_name in self.responses.keys():  
                raise RuntimeError("ResponsesController: Response name {} already exists.".format(response_name))

            response_analysis_name = response_settings["settings"]["analysis_name"].GetString()
            response_analysis = None
            # check for analysis
            if response_analysis_name !="":
                response_analysis = analyses_controller.GetAnalysis(response_analysis_name)

            # check that  design_types & design_model_parts havce the same size
            evaluate_model_parts_name_list = response_settings["settings"]["evaluate_model_parts"].GetStringArray()
            design_model_parts_name_list = response_settings["settings"]["design_model_parts"].GetStringArray()

            # check for evaluate_model_parts     
            if len(evaluate_model_parts_name_list)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(evaluate_model_parts_name_list,True)
            # check for design_model_parts                         
            if len(design_model_parts_name_list)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(design_model_parts_name_list,True)

            response = response_function_factory.CreateResponseFunction(response_name,response_type,response_settings["settings"],response_analysis,model)
            self.responses[response_name] = response


    # --------------------------------------------------------------------------
    def Initialize(self):
        for response in self.responses.values():
            response.Initialize()

    # --------------------------------------------------------------------------
    def CheckIfResponseExists(self,response_name,raise_error=True):
        if not response_name in self.responses.keys():
            if raise_error:
                raise RuntimeError("ResponsesController:CheckIfResponseExists: Response name {} does not exist.".format(response_name))
            else: 
                return False
        else:
            return True

    # --------------------------------------------------------------------------
    def CheckIfResponsesExist(self,responses_name,raise_error=True):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CheckIfResponsesExist requires list of response names")
        
        if_exist = True
        for response_name in responses_name:
            if not response_name in self.responses.keys():
                if raise_error:
                    raise RuntimeError("ResponsesController:CheckIfResponsesExist: Response {} does not exist!".format(response_name))
                else:
                    if_exist = False
                    break

        return if_exist                    

    # --------------------------------------------------------------------------
    def CalculateResponseValue(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseValue: Try to calculate response {} which does not exist.".format(response_name))
        else:
            response_value = self.responses[response_name].CalculateValue()
        return response_value

    # --------------------------------------------------------------------------
    def CalculateResponsesValues(self):
        responses_values = []
        Logger.PrintInfo("ResponsesController:CalculateResponsesValues: Starting calculation of values of responses ")
        startTime = timer.time()          
        for name,response in self.responses.items():          
            response_value = response.CalculateValue()
            responses_values.append(response_value)
        Logger.PrintInfo("ResponsesController:CalculateResponsesValues: Time needed for all response values calculation ",round(timer.time() - startTime,2),"s")            
        return responses_values
    # --------------------------------------------------------------------------
    def CalculateResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseGradient: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            self.responses[response_name].CalculateGradients()

    # --------------------------------------------------------------------------
    def CalculateResponseGradient(self,response_name,response_var_type_var_name_dict):

        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseGradient: Try to calculate {} gradient of response {} which does not exist ".format(response_var_type_var_name_dict,response_name))

        if type(response_var_type_var_name_dict) is not dict or not bool(response_var_type_var_name_dict):
            raise RuntimeError("ResponsesController:CalculateResponseGradient: the inputs should be response name and a dict of a pair of design variable type and name ")
        
        self.responses[response_name].CalculateGradient(response_var_type_var_name_dict)

    # --------------------------------------------------------------------------
    def GetResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseGradients: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            return self.responses[response_name].GetGradients()

    # --------------------------------------------------------------------------
    def GetResponseGradient(self,response_name,response_var_type_var_name_dict):

        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseGradient: Try to calculate {} gradient of response {} which does not exist ".format(response_var_type_var_name_dict,response_name))

        if type(response_var_type_var_name_dict) is not dict or not bool(response_var_type_var_name_dict):
            raise RuntimeError("ResponsesController:GetResponseGradient: the inputs should be response name and a dict of a pair of design variable type and name ")
        
        return self.responses[response_name].GetGradient(response_var_type_var_name_dict) 




            

               

