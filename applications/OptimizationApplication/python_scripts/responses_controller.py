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
import KratosMultiphysics.OptimizationApplication.responses.analysis_based_response_function_factory as analysis_based_response_function_factory
import KratosMultiphysics.OptimizationApplication.responses.analysis_free_response_function_factory as analysis_free_response_function_factory
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

        required_settings = KM.Parameters("""
        {
            "name"                : "REPONSE_NAME",
            "type"                : "RESPONSE_TYPE",
            "settings"                : {
                "evaluated_model_parts": [],
                "controlled_model_parts": [],
                "control_types": []                
            }
        }""")

      
        for itr in range(self.reponses_settings.size()):
            for key in required_settings.keys():
                if not self.reponses_settings[itr].Has(key):
                    raise RuntimeError("ResponsesController: Required entry '{}' missing in 'response Nr.{}'!".format(key,itr+1))     

            for key in required_settings["settings"].keys():
                if not self.reponses_settings[itr]["settings"].Has(key):
                    raise RuntimeError("ResponsesController: Required entry '{}' missing in settings of 'response Nr.{}'!".format(key,itr+1)) 


        analysis_based_responses = []
        structural_responses = ["strain_energy", "eigenfrequency", "adjoint_local_stress", "adjoint_max_stress"]        
        analysis_based_responses.extend(structural_responses)

        analysis_free_responses = []
        shape_opt_responses = ["plane_based_packaging","mesh_based_packaging","surface_normal_shape_change","face_angle","airfoil_chord_length",
                               "airfoil_perimeter","total_volume"]                                
        analysis_free_responses.extend(shape_opt_responses)

        supported_responses = []
        supported_responses.extend(analysis_based_responses) 
        supported_responses.extend(analysis_free_responses)


        self.responses = {}
        self.responses_analyses = {}
        self.responses_evaluated_model_parts = {}
        self.responses_controlled_model_parts = {}
        self.responses_control_types = {}
        for itr in range(self.reponses_settings.size()):
            response_settings = self.reponses_settings[itr]
            response_name = response_settings["name"].GetString()            
            response_type = response_settings["type"].GetString()

            # check for name
            if  response_name in self.responses.keys():  
                raise RuntimeError("ResponsesController: Response name {} already exists.".format(response_name))

            # check for evaluated_model_parts
            evaluated_model_parts = response_settings["settings"]["evaluated_model_parts"].GetStringArray()                  
            if len(evaluated_model_parts)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(evaluated_model_parts,True)
            else:
                raise RuntimeError("ResponsesController: 'evaluated_model_parts' of response {} can not be empty.".format(response_name))

            self.responses_evaluated_model_parts[response_name] = evaluated_model_parts
            
            # check for controlled_model_parts   
            controlled_model_parts = response_settings["settings"]["controlled_model_parts"].GetStringArray()                      
            if len(controlled_model_parts)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(controlled_model_parts,True)
            else:
                raise RuntimeError("ResponsesController: 'controlled_model_parts' of response {} can not be empty.".format(response_name))                
            self.responses_controlled_model_parts[response_name] = controlled_model_parts

            # check if controlled_model_parts and evaluated_model_parts have the same root model parts
            first_root_evaluated_model_part = evaluated_model_parts[0].split(".")[0]

            for model_part in evaluated_model_parts:
                if model_part.split(".")[0] != first_root_evaluated_model_part:
                    raise RuntimeError("ResponsesController: 'evaluated_model_parts' of response {} must have the same root model part.".format(response_name))                     

            first_root_controlled_model_part = controlled_model_parts[0].split(".")[0]
            for model_part in controlled_model_parts:
                if model_part.split(".")[0] != first_root_controlled_model_part:
                    raise RuntimeError("ResponsesController: 'controlled_model_parts' of response {} must have the same root model part.".format(response_name))  

            if first_root_controlled_model_part != first_root_evaluated_model_part:
                raise RuntimeError("ResponsesController: 'controlled_model_parts' and 'evaluated_model_parts' of response {} must have the same root model parts.".format(response_name))                  

            # check for control_types
            control_types = response_settings["settings"]["control_types"].GetStringArray()
            if len(control_types) != len(controlled_model_parts):
                raise RuntimeError("ResponsesController: 'control_types' of response {} must be of the same size as controlled_model_parts .".format(response_name)) 
            self.responses_control_types[response_name] = control_types              

            response = None
            if response_type in analysis_based_responses:
                if not response_settings["settings"].Has("analysis_name"):
                    raise RuntimeError("ResponsesController: Response {} of type {} requires analysis with analysis_name entry in settings".format(response_name,response_type)) 
                response_analysis_name = response_settings["settings"]["analysis_name"].GetString()
                response_analysis = analyses_controller.GetAnalysis(response_analysis_name)
                response = analysis_based_response_function_factory.CreateResponseFunction(response_name,response_type,response_settings["settings"],response_analysis,model)
                self.responses_analyses[response_name] = response_analysis_name
            elif response_type in analysis_free_responses:
                response = analysis_free_response_function_factory.CreateResponseFunction(response_name,response_type,response_settings["settings"],model)
                self.responses_analyses[response_name] = None
            else:
                raise RuntimeError("ResponsesController: Response {} of type {} is not supported, supported responses are {}".format(response_name,response_type,supported_responses)) 
            
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
    def CalculateResponsesValue(self,responses_name):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CalculateResponsesValue requires list of response names")

        responses_value = []
        for response_name in responses_name:
            if not response_name in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponsesValue: Try to calculate response {} which does not exist.".format(response_name))            
            response_value = self.responses[response_name].CalculateValue()
            responses_value.append(response_value)

        return responses_value        

    # --------------------------------------------------------------------------
    def CalculateAllValues(self):
        responses_values = []
        Logger.PrintInfo("ResponsesController:CalculateAllValues: Starting calculation of values of responses ")
        startTime = timer.time()          
        for name,response in self.responses.items():          
            response_value = response.CalculateValue()
            responses_values.append(response_value)
        Logger.PrintInfo("ResponsesController:CalculateAllValues: Time needed for all response values calculation ",round(timer.time() - startTime,2),"s")            
        return responses_values
    # --------------------------------------------------------------------------
    def CalculateResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseGradient: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            self.responses[response_name].CalculateGradients()

    # --------------------------------------------------------------------------
    def CalculateResponsesGradients(self,responses_name):

        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CalculateResponsesValue requires list of response names")

        for response_name in responses_name:
            if not response_name in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponsesGradients: Try to calculate gradients of response {} which does not exist.".format(response_name))            
            self.responses[response_name].CalculateGradients()

    # --------------------------------------------------------------------------
    def GetResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseGradients: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            return self.responses[response_name].GetGradients()

    # --------------------------------------------------------------------------
    def GetResponseAnalysis(self,response_name):

        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseAnalysis: Try to get analysis of response {} which does not exist ".format(response_name))
        
        return self.responses_analyses[response_name]
    # --------------------------------------------------------------------------
    def GetResponsesAnalyses(self,responses_name):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:GetResponsesAnalyses requires list of response names")        

        analyses_list = []
        for response_name in responses_name:
            analyses_list.append(self.responses_analyses[response_name])
        
        return analyses_list     


            

               

