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

# Importing the analysis base class
from KratosMultiphysics.analysis_stage import AnalysisStage

# Additional imports
import KratosMultiphysics.OptimizationApplication.responses.response_function_factory as response_function_factory

import time as time

# ==============================================================================
def CreateController(reponses_settings,model,analyzers_controller):
    return ResponsesController(reponses_settings,model,analyzers_controller)

# ==============================================================================
class ResponsesController:
    # --------------------------------------------------------------------------
    def __init__(self,reponses_settings,model,analyzers_controller):
        
        self.reponses_settings = reponses_settings
        self.analyzers_controller = analyzers_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "name"                : "REPONSE_NAME",
            "type"                : "RESPONSE_TYPE",
            "settings"                : {
                "analyzer_name"   : "RESPONSE_ANALYZER_NAME",
                "model_parts": []
            }
        }""")


        for itr in range(self.reponses_settings.size()):
            for key in default_settings.keys():
                if not self.reponses_settings[itr].Has(key):
                    raise RuntimeError("ResponsesController: Required setting '{}' missing in 'response Nr.{}'!".format(key,itr+1))  
            self.reponses_settings[itr].ValidateAndAssignDefaults(default_settings)
  
        self.responses = {}


        for itr in range(self.reponses_settings.size()):
            response_settings = self.reponses_settings[itr]
            response_name = response_settings["name"].GetString()            
            response_type = response_settings["type"].GetString()
            response_analyzer_name = ""
            if response_settings["settings"].Has("analyzer_name"):
                response_analyzer_name = response_settings["settings"]["analyzer_name"].GetString()
                response_analyzer = analyzers_controller.GetAnalysis(response_analyzer_name)
            else:
                response_analyzer = None        

            response = response_function_factory.CreateResponseFunction(response_type,response_settings["settings"],response_analyzer,model)
            self.responses[response_name] = response

            # if response_type in csm_response_functions:
            #     if csm_response_factory is None:
            #         raise RuntimeError("ResponsesController: Response function {} requires StructuralMechanicsApplication.".format(response_name))

        #         csm_response_settings = response_settings["settings"]               
        #         csm_response_settings.AddEmptyValue("type").SetString(response_type)

        #         if not self.analyzers_controller.CheckIfAnalysisExists(response_analyzer_name):
        #             raise RuntimeError("ResponsesController: Response {} requires analysis {} which does not exist!".format(response_name,response_analyzer_name))
        #         self.responses[response_name] = csm_response_factory.CreateResponseFunction(response_id, csm_response_settings, self.model)
        #     else:
        #         raise NameError("The response function '{}' of type '{}' is not available.".format(response_id, response_type ))

    # --------------------------------------------------------------------------
    def Initialize(self):
        pass           


            

               

