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

# Additional imports

try:
    from KratosMultiphysics.ShapeOptimizationApplication.response_functions import response_function_factory as sho_response_factory
except ImportError:
    sho_response_factory = None
try:
    from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory as csm_response_factory
except ImportError:
    csm_response_factory = None
try:
    from KratosMultiphysics.ConvectionDiffusionApplication.response_functions import convection_diffusion_response_function_factory as convdiff_response_factory
except ImportError:
    convdiff_response_factory = None
try:
    from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response_function_factory as potential_flow_response_factory
except ImportError:
    potential_flow_response_factory = None

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
            "response_name"                : "REPONSE_NAME",
            "response_type"                : "RESPONSE_TYPE",
            "response_analyzer_name"   : "RESPONSE_ANALYZER_NAME",
            "response_model_parts": [],
            "response_settings"                : {}
        }""")


        for itr in range(self.reponses_settings.size()):
            for key in default_settings.keys():
                if not self.reponses_settings[itr].Has(key):
                    raise RuntimeError("ResponsesController: Required setting '{}' missing in 'response Nr.{}'!".format(key,itr+1))  
            self.reponses_settings[itr].ValidateAndAssignDefaults(default_settings)
  
        self.responses = {}


        sho_response_functions = [
            "plane_based_packaging",
            "mesh_based_packaging",
            "surface_normal_shape_change",
            "face_angle",
            "airfoil_angle_of_attack",
            "airfoil_chord_length",
            "airfoil_perimeter"
        ]
        csm_response_functions = ["strain_energy", "mass", "eigenfrequency", "adjoint_local_stress", "adjoint_max_stress"]
        cps_response_functions = ["adjoint_lift_potential_jump", "stochastic_adjoint_lift_potential_jump"]
        convdiff_response_functions = ["point_temperature"]


        for itr in range(self.reponses_settings.size()):
            response_settings = self.reponses_settings[itr]
            response_id = itr+1
            response_name = response_settings["response_name"].GetString()            
            response_type = response_settings["response_type"].GetString()
            response_analyzer_name = response_settings["response_analyzer_name"].GetString()            

            if response_type in csm_response_functions:
                csm_response_settings = response_settings["response_settings"]               
                csm_response_settings.AddEmptyValue("response_type").SetString(response_type)
                if csm_response_factory is None:
                    raise RuntimeError("ResponsesController: Response function {} requires StructuralMechanicsApplication.".format(response_name))
                if not self.analyzers_controller.CheckIfAnalysisExists(response_analyzer_name):
                    raise RuntimeError("ResponsesController: Response {} requires analysis {} which does not exist!".format(response_name,response_analyzer_name))
                self.responses[response_name] = csm_response_factory.CreateResponseFunction(response_id, csm_response_settings, self.model)
            else:
                raise NameError("The response function '{}' of type '{}' is not available.".format(response_id, response_type ))

    # --------------------------------------------------------------------------
    def Initialize(self):
        pass           


            

               

