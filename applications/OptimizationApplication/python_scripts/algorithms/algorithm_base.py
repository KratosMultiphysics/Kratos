# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ===============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics import Parameters, Logger
# ==============================================================================
class OptimizationAlgorithm:
    def __init__(self,name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):

        self.name = name
        self.opt_settings =  opt_settings
        self.model=model
        self.model_parts_controller = model_parts_controller
        self.analyses_controller = analyses_controller
        self.responses_controller = responses_controller
        self.controls_controller = controls_controller

        # objectives
        self.objectives = self.opt_settings["objectives"].GetStringArray()
        self.objectives_weights = self.opt_settings["objectives_weights"].GetVector()

        # constraints
        self.constraints = self.opt_settings["constraints"].GetStringArray()

        # all responses
        self.responses = self.objectives + self.constraints
        self.responses_type = {}
        self.responses_variables = {}
        for response in self.responses:
            self.responses_type[response] = self.responses_controller.GetResponseType(response,False)
            self.responses_variables[response] = self.responses_controller.GetResponseVariableName(response,False)



        # now extract analyses belong to the objectives
        self.analyses = self.responses_controller.GetResponsesAnalyses(self.responses)          

        # controls
        self.controls = opt_settings["controls"].GetStringArray()

        self.controls_reponses = {}
        self.controls_types = {}
        self.controls_controlling_objects = {}
        for control in self.controls:
            self.controls_types[control] = self.controls_controller.GetControlType(control)
            self.controls_controlling_objects[control] = self.controls_controller.GetControlControllingObjects(control)
            responses = self.responses_controller.GetResponses([self.controls_types[control]],self.controls_controlling_objects[control])
            if not len(responses)>0:
                raise RuntimeError("OptimizationAlgorithm: control {} is not associated to any of responses".format(control))
            if not set(responses) <= set(self.responses):
                raise RuntimeError("OptimizationAlgorithm: control {} is associated to responses {} which are not part of objectives + constraints list {} ".format(control,responses,self.responses))
            self.controls_reponses[control] = responses 

        
        # add required gradient data fields for controls 
        self.supported_control_types_variables_name = self.controls_controller.GetSupportedControlTypesVariablesName()
        for control in self.controls:
            control_type = self.controls_types[control]
            control_variable_name =  self.supported_control_types_variables_name[control_type]
            control_controlling_objects = self.controls_controlling_objects[control]
            control_controlling_root_model_parts = self.model_parts_controller.GetRootModelParts(control_controlling_objects)
            control_responses = self.controls_reponses[control]
            for response in control_responses:
                response_variable_name = self.responses_variables[response]
                response_control_gradient_field = "D_"+response_variable_name+"_D_"+control_variable_name
                for root_model in control_controlling_root_model_parts:
                    root_model.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(response_control_gradient_field))

        Logger.PrintInfo("::[OptimizationAlgorithm]:: ", "Variables ADDED")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> InitializeOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> RunOptimizationLoop << .")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        raise RuntimeError("Algorithm base class is called. Please check your implementation of the function >> FinalizeOptimizationLoop << .")

# ==============================================================================
