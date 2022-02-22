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
import KratosMultiphysics as KM 
import KratosMultiphysics.OptimizationApplication.controls.shape.explicit_vertex_morphing as evm

import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, math
import copy

# ==============================================================================
def CreateController(controls_settings,model,model_parts_controller):
    return ControlsController(controls_settings,model,model_parts_controller)

# ==============================================================================
class ControlsController:
    # --------------------------------------------------------------------------
    def __init__(self,controls_settings,model,model_parts_controller):
        
        self.controls_settings = controls_settings
        self.model_parts_controller = model_parts_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "name"       : "CONTROL_NAME",
            "type"       : "CONTROL_TYPE",
            "technique"  : "CONTROL_TECHNIQUE",
            "variables_names"  : [],
            "settings"       : {}
        }""")

        # pars the parameters
        for itr in range(self.controls_settings.size()):
            for key in default_settings.keys():
                if not self.controls_settings[itr].Has(key):
                    raise RuntimeError("ControlsController: Required setting '{}' missing in 'control Nr.{}'!".format(key,itr+1))  
            self.controls_settings[itr].ValidateAndAssignDefaults(default_settings)

        self.supported_control_types = ["shape","topology","thickness"]
        # sanity checks
        self.controls_types_vars_dict = {"shape":[],"topology":[],"thickness":[]}
        self.controls = {}
        for itr in range(self.controls_settings.size()):
            control_settings = self.controls_settings[itr]
            control_name = control_settings["name"].GetString()            
            control_type = control_settings["type"].GetString()
            control_technique = control_settings["technique"].GetString()
            control_variables_list = control_settings["variables_names"].GetStringArray()
            # check control type
            if not control_type in self.supported_control_types:
                raise RuntimeError("ControlsController: control type '{}' is not supported!".format(control_type))
            # check for repetitious control over variables
            control_variables_list_set = set(control_variables_list)
            if len(control_variables_list_set) != len(control_variables_list):
                raise RuntimeError("ControlsController: control '{}' has duplicated control variables!".format(control_name))
            if control_variables_list_set.issubset(set(self.controls_types_vars_dict[control_type])):
                raise RuntimeError("ControlsController: there are duplicated control variables")   
            if control_name in self.controls.keys():
                raise RuntimeError("ControlsController: there are duplicated control names") 

            # now checks passed and create the control
            model_parts_list = []
            for control_var in control_variables_list:
                if not self.model.HasModelPart(control_var):
                    raise RuntimeError("ControlsController: model part {} in 'control Nr.{}' does not exist in model_parts list!".format(control_var,itr+1)) 
                model_parts_list.append(self.model.GetModelPart(control_var))  
            if control_technique == "explicit_vertex_morphing":
                control = evm.ExplicitVertexMorphing(model_parts_list,control_settings["settings"])                          

            self.controls[control_name] = control
            self.controls_types_vars_dict[control_type].extend(control_variables_list)
                             

    # --------------------------------------------------------------------------
    def Initialize(self):
        for key,value in self.controls.items():
            value.Initialize()

    # --------------------------------------------------------------------------
    def CheckIfControlExists(self,control_name):
        exists = False
        for control_settings in self.controls_settings:
            if control_name == control_settings["name"].GetString():
                exists = True
                break
        return exists            

               

