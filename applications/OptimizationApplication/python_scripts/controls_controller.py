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
            "settings"       : {
                "technique"  : "CONTROL_TECHNIQUE",
                "technique_settings"       : {},
                "variables_name"  : []                
            }
        }""")


        for itr in range(self.controls_settings.size()):
            for key in default_settings.keys():
                if not self.controls_settings[itr].Has(key):
                    raise RuntimeError("ControlsController: Required setting '{}' missing in 'control Nr.{}'!".format(key,itr+1))  
            self.controls_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["settings"].keys():
                if key == "technique_settings":
                    continue
                if not self.controls_settings[itr]["settings"].Has(key):
                    raise RuntimeError("ControlsController: Required setting '{}' missing in 'settings' of 'control Nr.{}' !".format(key,itr+1))             
            self.controls_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])  


        self.supported_control_types_techniques = {"shape":["explicit_vertex_morphing"],"topology":[],"thickness":[]}
        # sanity checks
        self.controls_types_vars_dict = {"shape":[],"topology":[],"thickness":[]}
        self.controls = {}
        for itr in range(self.controls_settings.size()):
            control_settings = self.controls_settings[itr]
            control_name = control_settings["name"].GetString()            
            control_type = control_settings["type"].GetString()
            control_technique = control_settings["settings"]["technique"].GetString()
            control_variables_name_list = control_settings["settings"]["variables_name"].GetStringArray()
            # check control type
            if not control_type in self.supported_control_types_techniques.keys():
                raise RuntimeError("ControlsController: control type '{}' in control {} is not supported!, we support {} ".format(control_type,control_name,self.supported_control_types_techniques))
            # check control technique
            if not control_technique in self.supported_control_types_techniques[control_type]:
                raise RuntimeError("ControlsController: control technique '{}' for control type {} in control {} is not supported!, we support {}".format(control_technique,control_type,control_name,self.supported_control_types_techniques))

            # check for repetitious control over variables
            control_variables_name_list_set = set(control_variables_name_list)
            if len(control_variables_name_list_set) != len(control_variables_name_list):
                raise RuntimeError("ControlsController: control '{}' has duplicated control variables!".format(control_name))
            if control_variables_name_list_set.issubset(set(self.controls_types_vars_dict[control_type])):
                raise RuntimeError("ControlsController: there are duplicated {} control over {}!".format(control_type,control_variables_name_list))
            if control_name in self.controls.keys():
                raise RuntimeError("ControlsController: control name {} is duplicated !",control_name) 

            # now checks passed and create the control
            if control_type == "shape":
                # check if root model parts exist
                self.model_parts_controller.CheckIfRootModelPartsExist(control_variables_name_list,True)
                if control_technique == "explicit_vertex_morphing":
                    control = evm.ExplicitVertexMorphing(control_name,model,control_variables_name_list,control_settings["settings"]["technique_settings"])                          

            self.controls[control_name] = control
            self.controls_types_vars_dict[control_type].extend(control_variables_name_list)

                             

    # --------------------------------------------------------------------------
    def Initialize(self):
        for key,value in self.controls.items():
            value.Initialize()
           

               

