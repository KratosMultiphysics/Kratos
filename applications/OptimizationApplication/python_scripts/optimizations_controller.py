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

import time as timer

# ==============================================================================
def CreateController(optimizations_settings,model,controls_controller,responses_controller):
    return OptimizationsController(optimizations_settings,model,controls_controller,responses_controller)

# ==============================================================================
class OptimizationsController:
    # --------------------------------------------------------------------------
    def __init__(self,optimizations_settings,model,controls_controller,responses_controller):
        
        self.optimizations_settings = optimizations_settings
        self.controls_controller = controls_controller
        self.responses_controller = responses_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "name"                : "OPT_NAME",
            "type"                : "OPT_TYPE",
            "settings"                : {
                "objectives": [],
                "objectives_weights": [],
                "constraints": [],
                "constraints_types": [],
                "constraints_ref_values": [],
                "controls": [],
                "controls_lower_bounds": [],
                "controls_lower_bounds_values": [],
                "controls_upper_bounds": [],
                "controls_upper_bounds_values": [],                
                "algorithm" : "ALG_NAME",
                "algorithm_settings" : {}
            }
        }""")

        if not self.optimizations_settings.size()>0:
            raise RuntimeError("OptimizationsController: optimizations list in optimizer's parameters can not be empty !") 

        for itr in range(self.optimizations_settings.size()):
            for key in default_settings.keys():
                if not self.optimizations_settings[itr].Has(key):
                    raise RuntimeError("OptimizationsController: Required entry '{}' missing in 'optimization Nr.{}'!".format(key,itr+1))  
            self.optimizations_settings[itr].ValidateAndAssignDefaults(default_settings)
            opt_name = self.optimizations_settings[itr]["name"].GetString()
            for key in default_settings["settings"].keys():
                if not self.optimizations_settings[itr]["settings"].Has(key):
                    raise RuntimeError("OptimizationsController: Required entry '{}' missing in settings of optimization '{}' !".format(key,opt_name))
            self.optimizations_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])  


        self.optimizations = {}
        self.optimizations_types={}
        self.optimizations_objectives={}
        self.optimizations_objectives_weights={}
        self.optimizations_constraints={}
        self.optimizations_constraints_type={}
        self.optimizations_constraints_ref_values={}
        self.optimizations_controls={}
        self.optimizations_controls_bounds={}
        self.optimizations_algorithm={}
        self.supported_opt_types = ["gradient_based"]
        for itr in range(self.optimizations_settings.size()):
            opt_settings = self.optimizations_settings[itr]
            opt_name = opt_settings["name"].GetString()
            opt_type = opt_settings["type"].GetString()   
            # check for name
            if opt_name in self.optimizations.keys():  
                raise RuntimeError("OptimizationsController: Optimization name '{}' already exists.".format(opt_name))       
            # check for type
            if not opt_type in self.supported_opt_types:  
                raise RuntimeError("OptimizationsController: Optimization type '{}' is not supported, supprted types {}.".format(opt_type,self.supported_opt_types))                  
            self.optimizations_types[opt_name]=opt_type
            # check for objectives
            objectives_names = opt_settings["settings"]["objectives"].GetStringArray()
            if not len(objectives_names)>0:  
                raise RuntimeError("OptimizationsController: Objectives list of optimization '{}' can not be empty.".format(opt_name))   
            if len(set(objectives_names)) != len(objectives_names):
                raise RuntimeError("OptimizationsController: Objectives list of optimization '{}' has duplicate response names .".format(opt_name))
            self.responses_controller.CheckIfResponsesExist(objectives_names)
            self.optimizations_objectives[opt_name]=objectives_names
            # check for objectives weights
            if not opt_settings["settings"]["objectives_weights"].IsVector():  
                raise RuntimeError("OptimizationsController:'objectives_weights' of optimization '{}' should be vector of doubles.".format(opt_name)) 
            objectives_weights = opt_settings["settings"]["objectives_weights"].GetVector()     
            if not len(objectives_weights) == len(objectives_names):
                raise RuntimeError("OptimizationsController:'objectives_weights' of optimization '{}' should be of the same size of objectives list.".format(opt_name))
            self.optimizations_objectives_weights[opt_name]=objectives_weights

            # check for constraints
            constraints_names = opt_settings["settings"]["constraints"].GetStringArray()
            if len(constraints_names)>0:
                if len(set(constraints_names)) != len(constraints_names):
                    raise RuntimeError("OptimizationsController: Constraint list of optimization '{}' has duplicate response names .".format(opt_name))
                self.responses_controller.CheckIfResponsesExist(constraints_names)
                for constraints_name in constraints_names:
                    if constraints_name in objectives_names:
                        raise RuntimeError("OptimizationsController: Response {} in optimization {} is used as both objective and constraint.".format(constraints_name,opt_name))
                self.optimizations_constraints[opt_name]=constraints_names

                constraints_types = opt_settings["settings"]["constraints_types"].GetStringArray()  
                if not len(constraints_types) == len(constraints_names):
                    raise RuntimeError("OptimizationsController:'constraints_types' of optimization '{}' should be of the same size of constraint list.".format(opt_name))
                for index, type in enumerate(constraints_types):
                    if not type in ["equality","inequality"]: 
                        raise RuntimeError("OptimizationsController: constraint type {} of constraint {} of optimization '{}' should be either 'equality' or 'inequality'.".format(type,constraints_names[index],opt_name,))
                self.optimizations_constraints_type[opt_name]=constraints_types                

                if not opt_settings["settings"]["constraints_ref_values"].IsVector():  
                    raise RuntimeError("OptimizationsController:'constraints_ref_values' of optimization '{}' should be vector of doubles.".format(opt_name)) 
                constraints_ref_values = opt_settings["settings"]["constraints_ref_values"].GetVector() 
                if not len(constraints_ref_values) == len(constraints_names):
                    raise RuntimeError("OptimizationsController:'constraints_ref_values' of optimization '{}' should be of the same size of constraint list.".format(opt_name))
                self.optimizations_constraints_ref_values[opt_name]=constraints_ref_values
            
            # check for controls
            controls_names = opt_settings["settings"]["controls"].GetStringArray()
            if not len(controls_names)>0:  
                raise RuntimeError("OptimizationsController: Controls list of optimization '{}' can not be empty.".format(opt_name))   
            if len(set(controls_names)) != len(controls_names):
                raise RuntimeError("OptimizationsController: Controls list of optimization '{}' has duplicate control names .".format(opt_name))
            self.controls_controller.CheckIfControlsExist(controls_names)
            self.optimizations_controls[opt_name]=controls_names

            # if not opt_settings["settings"]["controls_bounds"].IsMatrix():  
            #     raise RuntimeError("OptimizationsController:'controls_bounds' of optimization '{}' should be matrix of lower and upper bounds.".format(opt_name)) 
            # controls_bounds = opt_settings["settings"]["controls_bounds"].GetMatrix()
            # print(len(controls_bounds))
            # gjhgjg
            # if len(controls_bounds)!=len(controls_names):
            #     raise RuntimeError("OptimizationsController:'controls_bounds' of optimization '{}' should be of the same size of control list.".format(opt_name))


    # --------------------------------------------------------------------------
    def Initialize(self):
        pass




            

               

