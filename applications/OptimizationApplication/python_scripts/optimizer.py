# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# Kratos Core and Apps
import KratosMultiphysics as KM
from KratosMultiphysics.OptimizationApplication import model_parts_controller
from KratosMultiphysics.OptimizationApplication import analyses_controller
from KratosMultiphysics.OptimizationApplication import responses_controller
from KratosMultiphysics.OptimizationApplication import controls_controller
from KratosMultiphysics.OptimizationApplication import optimizations_controller

# additional imports


## Purely for backward compatibility, should be removed soon.
def CreateOptimizer(optimization_settings,model):
    return Optimizer(optimization_settings,model)

def Create(optimization_settings,model):
    return Optimizer(optimization_settings,model)

# ==============================================================================
class Optimizer:
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings,model):
        self._ValidateSettings(optimization_settings)
        self.optimization_settings = optimization_settings
        self.model_parts_controller = model_parts_controller.CreateController(optimization_settings["model_parts"],model)
        self.analyses_controller = analyses_controller.CreateController(optimization_settings["analyses"],model,self.model_parts_controller)
        self.responses_controller = responses_controller.CreateController(optimization_settings["responses"],model,self.model_parts_controller,self.analyses_controller)
        self.controls_controller = controls_controller.CreateController(optimization_settings["controls"],model,self.model_parts_controller)
        self.optimizations_controller = optimizations_controller.CreateController(optimization_settings["optimizations"],model,self.model_parts_controller,self.analyses_controller,self.responses_controller,self.controls_controller)


    def _ValidateSettings(self, optimization_settings):
        self._ValidateTopLevelSettings(optimization_settings)
    # ------------------------------------------------------------------------------
    def _ValidateTopLevelSettings(self, optimization_settings):
        default_settings = KM.Parameters("""
        {
            "model_parts" : [ ],
            "analyses" : [ ],
            "responses" : [ ],            
            "controls" : [ ],
            "optimizations" : [ ]
        }""")

        for key in default_settings.keys():
            if not optimization_settings.Has(key):
                raise RuntimeError("Optimizer: Required setting '{}' missing in optimizer's input parameters !".format(key))

        optimization_settings.ValidateAndAssignDefaults(default_settings)

       
