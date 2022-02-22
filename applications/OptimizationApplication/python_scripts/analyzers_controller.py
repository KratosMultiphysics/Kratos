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

try:
    from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
except ImportError:
    StructuralMechanicsAnalysis = None
try:
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis
except ImportError:
    ConvectionDiffusionAnalysis = None
try:
    from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
except ImportError:
    PotentialFlowAnalysis = None

import KratosMultiphysics.kratos_utilities as kratos_utilities
import csv, math
import copy

# ==============================================================================
def CreateController(analyzers_settings,model,model_parts_controller):
    return AnalyzersController(analyzers_settings,model,model_parts_controller)

# ==============================================================================
class AnalyzersController:
    # --------------------------------------------------------------------------
    def __init__(self,analyzers_settings,model,model_parts_controller):
        
        self.analyzers_settings = analyzers_settings
        self.model_parts_controller = model_parts_controller
        self.model = model

        default_settings = KM.Parameters("""
        {
            "analysis_name"       : "ANALYSIS_NAME",
            "analysis_type"       : "ANALYSIS_TYPE",
            "analysis_model_part_name"       : "MODEL_PART_NAME",
            "analysis_settings"              : {
                "input_filename" : "ANALYSIS_SETTINGS_FILENAME"
            }
        }""")

        for itr in range(self.analyzers_settings.size()):
            for key in default_settings.keys():
                if not self.analyzers_settings[itr].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis Nr.{}'!".format(key,itr+1))  
            self.analyzers_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["analysis_settings"].keys():
                if not self.analyzers_settings[itr]["analysis_settings"].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis_settings' of 'analysis Nr.{}' !".format(key,itr+1))             
            self.analyzers_settings[itr]["analysis_settings"].ValidateAndAssignDefaults(default_settings["analysis_settings"])   

        self.analyzers = {}

    # --------------------------------------------------------------------------
    def Initialize(self):
        for analyzer_settings in self.analyzers_settings:
            analyzer_name = analyzer_settings["analysis_name"].GetString()            
            analyzer_type = analyzer_settings["analysis_type"].GetString()
            model_part_name = analyzer_settings["analysis_model_part_name"].GetString()
            input_file_name = analyzer_settings["analysis_settings"]["input_filename"].GetString() 
            with open(input_file_name) as parameters_file:
                AnalysisProjectParameters = KM.Parameters(parameters_file.read())         
            if  analyzer_type == "StructuralMechanicsAnalysis":
                if StructuralMechanicsAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires StructuralMechanicsApplication.".format(analyzer_name))                 
                if not self.model_parts_controller.CheckIfModelPartExists(model_part_name):
                    raise RuntimeError("AnalyzersController: Analysis {} requires model_part {} which does not exist!".format(analyzer_name,model_part_name))
                analyzer = StructuralMechanicsAnalysis(self.model, AnalysisProjectParameters) 
            self.analyzers[analyzer_name] = analyzer

    # --------------------------------------------------------------------------
    def CheckIfAnalysisExists(self,analysis_name):
        exists = False
        for analyzer_settings in self.analyzers_settings:
            if analysis_name == analyzer_settings["name"].GetString():
                exists = True
                break
        return exists            

               

