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
            "name"       : "ANALYSIS_NAME",
            "type"       : "ANALYSIS_TYPE",            
            "settings"              : {
                "input_filename" : "ANALYSIS_SETTINGS_FILENAME",
                "model_part_name"       : "MODEL_PART_NAME"
            }
        }""")

        for itr in range(self.analyzers_settings.size()):
            for key in default_settings.keys():
                if not self.analyzers_settings[itr].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis Nr.{}'!".format(key,itr+1))  
            self.analyzers_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["settings"].keys():
                if not self.analyzers_settings[itr]["settings"].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis_settings' of 'analysis Nr.{}' !".format(key,itr+1))             
            self.analyzers_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])   

        self.analyzers = {}
        self.analyzers_types_model_parts_dict = {"StructuralMechanicsAnalysis":[],"ConvectionDiffusionAnalysis":[],"PotentialFlowAnalysis":[]}
        for analyzer_settings in self.analyzers_settings:
            analyzer_name = analyzer_settings["name"].GetString()            
            analyzer_type = analyzer_settings["type"].GetString()
            model_part_name = analyzer_settings["settings"]["model_part_name"].GetString()
            input_file_name = analyzer_settings["settings"]["input_filename"].GetString() 
            with open(input_file_name) as parameters_file:
                AnalysisProjectParameters = KM.Parameters(parameters_file.read())

            # check types
            if  not analyzer_type in self.analyzers_types_model_parts_dict.keys():  
                raise RuntimeError("AnalyzersController: Analysis {} is not supported.".format(analyzer_type))
            # check names
            if  analyzer_name in self.analyzers.keys():  
                raise RuntimeError("AnalyzersController: Analysis name {} already exists.".format(analyzer_name))
            # check if model part exists
            if not self.model_parts_controller.CheckIfModelPartExists(model_part_name):
                raise RuntimeError("AnalyzersController: Analysis {} requires model_part {} which does not exist!".format(analyzer_name,model_part_name))
            # check for duplicated analysis for the same model part
            if  model_part_name in self.analyzers_types_model_parts_dict[analyzer_type]:
                raise RuntimeError("AnalyzersController: Analysis {} for model part {} already exists.".format(analyzer_type,model_part_name))

            if  analyzer_type == "StructuralMechanicsAnalysis":
                if StructuralMechanicsAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires StructuralMechanicsApplication.".format(analyzer_name))                 
                analyzer = StructuralMechanicsAnalysis(self.model, AnalysisProjectParameters)
                self.analyzers_types_model_parts_dict[analyzer_type] = model_part_name 

            self.analyzers[analyzer_name] = analyzer
                   

    # --------------------------------------------------------------------------
    def Initialize(self):
        for anal in self.analyzers.values():
            anal.Initialize()

    # --------------------------------------------------------------------------
    def CheckIfAnalysisExists(self,analysis_name):
        exists = False
        for analyzer_settings in self.analyzers_settings:
            if analysis_name == analyzer_settings["name"].GetString():
                exists = True
                break
        return exists            

               

