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
from KratosMultiphysics import Parameters, Logger
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
import time as timer

# ==============================================================================
def CreateController(analyses_settings,model,model_parts_controller):
    return AnalyzersController(analyses_settings,model,model_parts_controller)

# ==============================================================================
class AnalyzersController:
    # --------------------------------------------------------------------------
    def __init__(self,analyses_settings,model,model_parts_controller):
        
        self.analyses_settings = analyses_settings
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

        for itr in range(self.analyses_settings.size()):
            for key in default_settings.keys():
                if not self.analyses_settings[itr].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis Nr.{}'!".format(key,itr+1))  
            self.analyses_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["settings"].keys():
                if not self.analyses_settings[itr]["settings"].Has(key):
                    raise RuntimeError("AnalyzersController: Required setting '{}' missing in 'analysis_settings' of 'analysis Nr.{}' !".format(key,itr+1))             
            self.analyses_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])   

        self.analyses = {}
        self.analyses_root_model_parts = {}
        self.analyses_types_model_parts_dict = {"StructuralMechanicsAnalysis":[],"ConvectionDiffusionAnalysis":[],"PotentialFlowAnalysis":[]}
        for analysis_settings in self.analyses_settings:
            analysis_name = analysis_settings["name"].GetString()            
            analysis_type = analysis_settings["type"].GetString()
            model_part_name = analysis_settings["settings"]["model_part_name"].GetString()
            input_file_name = analysis_settings["settings"]["input_filename"].GetString() 
            with open(input_file_name) as parameters_file:
                AnalysisProjectParameters = KM.Parameters(parameters_file.read())

            # check types
            if  not analysis_type in self.analyses_types_model_parts_dict.keys():  
                raise RuntimeError("AnalyzersController: Analysis {} is not supported, \n supported types are {}.".format(analysis_type,self.analyses_types_model_parts_dict.keys()))
            # check names
            if  analysis_name in self.analyses.keys():  
                raise RuntimeError("AnalyzersController: Analysis name {} already exists.".format(analysis_name))
            # check if root model part exists
            if not self.model_parts_controller.CheckIfRootModelPartExists(model_part_name):
                raise RuntimeError("AnalyzersController: Analysis {} requires model_part {} which does not exist!".format(analysis_name,model_part_name))
            # check for duplicated analysis for the same model part
            if  model_part_name in self.analyses_types_model_parts_dict[analysis_type]:
                raise RuntimeError("AnalyzersController: Analysis {} for model part {} already exists.".format(analysis_type,model_part_name))

            analysis = None
            if  analysis_type == "StructuralMechanicsAnalysis":
                if StructuralMechanicsAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires StructuralMechanicsApplication.".format(analysis_name))                 
                analysis = StructuralMechanicsAnalysis(self.model, AnalysisProjectParameters)
                self.analyses_types_model_parts_dict[analysis_type] = model_part_name 
            elif  analysis_type == "ConvectionDiffusionAnalysis":
                if ConvectionDiffusionAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires ConvectionDiffusionAnalysis.".format(analysis_name))                 
                analysis = ConvectionDiffusionAnalysis(self.model, AnalysisProjectParameters)
                self.analyses_types_model_parts_dict[analysis_type] = model_part_name       
            elif  analysis_type == "PotentialFlowAnalysis":
                if PotentialFlowAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires PotentialFlowAnalysis.".format(analysis_name))                 
                analysis = PotentialFlowAnalysis(self.model, AnalysisProjectParameters)
                self.analyses_types_model_parts_dict[analysis_type] = model_part_name
                
            self.analyses[analysis_name] = analysis
            self.analyses_root_model_parts[analysis_name] = self.model_parts_controller.GetRootModelPart(model_part_name)
                   

    # --------------------------------------------------------------------------
    def Initialize(self):
        for anal in self.analyses.values():
            anal.Initialize()

    # --------------------------------------------------------------------------
    def CheckIfAnalysisExists(self,analysis_name):
        exists = False
        for analysis_settings in self.analyses_settings:
            if analysis_name == analysis_settings["name"].GetString():
                exists = True
                break
        return exists      

    # --------------------------------------------------------------------------
    def GetAnalysis(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            return self.analyses[analysis_name]
    # --------------------------------------------------------------------------
    def GetAnalysisModelPartName(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            for analysis_settings in self.analyses_settings:
                if analysis_settings["name"].GetString() == analysis_name:
                    model_part_name = analysis_settings["settings"]["model_part_name"].GetString() 
                    return model_part_name        

    # --------------------------------------------------------------------------
    def GetAnalysisModelPart(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            for analysis_settings in self.analyses_settings:
                if analysis_settings["name"].GetString() == analysis_name:
                    model_part_name = analysis_settings["settings"]["model_part_name"].GetString() 
                    return self.model.GetModelPart(model_part_name) 

    # --------------------------------------------------------------------------
    def RunAnalysis(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalyzersController: Try to run an analysis {} which does not exist.".format(analysis_name))
        else:
            Logger.PrintInfo("AnalyzersController", "Starting analysis ", analysis_name)
            startTime = timer.time()
            analysis = self.analyses[analysis_name]
            root_model_part = self.analyses_root_model_parts[analysis_name]
            time_before_analysis = root_model_part.ProcessInfo.GetValue(KM.TIME)
            step_before_analysis = root_model_part.ProcessInfo.GetValue(KM.STEP)
            delta_time_before_analysis = root_model_part.ProcessInfo.GetValue(KM.DELTA_TIME)

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis-1)
            root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis-1)
            root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

            analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
            analysis.InitializeSolutionStep()            
            analysis._GetSolver().Predict()
            analysis._GetSolver().SolveSolutionStep()

            analysis.FinalizeSolutionStep()
            analysis.OutputSolutionStep()            

            # Clear results or modifications on model part
            root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis)
            root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis)
            root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_analysis)
            Logger.PrintInfo("AnalyzersController", "Time needed for solving the analysis",round(timer.time() - startTime,2),"s")                   

    # --------------------------------------------------------------------------
    def RunAll(self):

        for name,analysis in self.analyses.items():
            Logger.PrintInfo("AnalyzersController", "Starting analysis ", name)
            startTime = timer.time()            
            root_model_part = self.analyses_root_model_parts[name]
            time_before_analysis = root_model_part.ProcessInfo.GetValue(KM.TIME)
            step_before_analysis = root_model_part.ProcessInfo.GetValue(KM.STEP)
            delta_time_before_analysis = root_model_part.ProcessInfo.GetValue(KM.DELTA_TIME)

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis-1)
            root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis-1)
            root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

            analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
            analysis.InitializeSolutionStep()            
            analysis._GetSolver().Predict()
            analysis._GetSolver().SolveSolutionStep()

            analysis.FinalizeSolutionStep()
            analysis.OutputSolutionStep()            

            # Clear results or modifications on model part
            root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis)
            root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis)
            root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_analysis)
            Logger.PrintInfo("AnalyzersController", "Time needed for solving the analysis",round(timer.time() - startTime,2),"s")                        

               

