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
        self.analyzers_root_model_parts = {}
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
                raise RuntimeError("AnalyzersController: Analysis {} is not supported, \n supported types are {}.".format(analyzer_type,self.analyzers_types_model_parts_dict.keys()))
            # check names
            if  analyzer_name in self.analyzers.keys():  
                raise RuntimeError("AnalyzersController: Analysis name {} already exists.".format(analyzer_name))
            # check if root model part exists
            if not self.model_parts_controller.CheckIfRootModelPartExists(model_part_name):
                raise RuntimeError("AnalyzersController: Analysis {} requires model_part {} which does not exist!".format(analyzer_name,model_part_name))
            # check for duplicated analysis for the same model part
            if  model_part_name in self.analyzers_types_model_parts_dict[analyzer_type]:
                raise RuntimeError("AnalyzersController: Analysis {} for model part {} already exists.".format(analyzer_type,model_part_name))

            analyzer = None
            if  analyzer_type == "StructuralMechanicsAnalysis":
                if StructuralMechanicsAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires StructuralMechanicsApplication.".format(analyzer_name))                 
                analyzer = StructuralMechanicsAnalysis(self.model, AnalysisProjectParameters)
                self.analyzers_types_model_parts_dict[analyzer_type] = model_part_name 
            elif  analyzer_type == "ConvectionDiffusionAnalysis":
                if ConvectionDiffusionAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires ConvectionDiffusionAnalysis.".format(analyzer_name))                 
                analyzer = ConvectionDiffusionAnalysis(self.model, AnalysisProjectParameters)
                self.analyzers_types_model_parts_dict[analyzer_type] = model_part_name       
            elif  analyzer_type == "PotentialFlowAnalysis":
                if PotentialFlowAnalysis is None:
                    raise RuntimeError("AnalyzersController: Analysis {} requires PotentialFlowAnalysis.".format(analyzer_name))                 
                analyzer = PotentialFlowAnalysis(self.model, AnalysisProjectParameters)
                self.analyzers_types_model_parts_dict[analyzer_type] = model_part_name
                
            self.analyzers[analyzer_name] = analyzer
            self.analyzers_root_model_parts[analyzer_name] = self.model_parts_controller.GetRootModelPart(model_part_name)
                   

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

    # --------------------------------------------------------------------------
    def GetAnalysis(self,analysis_name):
        if not analysis_name in self.analyzers.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            return self.analyzers[analysis_name]
    # --------------------------------------------------------------------------
    def GetAnalysisModelPartName(self,analysis_name):
        if not analysis_name in self.analyzers.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            for analyzer_settings in self.analyzers_settings:
                if analyzer_settings["name"].GetString() == analysis_name:
                    model_part_name = analyzer_settings["settings"]["model_part_name"].GetString() 
                    return model_part_name        

    # --------------------------------------------------------------------------
    def GetAnalysisModelPart(self,analysis_name):
        if not analysis_name in self.analyzers.keys():
            raise RuntimeError("AnalyzersController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            for analyzer_settings in self.analyzers_settings:
                if analyzer_settings["name"].GetString() == analysis_name:
                    model_part_name = analyzer_settings["settings"]["model_part_name"].GetString() 
                    return self.model.GetModelPart(model_part_name) 

    # --------------------------------------------------------------------------
    def RunAnalysis(self,analysis_name):
        if not analysis_name in self.analyzers.keys():
            raise RuntimeError("AnalyzersController: Try to run an analysis {} which does not exist.".format(analysis_name))
        else:
            Logger.PrintInfo("AnalyzersController", "Starting analysis ", analysis_name)
            startTime = timer.time()
            analysis = self.analyzers[analysis_name]
            root_model_part = self.analyzers_root_model_parts[analysis_name]
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

        for name,analysis in self.analyzers.items():
            Logger.PrintInfo("AnalyzersController", "Starting analysis ", name)
            startTime = timer.time()            
            root_model_part = self.analyzers_root_model_parts[name]
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

               

