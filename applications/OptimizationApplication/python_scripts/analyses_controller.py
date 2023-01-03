# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# additional imports
import KratosMultiphysics as KM
from KratosMultiphysics import Logger
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

import time as timer

# ==============================================================================
def CreateController(analyses_settings,model,model_parts_controller):
    return AnalysesController(analyses_settings,model,model_parts_controller)

# ==============================================================================
class AnalysesController:
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
                "model_parts"       : []
            }
        }""")

        for itr in range(self.analyses_settings.size()):
            for key in default_settings.keys():
                if not self.analyses_settings[itr].Has(key):
                    raise RuntimeError("AnalysesController: Required setting '{}' missing in 'analysis Nr.{}'!".format(key,itr+1))
            self.analyses_settings[itr].ValidateAndAssignDefaults(default_settings)
            for key in default_settings["settings"].keys():
                if not self.analyses_settings[itr]["settings"].Has(key):
                    raise RuntimeError("AnalysesController: Required setting '{}' missing in 'analysis_settings' of 'analysis Nr.{}' !".format(key,itr+1))
            self.analyses_settings[itr]["settings"].ValidateAndAssignDefaults(default_settings["settings"])

        self.analyses = {}
        self.analyses_root_model_parts = {}
        self.analyses_types_model_parts_dict = {"StructuralMechanicsAnalysis":[],"ConvectionDiffusionAnalysis":[],"PotentialFlowAnalysis":[]}
        for analysis_settings in self.analyses_settings:
            analysis_name = analysis_settings["name"].GetString()
            analysis_type = analysis_settings["type"].GetString()
            model_parts_name_list = analysis_settings["settings"]["model_parts"].GetStringArray()
            input_file_name = analysis_settings["settings"]["input_filename"].GetString()
            with open(input_file_name) as parameters_file:
                AnalysisProjectParameters = KM.Parameters(parameters_file.read())

            # check types
            if  not analysis_type in self.analyses_types_model_parts_dict.keys():
                raise RuntimeError("AnalysesController: Analysis {} is not supported, \n supported types are {}.".format(analysis_type,self.analyses_types_model_parts_dict.keys()))
            # check names
            if  analysis_name in self.analyses.keys():
                raise RuntimeError("AnalysesController: Analysis name {} already exists.".format(analysis_name))
            # check if root model part exists
            self.model_parts_controller.CheckIfRootModelPartsExist(model_parts_name_list,True)

            analysis = None
            if  analysis_type == "StructuralMechanicsAnalysis":
                if StructuralMechanicsAnalysis is None:
                    raise RuntimeError("AnalysesController: Analysis {} requires StructuralMechanicsApplication.".format(analysis_name))
                analysis = StructuralMechanicsAnalysis(self.model, AnalysisProjectParameters)
            elif  analysis_type == "ConvectionDiffusionAnalysis":
                if ConvectionDiffusionAnalysis is None:
                    raise RuntimeError("AnalysesController: Analysis {} requires ConvectionDiffusionAnalysis.".format(analysis_name))
                analysis = ConvectionDiffusionAnalysis(self.model, AnalysisProjectParameters)
            elif  analysis_type == "PotentialFlowAnalysis":
                if PotentialFlowAnalysis is None:
                    raise RuntimeError("AnalysesController: Analysis {} requires PotentialFlowAnalysis.".format(analysis_name))
                analysis = PotentialFlowAnalysis(self.model, AnalysisProjectParameters)

            self.analyses[analysis_name] = analysis
            self.analyses_types_model_parts_dict[analysis_type] = model_parts_name_list
            self.analyses_root_model_parts[analysis_name] = self.model_parts_controller.GetRootModelParts(model_parts_name_list)


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
            raise RuntimeError("AnalysesController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            return self.analyses[analysis_name]
    # --------------------------------------------------------------------------
    def GetAnalysisModelPartsName(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalysesController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            for analysis_settings in self.analyses_settings:
                if analysis_settings["name"].GetString() == analysis_name:
                    return analysis_settings["settings"]["model_parts"].GetStringArray()

    # --------------------------------------------------------------------------
    def GetAnalysisModelParts(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalysesController: Try to get an analysis {} which does not exist.".format(analysis_name))
        else:
            return self.analyses_root_model_parts[analysis_name]

    # --------------------------------------------------------------------------
    def RunAnalysis(self,analysis_name):
        if not analysis_name in self.analyses.keys():
            raise RuntimeError("AnalysesController: Try to run an analysis {} which does not exist.".format(analysis_name))
        else:
            Logger.PrintInfo("AnalysesController", "Starting analysis ", analysis_name)
            startTime = timer.time()
            analysis = self.analyses[analysis_name]
            root_model_parts = self.analyses_root_model_parts[analysis_name]
            time_before_analysis = []
            step_before_analysis = []
            delta_time_before_analysis = []
            for root_model_part in root_model_parts:
                time_before_analysis.append(root_model_part.ProcessInfo.GetValue(KM.TIME))
                step_before_analysis.append(root_model_part.ProcessInfo.GetValue(KM.STEP))
                delta_time_before_analysis.append(root_model_part.ProcessInfo.GetValue(KM.DELTA_TIME))

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            for index in range(len(root_model_parts)):
                root_model_part = root_model_parts[index]
                root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis[index]-1)
                root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis[index]-1)
                root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

            analysis.time = analysis._GetSolver().AdvanceInTime(analysis.time)
            analysis.InitializeSolutionStep()
            analysis._GetSolver().Predict()
            analysis._GetSolver().SolveSolutionStep()

            analysis.FinalizeSolutionStep()
            analysis.OutputSolutionStep()

            # Clear results or modifications on model parts
            for index in range(len(root_model_parts)):
                root_model_part = root_model_parts[index]
                root_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis[index])
                root_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis[index])
                root_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_analysis[index])
            Logger.PrintInfo("AnalysesController", "Time needed for solving the analysis",round(timer.time() - startTime,2),"s")

    # --------------------------------------------------------------------------
    def RunAll(self):
        for name in self.analyses.keys():
            self.RunAnalysis(name)

    # --------------------------------------------------------------------------
    def RunAnalyses(self,analyses_name):
        if type(analyses_name) is not list:
            raise RuntimeError("AnalysesController:RunAnalyses requires list of analysis names")
        for analysis_name in analyses_name:
            self.RunAnalysis(analysis_name)



