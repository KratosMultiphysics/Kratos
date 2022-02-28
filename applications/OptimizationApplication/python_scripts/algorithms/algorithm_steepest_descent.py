# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
import KratosMultiphysics.OptimizationApplication as KOPT

# Additional imports
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer


# ==============================================================================
class AlgorithmSteepestDescent(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self,name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):
        default_algorithm_settings = KM.Parameters("""
        {
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3
        }""")

        self.name = name
        self.opt_settings =  opt_settings
        self.opt_settings["algorithm_settings"].RecursivelyValidateAndAssignDefaults(default_algorithm_settings)
        self.algorithm_settings = self.opt_settings["algorithm_settings"]
        self.model=model
        self.model_parts_controller = model_parts_controller
        self.analyses_controller = analyses_controller
        self.responses_controller = responses_controller
        self.controls_controller = controls_controller

        self.objectives = self.opt_settings["objectives"].GetStringArray()
        self.objectives_weights = self.opt_settings["objectives_weights"].GetVector()


    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt()
        # now extract analyses belong to the objectives
        self.analyses = self.responses_controller.GetResponsesAnalyses(self.objectives)        
    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):

        timer = Timer()
        timer.StartTimer()        

        for self.optimization_iteration in range(1,self.max_iterations+1):
            KM.Logger.Print("")
            KM.Logger.Print("===============================================================================")
            KM.Logger.PrintInfo("Opt", "",timer.GetTimeStamp(), ": Starting optimization iteration ",self.optimization_iteration)
            KM.Logger.Print("===============================================================================\n")

            self.model_parts_controller.UpdateTimeStep(self.optimization_iteration)

            self.analyses_controller.RunAnalyses(self.analyses)

            self.responses_controller.CalculateResponsesValue(self.objectives)
            self.responses_controller.CalculateResponsesGradients(self.objectives)


        #     timer.StartNewLap()

        #     self.__initializeNewShape()

        #     self.__analyzeShape()

        #     if self.line_search_type == "adaptive_stepping" and self.optimization_iteration > 1:
        #         self.__adjustStepSize()

        #     self.__computeShapeUpdate()

        #     self.__logCurrentOptimizationStep()

            KM.Logger.Print("")
            KM.Logger.PrintInfo("Opt", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            KM.Logger.PrintInfo("Opt", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

        #     if self.__isAlgorithmConverged():
        #         break
        #     else:
        #         self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        pass
        # self.data_logger.FinalizeDataLogging()
        # self.analyzer.FinalizeAfterOptimizationLoop()

# ==============================================================================
