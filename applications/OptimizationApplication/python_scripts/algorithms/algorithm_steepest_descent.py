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
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics import Parameters, Logger

# Additional imports
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.OptimizationApplication.utilities.timer import Timer



# ==============================================================================
class AlgorithmSteepestDescent(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self,name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):
        super().__init__(name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller)
        default_algorithm_settings = Parameters("""
        {
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3
        }""")

        self.opt_settings["algorithm_settings"].RecursivelyValidateAndAssignDefaults(default_algorithm_settings)
        self.algorithm_settings = self.opt_settings["algorithm_settings"]
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt()

        self.opt_algorithm = KOA.AlgorithmSteepestDescent(self.name,self.model,self.opt_parameters)

        Logger.PrintInfo("::[AlgorithmSteepestDescent]:: ", "Construction finished")


    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        super().InitializeOptimizationLoop()
        self.opt_algorithm.Initialize()
        self._InitializeCSVLogger()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):

        timer = Timer()
        timer.StartTimer()        

        for self.optimization_iteration in range(1,self.max_iterations+1):
            Logger.Print("")
            Logger.Print("===============================================================================")
            Logger.PrintInfo("AlgorithmSteepestDescent", "",timer.GetTimeStamp(), ": Starting optimization iteration ",self.optimization_iteration)
            Logger.Print("===============================================================================\n")

            self.model_parts_controller.UpdateTimeStep(self.optimization_iteration)
            self.opt_parameters["opt_itr"].SetInt(self.optimization_iteration)

            # update controls
            for control in self.controls:
                self.controls_controller.UpdateControl(control,False)            


            # evaluate all analysis-based responses
            for analysis in self.analyses_responses.keys():
                self.analyses_controller.RunAnalysis(analysis)
                for response in self.analyses_responses[analysis]:
                    response_value = self.responses_controller.CalculateResponseValue(response)
                    self.SetResponseValue(response,response_value)
                    self.responses_controller.CalculateResponseGradientsForTypesAndObjects(response,self.responses_control_types[response],self.responses_controlled_objects[response])

            # evaluate all analysis-free responses  
            for response in self.analysis_free_responses:
                response_value = self.responses_controller.CalculateResponseValue(response)
                self.SetResponseValue(response,response_value)
                self.responses_controller.CalculateResponseGradientsForTypesAndObjects(response,self.responses_control_types[response],self.responses_controlled_objects[response])   

            self._WriteCurrentResponseValuesToCSVFile()

            # compute objective values
            responses_value = self.responses_controller.CalculateResponsesValue(self.objectives)
            for response,value in responses_value.items():
                Logger.Print("  ========================= objective ",response,": ",value," =========================  ")

            # calculate gradients
            for response in self.responses_controlled_objects.keys():
                self.responses_controller.CalculateResponseGradientsForTypesAndObjects(response,self.responses_control_types[response],self.responses_controlled_objects[response])

            # calculate control gradients
            for control in self.controls_response_gradient_names.keys():
                for itr in range(len(self.controls_response_gradient_names[control])):
                    self.controls_controller.MapControlFirstDerivative(control,KM.KratosGlobals.GetVariable(self.controls_response_gradient_names[control][itr]),KM.KratosGlobals.GetVariable(self.controls_response_control_gradient_names[control][itr]),False)

            # calcuate 
            self.opt_algorithm.CalculateSolutionStep() 

            # compute controls
            for control in self.controls:
                self.controls_controller.ComputeControl(control,False)

            # now output 
            for control_model_part,vtkIO in self.root_model_parts_vtkIOs.items():
                root_control_model_part = self.model_parts_controller.GetRootModelPart(control_model_part)
                OriginalTime = root_control_model_part.ProcessInfo[KM.TIME]
                root_control_model_part.ProcessInfo[KM.TIME] = self.optimization_iteration

                vtkIO.ExecuteInitializeSolutionStep()
                if(vtkIO.IsOutputStep()):
                    vtkIO.PrintOutput()
                vtkIO.ExecuteFinalizeSolutionStep()

                root_control_model_part.ProcessInfo[KM.TIME] = OriginalTime
                            
            Logger.Print("")
            Logger.PrintInfo("AlgorithmSteepestDescent", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            Logger.PrintInfo("AlgorithmSteepestDescent", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

        self.FinalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        super().FinalizeOptimizationLoop()

# ==============================================================================
