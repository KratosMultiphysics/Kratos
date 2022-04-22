# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ===============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM
from KratosMultiphysics.LinearSolversApplication import dense_linear_solver_factory
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics import Parameters, Logger

# Additional imports
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_base import OptimizationAlgorithm
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer

import csv, os

# ==============================================================================
class AlgorithmGradientProjection(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self,name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller):
        super().__init__(name,opt_settings,model,model_parts_controller,analyses_controller,responses_controller,controls_controller)
        default_algorithm_settings = Parameters("""
        {
            "max_iterations"     : 100,
            "relative_tolerance" : 1e-3,
            "alpha": 0.3
        }""")

        self.opt_settings["algorithm_settings"].RecursivelyValidateAndAssignDefaults(default_algorithm_settings)
        self.algorithm_settings = self.opt_settings["algorithm_settings"]
        self.max_iterations = self.algorithm_settings["max_iterations"].GetInt()
        self.opt_parameters.AddDouble("alpha",self.opt_settings["algorithm_settings"]["alpha"].GetDouble())

        # check if constraint list is empty or not 
        if not len(self.constraints)>0:
            raise RuntimeError("AlgorithmGradientProjection:__init__: constraints list can not be empty !")  

        self.lin_solver_settings = KM.Parameters('{ "solver_type" : "LinearSolversApplication.dense_col_piv_householder_qr" }')
        self.lin_solver = dense_linear_solver_factory.ConstructSolver(self.lin_solver_settings)          

        self.opt_algorithm = KOA.AlgorithmGradientProjection(self.name,self.model,self.lin_solver,self.opt_parameters)

        Logger.PrintInfo("::[AlgorithmGradientProjection]:: ", "Construction finished")


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
            Logger.PrintInfo("AlgorithmGradientProjection", "",timer.GetTimeStamp(), ": Starting optimization iteration ",self.optimization_iteration)
            Logger.Print("===============================================================================\n")

            self.model_parts_controller.UpdateTimeStep(self.optimization_iteration)
            self.opt_parameters["opt_itr"].SetInt(self.optimization_iteration)

            # update controls
            for control in self.controls:
                self.controls_controller.UpdateControl(control,False)            

            # do all analyses
            self.analyses_controller.RunAnalyses(self.analyses)

            # compute objectives values
            objectives_values = self.responses_controller.CalculateResponsesValue(self.objectives)
            self.objectives_hist.append(list(objectives_values.values()))
            for name,value in objectives_values.items():
                Logger.Print("  ========================= objective ",name,": ",value," =========================  ")
                # set value
                for objective in self.opt_parameters["objectives"]:
                    if objective["name"].GetString() == name:
                        objective["value"].SetDouble(value)
        
            # compute constraint values
            constraints_values = self.responses_controller.CalculateResponsesValue(self.constraints)
            self.constraints_hist.append(list(constraints_values.values()))
            for name,value in constraints_values.items():
                Logger.Print("  ========================= Constraint ",name,": ",value," =========================  ")
                # set value
                for constraint in self.opt_parameters["constraints"]:
                    if constraint["name"].GetString() == name:
                        constraint["value"].SetDouble(value)        
            self._WriteCurrentResponseValuesToCSVFile()

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
            Logger.PrintInfo("AlgorithmGradientProjection", "Time needed for current optimization step = ", timer.GetLapTime(), "s")
            Logger.PrintInfo("AlgorithmGradientProjection", "Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

        #     if self.__isAlgorithmConverged():
        #         break
        #     else:
        #         self.__determineAbsoluteChanges()

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        super().FinalizeOptimizationLoop()

    # --------------------------------------------------------------------------
    def _InitializeCSVLogger(self):
        self.complete_log_file_name = "optimization_log.csv"
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            for itr in range(len(self.objectives)):
                row.append("{:>4s}".format("f: "+str(self.objectives[itr])))
                row.append("{:>4s}".format("abs[%]"))
                row.append("{:>4s}".format("rel[%]"))

            for itr in range(len(self.constraints)):
                row.append("{:>4s}".format("c: "+str(self.constraints[itr])))
                row.append("{:>4s}".format("ref_val "))
                row.append("{:>4s}".format("ref_diff[%]"))

            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentResponseValuesToCSVFile( self ):
        
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.optimization_iteration))
            if self.optimization_iteration>1:                
                objectives_initial_values = self.objectives_hist[0]
                objectives_previous_values = self.objectives_hist[self.optimization_iteration-2]
                objectives_current_values = self.objectives_hist[self.optimization_iteration-1]
                for obj_i in range(len(objectives_current_values)):
                    current_val = objectives_current_values[obj_i]
                    rel_change = 100 * (current_val-objectives_previous_values[obj_i])/objectives_previous_values[obj_i]
                    abs_change = 100 * (current_val-objectives_initial_values[obj_i])/objectives_initial_values[obj_i]
                    row.append(" {:> .5E}".format(current_val))
                    row.append(" {:> .5E}".format(abs_change))
                    row.append(" {:> .5E}".format(rel_change))
                    
            else:
                objectives_initial_values = self.objectives_hist[self.optimization_iteration-1]                
                for obj_i in range(len(objectives_initial_values)):
                    initial_val = objectives_initial_values[obj_i]
                    row.append(" {:> .5E}".format(initial_val))
                    row.append(" {:> .5E}".format(0))
                    row.append(" {:> .5E}".format(0))

            constraints_current_values = self.constraints_hist[self.optimization_iteration-1]
            for cons_i in range(len(constraints_current_values)):
                current_val = constraints_current_values[cons_i]
                ref_val = self.constraints_ref_values[cons_i]
                abs_change = 100 * abs(current_val-ref_val)/abs(ref_val)
                row.append(" {:> .5E}".format(current_val))
                row.append(" {:> .5E}".format(ref_val))
                row.append(" {:> .5E}".format(abs_change))

            historyWriter.writerow(row)





# ==============================================================================
