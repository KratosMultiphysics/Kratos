from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem    as MainCouplingFemDem
import KratosMultiphysics.FemToDemApplication.MainPFEM_for_coupling as MainPFEM_for_coupling

def Wait():
    input("PFEM-FEMDEM -> Press Something")

def KratosPrintInfo(message):
    """This function prints info on screen
    """
    KM.Logger.Print(message, label="")
    KM.Logger.Flush()
#============================================================================================================================
class MainCouplingPfemFemDem_Solution:
#============================================================================================================================

    def __init__(self, Model, PFEMparameters):
        # Initialize solutions of the FEMDEM and PFEM
        self.FEMDEM_Solution = MainCouplingFemDem.MainCoupledFemDem_Solution(Model)
        self.FEMDEM_Solution.Initialize()

        self.PFEM_Solution = MainPFEM_for_coupling.MainPFEM_for_coupling_solution(Model, 
                                                                                  self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                                  PFEMparameters)
        KratosPrintInfo("   ___                  _          _            _ _   _         ___  ___  __        ")
        KratosPrintInfo("  / __\___  _   _ _ __ | | ___  __| | __      _(_) |_| |__     / _ \/ __\/__\/\/\   ")
        KratosPrintInfo(" / /  / _ \| | | | '_ \| |/ _ \/ _` | \ \ /\ / / | __| '_ \   / /_)/ _\ /_\ /    \  ")
        KratosPrintInfo("/ /__| (_) | |_| | |_) | |  __/ (_| |  \ V  V /| | |_| | | | / ___/ /  //__/ /\/\ \ ")
        KratosPrintInfo("\____/\___/ \__,_| .__/|_|\___|\__,_|   \_/\_/ |_|\__|_| |_| \/   \/   \__/\/    \/ ")
        KratosPrintInfo("")

#============================================================================================================================
    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):
        # Now we storage the FEM modelpart in the PFEM_Solution
        self.PFEM_Solution.Initialize()

        # We copy the output params in the PFEM
        self.PFEM_Solution.graphical_output = self.PFEM_Solution.SetCustomGraphicalOutput(self.FEMDEM_Solution.FEM_Solution.ProjectParameters)
        self.PFEM_Solution.GraphicalOutputExecuteInitialize()

#============================================================================================================================
    def RunMainTemporalLoop(self):
        self.RunBeforeSolutionLoopFEMDEM()
        while self.FEMDEM_Solution.FEM_Solution.time <= self.FEMDEM_Solution.FEM_Solution.end_time:
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def SolveSolutionStep(self):
        # It's necessary to Fix in order to maintain the FEMDEM Kinematics
        self.FixNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

        KratosPrintInfo("==============================================")
        KratosPrintInfo("==== SOLVING PFEM PART OF THE CALCULATION ====")
        KratosPrintInfo("==============================================")
        self.SolveSolutionStepPFEM()

        # Now we Free the nodes to be calculated by the FEMDEM
        self.FreeNodesModelPart(self.FEMDEM_Solution.FEM_Solution.main_model_part)

        # Transfer pressure forces
        self.RegenerateAndUpdatePFEMPressureConditions()

        KratosPrintInfo("=================================================")
        KratosPrintInfo("==== SOLVING FEM-DEM PART OF THE CALCULATION ====")
        KratosPrintInfo("=================================================")
        self.SolveSolutionStepFEMDEM()

#============================================================================================================================
    def SolveSolutionStepPFEM(self):
        is_converged = self.PFEM_Solution._GetSolver().SolveSolutionStep()
        # self.PFEM_Solution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

#============================================================================================================================
    def SolveSolutionStepFEMDEM(self):
        self.FEMDEM_Solution.SolveSolutionStep()

#============================================================================================================================
    def RunBeforeSolutionLoopFEMDEM(self):
        # Solving the problem (time integration)
        self.FEMDEM_Solution.DEM_Solution.step           = 0
        self.FEMDEM_Solution.DEM_Solution.time           = 0.0
        self.FEMDEM_Solution.DEM_Solution.time_old_print = 0.0

        if self.FEMDEM_Solution.DoRemeshing:
            self.FEMDEM_Solution.RemeshingProcessMMG.ExecuteBeforeSolutionLoop()

#============================================================================================================================
    def FinalizeSolutionStep(self):
        self.PFEM_Solution.FinalizeSolutionStep()
        self.PFEM_Solution.OutputSolutionStep()
        self.FEMDEM_Solution.FinalizeSolutionStep()

#============================================================================================================================
    def Finalize(self):
        self.FEMDEM_Solution.Finalize()
        self.PFEM_Solution.Finalize()

#============================================================================================================================
    def InitializeSolutionStep(self):
        self.UpdateDeltaTimeInSolutions()
        self.PFEM_Solution.time = self.PFEM_Solution._GetSolver().AdvanceInTime(self.PFEM_Solution.time)
        self.PFEM_Solution.InitializeSolutionStep()
        self.PFEM_Solution._GetSolver().Predict()
        self.FEMDEM_Solution.InitializeSolutionStep()

        # We set the same delta time to all the solutions
        self.UpdateDeltaTimeInSolutions()

#============================================================================================================================
    def UpdateDeltaTimeInSolutions(self):
        FEM_delta_time = self.FEMDEM_Solution.FEM_Solution.delta_time
        self.PFEM_Solution.delta_time = FEM_delta_time
        self.PFEM_Solution.main_model_part.ProcessInfo[KM.DELTA_TIME] = FEM_delta_time

#============================================================================================================================
    def RegenerateAndUpdatePFEMPressureConditions(self):
        regenerate_cond_process = FEMDEM.RegeneratePfemPressureConditionsProcess3D(self.FEMDEM_Solution.FEM_Solution.main_model_part)
        regenerate_cond_process.Execute()
        update_cond_process = FEMDEM.UpdatePressureValuePfemConditionsProcess3D(self.FEMDEM_Solution.FEM_Solution.main_model_part)
        update_cond_process.Execute()


#============================================================================================================================
    def FixNodesModelPart(self, ModelPart):
        for node in ModelPart.Nodes:
            node.Fix(KM.VELOCITY_X)
            node.Fix(KM.VELOCITY_Y)
            node.Fix(KM.VELOCITY_Z)

    def FreeNodesModelPart(self, ModelPart):
        for node in ModelPart.Nodes:
            node.Free(KM.VELOCITY_X)
            node.Free(KM.VELOCITY_Y)
            node.Free(KM.VELOCITY_Z)

#============================================================================================================================