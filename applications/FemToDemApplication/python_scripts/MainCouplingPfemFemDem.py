from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem    as MainCouplingFemDem
import KratosMultiphysics.FemToDemApplication.MainPFEM_for_coupling as MainPFEM_for_coupling

def Wait():
    input("PFEM-FEMDEM -> Press Something")

def KratosPrintInfo(message):
    """This function prints info on screen
    """
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()
#============================================================================================================================
class MainCouplingPfemFemDem_Solution:
#============================================================================================================================

    def __init__(self, Model, PFEMparameters):
        # Initialize solutions of the FEMDEM and PFEM
        self.FEMDEM_Solution = MainCouplingFemDem.MainCoupledFemDem_Solution(Model)
        self.PFEM_Solution = MainPFEM_for_coupling.MainPFEM_for_coupling_solution(Model, PFEMparameters)
        # self.ImposeIdenticalTimeStepInApplications()

#============================================================================================================================
    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def Initialize(self):
        self.FEMDEM_Solution.Initialize()
        self.PFEM_Solution.Initialize()

#============================================================================================================================
    def RunMainTemporalLoop(self):
        self.RunBeforeSolutionLoopFEMDEM()
        while self.FEMDEM_Solution.FEM_Solution.time <= self.FEMDEM_Solution.FEM_Solution.end_time:
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def SolveSolutionStep(self):
        KratosPrintInfo("=============================================")
        KratosPrintInfo("=== SOLVING PFEM PART OF THE CALCULATION ====")
        KratosPrintInfo("=============================================")
        self.SolveSolutionStepPFEM()
        # transfer pressure forces -> TODO
        # self.TransferPressureForces()
        KratosPrintInfo("================================================")
        KratosPrintInfo("=== SOLVING FEM-DEM PART OF THE CALCULATION ====")
        KratosPrintInfo("================================================")
        self.SolveSolutionStepFEMDEM()

#============================================================================================================================
    def SolveSolutionStepPFEM(self):
        is_converged = self.PFEM_Solution._GetSolver().SolveSolutionStep()
        # self.PFEM_Solution.__CheckIfSolveSolutionStepReturnsAValue(is_converged)

#============================================================================================================================
    def SolveSolutionStepFEMDEM(self):
        # self.FEMDEM_Solution.InitializeSolutionStep()
        self.FEMDEM_Solution.SolveSolutionStep()
        # self.FEMDEM_Solution.FinalizeSolutionStep()

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
        self.PFEM_Solution.time = self.PFEM_Solution._GetSolver().AdvanceInTime(self.PFEM_Solution.time)
        self.PFEM_Solution.InitializeSolutionStep()
        self.PFEM_Solution._GetSolver().Predict()
        self.FEMDEM_Solution.InitializeSolutionStep()

        # We set the same delta time to all the solutions
        # self.PFEM_Solution.delta_time = self.FEMDEM_Solution.delta_time

#============================================================================================================================
