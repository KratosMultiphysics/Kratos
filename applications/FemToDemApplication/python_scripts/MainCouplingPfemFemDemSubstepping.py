from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDemSubstepping_for_PFEM_coupling as MainCouplingFemDemSubstepping_for_PFEM_coupling
import KratosMultiphysics.FemToDemApplication.MainPFEM_for_coupling as MainPFEM_for_coupling
import KratosMultiphysics.FemToDemApplication.MainCouplingPfemFemDem as MainCouplingPfemFemDem

def Wait():
    input("PFEM-FEMDEM -> Press Something")

def KratosPrintInfo(message):
    """This function prints info on screen
    """
    KM.Logger.Print("", message)
    KM.Logger.Flush()
    
#============================================================================================================================
class MainCouplingPfemFemDemSubstepping_Solution(MainCouplingPfemFemDem.MainCouplingPfemFemDem_Solution):
#============================================================================================================================

    def __init__(self, Model, PFEMparameters):
        # Initialize solutions of the FEMDEM and PFEM
        self.model = Model
        self.FEMDEM_Solution = MainCouplingFemDemSubstepping_for_PFEM_coupling.MainCoupledFemDemSubstepping_for_PFEM_coupling_Solution(Model)
        self.FEMDEM_Solution.is_slave = True
        self.FEMDEM_Solution.Initialize()

        self.PFEM_Solution = MainPFEM_for_coupling.MainPFEM_for_coupling_solution(Model, 
                                                                                  self.FEMDEM_Solution.FEM_Solution.main_model_part,
                                                                                  PFEMparameters)
        KratosPrintInfo("    ___                  _          _            _ _   _         ___  ___  __       "  + "\n" +
                       "    / __\___  _   _ _ __ | | ___  __| | __      _(_) |_| |__     / _ \/ __\/__\/\/\   " + "\n" +
                       "   / /  / _ \| | | | '_ \| |/ _ \/ _` | \ \ /\ / / | __| '_ \   / /_)/ _\ /_\ /    \  " + "\n" +
                       "  / /__| (_) | |_| | |_) | |  __/ (_| |  \ V  V /| | |_| | | | / ___/ /  //__/ /\/\ \ " + "\n" +
                       "  \____/\___/ \__,_| .__/|_|\___|\__,_|   \_/\_/ |_|\__|_| |_| \/   \/   \__/\/    \/ " + "\n")