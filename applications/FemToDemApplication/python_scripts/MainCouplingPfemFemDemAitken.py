from __future__ import print_function, absolute_import, division  #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.FemToDemApplication as FEMDEM
import KratosMultiphysics.PfemFluidDynamicsApplication as PFEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem_for_PFEM_coupling as MainCouplingFemDem_for_PFEM_coupling
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
class MainCouplingPfemFemDemAitken_Solution(MainCouplingPfemFemDem.MainCouplingPfemFemDem_Solution):
#============================================================================================================================

#============================================================================================================================
    def SolveSolutionStep(self):
        super(MainCouplingPfemFemDemAitken_Solution, self).SolveSolutionStep()