

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

import KratosMultiphysics.FemToDemApplication.MainFEM_for_PFEM_coupling as FEM
import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDemSubstepping as MainCouplingFemDemSubstepping
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDemSubstepping_for_PFEM_coupling_Solution(MainCouplingFemDemSubstepping.MainCoupledFemDemSubstepping_Solution):
#============================================================================================================================
    def __init__(self, Model):
        # Initialize solutions
        self.FEM_Solution = FEM.FEM_for_PFEM_coupling_Solution(Model)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model)

        self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False
        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
