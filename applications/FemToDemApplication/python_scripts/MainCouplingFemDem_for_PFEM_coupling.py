

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

import KratosMultiphysics.FemToDemApplication.MainFEM_for_PFEM_coupling as FEM
import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDem_for_PFEM_coupling_Solution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================
    def __init__(self, Model, path = ""):
        # Initialize solutions

        if path == "":
            DEMProjectParametersFile = open("ProjectParametersDEM.json", 'r')
        else:
            DEMProjectParametersFile = open(os.path.join(path, "ProjectParametersDEM.json"), 'r')
        DEM_project_parameters = KratosMultiphysics.Parameters(DEMProjectParametersFile.read())

        self.FEM_Solution = FEM.FEM_for_PFEM_coupling_Solution(Model, path)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, DEM_project_parameters)

        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False
