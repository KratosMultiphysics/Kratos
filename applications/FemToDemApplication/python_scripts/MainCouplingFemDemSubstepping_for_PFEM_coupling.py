

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem

import KratosMultiphysics.FemToDemApplication.MainFEM_for_PFEM_coupling as FEM
import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDemSubstepping as MainCouplingFemDemSubstepping
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os
import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.MeshingApplication.mmg_process as MMG

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDemSubstepping_for_PFEM_coupling_Solution(MainCouplingFemDemSubstepping.MainCoupledFemDemSubstepping_Solution):
#============================================================================================================================
    def __init__(self, Model):
        # Initialize solutions
        self.FEM_Solution = FEM.FEM_for_PFEM_coupling_Solution(Model)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model)

        # Initialize Remeshing files
        self.DoRemeshing = self.FEM_Solution.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
        if self.DoRemeshing:
            self.mmg_parameter_file = open("MMGParameters.json",'r')
            self.mmg_parameters = KratosMultiphysics.Parameters(self.mmg_parameter_file.read())
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)
        self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False
        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
