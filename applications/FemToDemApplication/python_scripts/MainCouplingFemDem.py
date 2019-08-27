from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import FEMDEMParticleCreatorDestructor as PCD
import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import CouplingFemDem
import math
import KratosMultiphysics.MeshingApplication as MeshingApplication

def Wait():
    input("Press Something")

class MainCoupledFemDem_Solution:

#============================================================================================================================
    def __init__(self, Model):

        # Initialize solutions
        self.FEM_Solution = FEM.FEM_for_coupling_Solution(Model)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model)

        # Initialize Remeshing files
        self.DoRemeshing = self.FEM_Solution.ProjectParameters["AMR_data"]["activate_AMR"].GetBool()
        if self.DoRemeshing:
            self.mmg_parameter_file = open("MMGParameters.json",'r')
            self.mmg_parameters = KratosMultiphysics.Parameters(self.mmg_parameter_file.read())
            self.RemeshingProcessMMG = MMG.MmgProcess(Model, self.mmg_parameters)

        self.InitializePlotsFiles()
        self.echo_level = 0

#============================================================================================================================
    def Run(self):

        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================