from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDemSubstepping_Solution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================

#============================================================================================================================
    def SolveSolutionStep(self):  # Method to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        self.ExecuteBeforeGeneratingDEM()
        self.GenerateDEM() # we create the new DEM of this time step
        self.ExecuteAfterGeneratingDEM()
        self.BeforeSolveDEMOperations()

        #### SOLVE DEM #########################################
        self.DEM_Solution.solver.Solve()
        ########################################################

        self.DEM_Solution.FinalizeSolutionStep()
        self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)

        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # DEM GiD print output
        self.PrintDEMResults()

        self.DEM_Solution.FinalizeTimeStep(self.DEM_Solution.time)

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()

        self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Update Coupled Postprocess file for Gid (post.lst)
        self.WritePostListFile()

        # Print required info
        self.PrintPlotsFiles()