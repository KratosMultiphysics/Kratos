from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import KratosMultiphysics.KratosUnittest as KratosUnittest

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================


#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.DEM_Solution.FinalizeSolutionStep()
        self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)

        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # DEM GiD print output
        # self.PrintDEMResults()

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()

        # self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Print required info
        # self.PrintPlotsFiles()
        
        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # write output results GiD: (frequency writing is controlled internally)
        # self.FEM_Solution.GraphicalOutputPrintOutput()

        # processes to be executed after writting the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        # if self.DoRemeshing:
        #      self.RemeshingProcessMMG.ExecuteFinalizeSolutionStep()
        
        # if not self.is_slave:
        #     self.PrintResults()

        # Check variables
        print("Checkeando jejejejje")


class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def small_strain(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/small_strain/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()