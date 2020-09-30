from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import shutil

def Wait():
    input("Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(MainCouplingFemDem.MainCoupledFemDem_Solution):
#============================================================================================================================
        
"""
Main file for small strain case test
"""

#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.DEM_Solution.FinalizeSolutionStep()
        self.DEM_Solution.solver._MoveAllMeshes(self.DEM_Solution.time, self.DEM_Solution.solver.dt)

        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()
        
        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # processes to be executed after writting the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        self.CheckControlValuesForTesting()

#============================================================================================================================
    def CheckControlValuesForTesting(self):
        KratosPrintInfo("geeeeeeeeeeee")

        for elem in self.FEM_Solution.main_model_part.Elements:
            # print(elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)) # se ven numeros
            damage = elem.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]
            KratosPrintInfo(str(damage))
        # Wait()


        






#============================================================================================================================
    def Finalize(self):
        super(MainCouplingFemDemForTestingSolution, self).Finalize()

        shutil.rmtree(self.FEM_Solution.problem_name + "_Graphs")
        shutil.rmtree(self.FEM_Solution.problem_name + "_MPI_results")
        shutil.rmtree(self.FEM_Solution.problem_name + "_Post_Files")
        shutil.rmtree(self.FEM_Solution.problem_name + "_Results_and_Data")
        shutil.rmtree("__pycache__")
        os.remove("PlotFile.txt")
        os.remove(self.FEM_Solution.problem_name + "_0.post.bin")
        os.remove(self.FEM_Solution.problem_name + ".post.lst")
        os.remove("tests.post.lst")

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