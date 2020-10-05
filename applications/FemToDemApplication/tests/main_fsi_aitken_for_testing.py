from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.FemToDemApplication.MainCouplingPfemFemDemAitken as MainCouplingPfemFemDemAitken
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import shutil

def Wait():
    input("Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

#============================================================================================================================
class MainCouplingPfemFemDemAitkenForTestingSolution(MainCouplingPfemFemDemAitken.MainCouplingPfemFemDemAitken_Solution):
#============================================================================================================================

#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.PFEM_Solution.FinalizeSolutionStep()
        self.PFEM_Solution.OutputSolutionStep()
        self.FEMDEM_Solution.FinalizeSolutionStep()
        KratosMultiphysics.PfemFluidDynamicsApplication.PostProcessUtilities().RebuildPostProcessModelPart(self.PFEM_Solution.post_process_model_part, self.PFEM_Solution.main_model_part)
        # self.PrintResults()

        self.CheckControlValuesForTesting()

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(dy))
        
        for node in self.FEMDEM_Solution.FEM_Solution.main_model_part.GetSubModelPart("testing_nodes").Nodes:
            KratosPrintInfo("hey")


#============================================================================================================================
    def Finalize(self):
        super(MainCouplingPfemFemDemAitkenForTestingSolution, self).Finalize()
        self.FEMDEM_Solution.FEM_Solution
        shutil.rmtree(self.FEMDEM_Solution.FEM_Solution.problem_name + "_Graphs")
        shutil.rmtree(self.FEMDEM_Solution.FEM_Solution.problem_name + "_MPI_results")
        shutil.rmtree(self.FEMDEM_Solution.FEM_Solution.problem_name + "_Post_Files")
        shutil.rmtree(self.FEMDEM_Solution.FEM_Solution.problem_name + "_Results_and_Data")
        # shutil.rmtree("__pycache__")
        os.remove("PlotFile.txt")
        os.remove(self.FEMDEM_Solution.FEM_Solution.problem_name + "_0.post.bin")
        os.remove(self.FEMDEM_Solution.FEM_Solution.problem_name + ".post.lst")
        os.remove("tests.post.lst")
        # os.remove("totalVolumeBeforeMeshing.txt")


class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def two_dimensional_fsi(self):

        with open("fsi_tests/wall_2d/PFEMProjectParameters.json",'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        
        model = KratosMultiphysics.Model()
        MainCouplingPfemFemDemAitkenForTestingSolution(model, parameters, "fsi_tests/wall_2d/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()