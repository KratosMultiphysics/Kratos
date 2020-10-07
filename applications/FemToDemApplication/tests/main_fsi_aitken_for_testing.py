
import KratosMultiphysics
from KratosMultiphysics import Logger
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
        KratosMultiphysics.PfemFluidDynamicsApplication.PostProcessUtilities().RebuildPostProcessModelPart(
                                                                                    self.PFEM_Solution.post_process_model_part,
                                                                                    self.PFEM_Solution.main_model_part)
        # self.PrintResults()

        self.CheckControlValuesForTesting()

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # message = "The obtained dx and vx is: " + str(dx) + " m  and " + str(vx) + " m/s." KratosPrintInfo(message)

        tol = 1e-5
        for node in self.FEMDEM_Solution.FEM_Solution.main_model_part.GetSubModelPart("testing_nodes").Nodes:
            # KratosPrintInfo("hey")
            if self.FEMDEM_Solution.FEM_Solution.step == 3:
                dx = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                ref_dx = 7.539979970827383e-06
                ref_vx = 0.01089025506774863
                if (dx - ref_dx) / ref_dx > tol or (vx - ref_vx) / ref_vx > tol:
                    raise ValueError('The computed displ or velocity at step = 3 is not correct')
            if self.FEMDEM_Solution.FEM_Solution.step == 5:
                dx = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                ref_dx = 7.605594669875214e-05
                ref_vx = 0.05646212212293625
                if (dx - ref_dx) / ref_dx > tol or (vx - ref_vx) / ref_vx > tol:
                    raise ValueError('The computed displ or velocity at step = 5 is not correct')
            if self.FEMDEM_Solution.FEM_Solution.step == 8:
                dx = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                ref_dx = 0.0001966109980103405
                ref_vx = 0.030178525378407893
                if (dx - ref_dx) / ref_dx > tol or (vx - ref_vx) / ref_vx > tol:
                    raise ValueError('The computed displ or velocity at step = 8 is not correct')
            if self.FEMDEM_Solution.FEM_Solution.step == 11:
                dx = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                vx = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                ref_dx = 0.0002784147973531396
                ref_vx = 0.047863965496993656
                if (dx - ref_dx) / ref_dx > tol or (vx - ref_vx) / ref_vx > tol:
                    raise ValueError('The computed displ or velocity at step = 11 is not correct')


#============================================================================================================================
    def Finalize(self):
        super(MainCouplingPfemFemDemAitkenForTestingSolution, self).Finalize()

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