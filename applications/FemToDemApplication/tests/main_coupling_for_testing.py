
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
    def CheckControlValuesForTesting(self): # KratosPrintInfo(str(dy))

        # Here we check the damage obtained at each FE
        for elem in self.FEM_Solution.main_model_part.Elements:
            damage = elem.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]
            if self.FEM_Solution.step == 26:
                if damage != 0.11526580049026725:
                    raise ValueError('The computed damage at step = 26 is not correct')
            elif self.FEM_Solution.step == 36:
                if damage != 0.4722648310044538:
                    raise ValueError('The computed damage at step = 36 is not correct')
            elif self.FEM_Solution.step == 46:
                if damage != 0.5600214207342531:
                    raise ValueError('The computed damage at step = 46 is not correct')
            elif self.FEM_Solution.step == 61:
                if damage != 0.5600214207342531:
                    raise ValueError('The computed damage at step = 61 is not correct')


        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(1)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 26:
            if dy != 1.971665114439254e-05:
                raise ValueError('The computed displacement at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if dy != 2.7299978508712653e-05:
                raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if dy != 2.578331303724928e-05:
                raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if dy != 1.4408321991404051e-05:
                raise ValueError('The computed displacement at step = 61 is not correct')

        # Here we check the stresses and strains at one FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        Sx = element.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED,           self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
        Ex = element.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]

        if self.FEM_Solution.step == 26:
            if Sx != 1441812.5386046136:
                raise ValueError('The computed stress at step = 26 is not correct')
            if Ex != 4.6561604584527234e-05:
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if Sx != 1190806.4343404802:
                raise ValueError('The computed stress at step = 36 is not correct')
            if Ex != 6.446991404011464e-05:
                raise ValueError('The computed strain at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if Sx != 937633.4336071612:
                raise ValueError('The computed stress at step = 46 is not correct')
            if Ex != 6.0888252148997134e-05:
                raise ValueError('The computed strain at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if Sx != 523971.6246628269:
                raise ValueError('The computed stress at step = 61 is not correct')
            if Ex != 3.4025787965616143e-05:
                raise ValueError('The computed strain at step = 61 is not correct')
        

#============================================================================================================================
    def Finalize(self):
        super(MainCouplingFemDemForTestingSolution, self).Finalize()

        shutil.rmtree(self.FEM_Solution.problem_name + "_Graphs")
        shutil.rmtree(self.FEM_Solution.problem_name + "_MPI_results")
        shutil.rmtree(self.FEM_Solution.problem_name + "_Post_Files")
        shutil.rmtree(self.FEM_Solution.problem_name + "_Results_and_Data")
        # shutil.rmtree("__pycache__")
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