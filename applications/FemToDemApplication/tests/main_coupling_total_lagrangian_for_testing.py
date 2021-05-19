
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
# import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import shutil

def Wait():
    input("Press Something")

def KratosPrintInfo(message):
    KratosMultiphysics.Logger.Print(message, label="")
    KratosMultiphysics.Logger.Flush()

def IsOk(value, ref, tol):
    if abs(value - ref) / ref > tol:
        return False
    return True
#============================================================================================================================
class MainCouplingFemDemForTestingSolution(main_coupling_for_testing.MainCouplingFemDemForTestingSolution):
#============================================================================================================================

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(dy))
        tol = 1e-6
        
        # Here we check the damage obtained at each FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        damage = element.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]
        if self.FEM_Solution.step == 26:
            if damage != 0.11529157494444009:
                raise ValueError('The computed damage at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if damage != 0.47228955086766256:
                raise ValueError('The computed damage at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if damage != 0.560045933531043:
                raise ValueError('The computed damage at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if damage != 0.560045933531043:
                raise ValueError('The computed damage at step = 61 is not correct')

        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(1)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 26:
            ref = 1.972140108200845e-05
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            ref = 2.730677936681159e-05
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            ref = 2.5780167460400023e-05
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            ref = 1.4423580956034607e-05
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 61 is not correct')

        # Here we check the stresses and strains at one FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        Sx = element.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED,           self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
        Ex = element.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]

        if self.FEM_Solution.step == 26:
            if not IsOk(Sx, 1441795.525579841, tol):
                raise ValueError('The computed stress at step = 26 is not correct')
            if not IsOk(Ex, 4.6562688575946254e-05, tol):
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if not IsOk(Sx, 1190781.654419593, tol):
                raise ValueError('The computed stress at step = 26 is not correct')
            if not IsOk(Ex, 6.447199222492372e-05, tol):
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 46:
            if not IsOk(Sx, 937612.7085723856, tol):
                raise ValueError('The computed stress at step = 26 is not correct')
            if not IsOk(Ex, 6.089010583876231e-05, tol):
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 61:
            if not IsOk(Sx, 523937.42082093155, tol):
                raise ValueError('The computed stress at step = 26 is not correct')
            if not IsOk(Ex, 3.4026366842798694e-05, tol):
                raise ValueError('The computed strain at step = 26 is not correct')
        


class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def total_lagrangian(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/total_lagrangian/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()