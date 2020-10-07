
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

#============================================================================================================================
class MainCouplingFemDemForTestingSolution(main_coupling_for_testing.MainCouplingFemDemForTestingSolution):
#============================================================================================================================

#============================================================================================================================
    def CheckControlValuesForTesting(self):  # KratosPrintInfo(str(dy))
        
        # Here we check the damage obtained at each FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        damage = element.CalculateOnIntegrationPoints(KratosFemDem.DAMAGE_ELEMENT, self.FEM_Solution.main_model_part.ProcessInfo)[0]
        if self.FEM_Solution.step == 26:
            if damage != 0.11529157494394626:
                raise ValueError('The computed damage at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if damage != 0.47228955086430224:
                raise ValueError('The computed damage at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if damage != 0.5600459335274628:
                raise ValueError('The computed damage at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if damage != 0.5600459335274628:
                raise ValueError('The computed damage at step = 61 is not correct')

        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(1)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 26:
            if dy != 1.972140108200514e-05:
                raise ValueError('The computed displacement at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if dy != 2.7306779366995683e-05:
                raise ValueError('The computed displacement at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if dy != 2.5780167460400344e-05:
                raise ValueError('The computed displacement at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if dy != 1.4423580956402616e-05:
                raise ValueError('The computed displacement at step = 61 is not correct')

        # Here we check the stresses and strains at one FE
        element = self.FEM_Solution.main_model_part.GetElement(1)
        Sx = element.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED,           self.FEM_Solution.main_model_part.ProcessInfo)[0][0]
        Ex = element.CalculateOnIntegrationPoints(KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, self.FEM_Solution.main_model_part.ProcessInfo)[0][0]

        if self.FEM_Solution.step == 26:
            if Sx != 1441795.5255806458:
                raise ValueError('The computed stress at step = 26 is not correct')
            if Ex != 4.6562688575946254e-05:
                raise ValueError('The computed strain at step = 26 is not correct')
        elif self.FEM_Solution.step == 36:
            if Sx != 1190781.6544263214:
                raise ValueError('The computed stress at step = 36 is not correct')
            if Ex != 6.447199222492372e-05:
                raise ValueError('The computed strain at step = 36 is not correct')
        elif self.FEM_Solution.step == 46:
            if Sx != 937612.7085800157:
                raise ValueError('The computed stress at step = 46 is not correct')
            if Ex != 6.089010583876231e-05:
                raise ValueError('The computed strain at step = 46 is not correct')
        elif self.FEM_Solution.step == 61:
            if Sx != 523937.4208251952:
                raise ValueError('The computed stress at step = 61 is not correct')
            if Ex != 3.4026366842798694e-05:
                raise ValueError('The computed strain at step = 61 is not correct')
        


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