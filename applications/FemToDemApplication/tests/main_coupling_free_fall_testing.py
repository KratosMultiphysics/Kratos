
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
# import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest

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


        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(39)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 20:
            if dy != -0.0466600448812516:
                raise ValueError('The computed displacement at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            if dy != -0.3004937948812539:
                raise ValueError('The computed displacement at step = 50 is not correct')
        elif self.FEM_Solution.step == 90:
            if dy != -0.977376625965171:
                raise ValueError('The computed displacement at step = 90 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 0.5507018756981786
            if (dy - ref) / ref > 1e-4:
                raise ValueError('The computed displacement at step = 140 is not correct')
        
        vy = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
        if self.FEM_Solution.step == 20:
            if vy != -0.9564750000005906:
                raise ValueError('The computed velocity at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            if vy != -2.4279750000009126:
                raise ValueError('The computed velocity at step = 50 is not correct')
        elif self.FEM_Solution.step == 90:
            if vy != -3.2313066452534627:
                raise ValueError('The computed velocity at step = 90 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 5.457634130207384
            if abs((vy - ref) / ref) > 1e-4:
                raise ValueError('The computed displacement at step = 140 is not correct')

class TestAnalytics(KratosUnittest.TestCase):
    
    def setUp(self):
        pass

    @classmethod
    def free_fall(self):
        model = KratosMultiphysics.Model()
        MainCouplingFemDemForTestingSolution(model, "small_tests/free_fall_contact/").Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()