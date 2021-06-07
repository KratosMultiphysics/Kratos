
import KratosMultiphysics
from KratosMultiphysics import Logger
import KratosMultiphysics.FemToDemApplication as KratosFemDem
# import KratosMultiphysics.FemToDemApplication.MainCouplingFemDem as MainCouplingFemDem
import main_coupling_for_testing
import KratosMultiphysics.KratosUnittest as KratosUnittest
import os

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

        tol = 1e-8
        # Here we check the vertical displacement of a node
        node = self.FEM_Solution.main_model_part.GetNode(39)
        dy = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
        if self.FEM_Solution.step == 20:
            ref = -0.0466600448812516
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            ref = -0.3004937948812539
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 50 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 0.6400295473530965
            if abs((dy-ref)/ref) > tol:
                raise ValueError('The computed displacement at step = 140 is not correct')

        vy = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
        if self.FEM_Solution.step == 20:
            ref = -0.9564750000005906
            if abs((vy-ref)/ref) > tol:
                raise ValueError('The computed velocity at step = 20 is not correct')
        elif self.FEM_Solution.step == 50:
            ref = -2.4279750000009126
            if abs((vy-ref)/ref) > tol:
                raise ValueError('The computed velocity at step = 50 is not correct')
        elif self.FEM_Solution.step == 140:
            ref = 5.61095572433878
            if abs((vy-ref)/ref) > tol:
                raise ValueError('The computed velocity at step = 140 is not correct')

class TestAnalytics(KratosUnittest.TestCase):

    def setUp(self):
        pass

    @classmethod
    def free_fall(self):
        model = KratosMultiphysics.Model()

        MainCouplingFemDemForTestingSolution(model, os.path.join(os.path.abspath(os.path.dirname(__file__)), "small_tests", "free_fall_contact")).Run()



if __name__ == "__main__":
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)
    KratosUnittest.main()