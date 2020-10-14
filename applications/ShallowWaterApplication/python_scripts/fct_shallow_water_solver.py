# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.stabilized_shallow_water_solver import StabilizedShallowWaterSolver

def CreateSolver(model, custom_settings):
    return FCTShallowWaterSolver(model, custom_settings)

class FCTShallowWaterSolver(StabilizedShallowWaterSolver):

    def Initialize(self):
        super().Initialize()
        self.main_model_part.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 2)
        fct_parameters = KM.Parameters("""{
            "rebuild_level"     : 0,
            "limiting_variable" : "HEIGHT"
        }""")
        self.fct_utility = SW.AlgebraicFluxCorrectionUtility(self.main_model_part, fct_parameters)

    def InitializeSolutionStep(self):
        print("yyyyyyyyyyyyyyyyyyyyyyyyyy")
        self.fct_utility.ExecuteInitializeLowOrderStep()
        print("zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz")
        super().InitializeSolutionStep()
        print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            print("aaaaaaaaaaaaaaaaaa")
            is_converged = self.solver.SolveSolutionStep()
            print("bbbbbbbbbbbbbbbbbbbb")
            self.fct_utility.ExecuteFinalizeLowOrderStep()
            print("ccccccccccccccccc")
            self.fct_utility.ExecuteInitializeHighOrderStep()
            print("ddddddddddddddddddddddddddddd")
            is_converged = self.solver.SolveSolutionStep()
            print("eeeeeeeeeeeeeeeeeee")
            return is_converged
        else:
            return True

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        print("hhhhhhhhhhhhhhhhhhhhhhhh")
        self.fct_utility.ExecuteFinalizeHighOrderStep()
        print("2hhhhhhhhhhhhhhhhhhhhhhh")
