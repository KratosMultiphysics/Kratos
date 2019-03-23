from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import swimming_DEM_procedures as SDP
import math
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
sys.path.insert(0, dir_path)
from swimming_DEM_analysis import SwimmingDEMAnalysis
from swimming_DEM_analysis import Say

class InterpolationTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(InterpolationTestAnalysis, self).__init__(model, varying_parameters)
        self.disperse_phase_solution.mdpas_folder_path = os.path.join(self.disperse_phase_solution.main_path, 'interpolation_tests')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def _CreateSolver(self):
        import interpolation_test_solver as sdem_solver
        return sdem_solver.InterpolationTestSolver(self.model,
                                                   self.project_parameters,
                                                   self.GetFieldUtility(),
                                                   self.fluid_solution._GetSolver(),
                                                   self.disperse_phase_solution._GetSolver(),
                                                   self.vars_man)

    def FinalizeSolutionStep(self):
        super(InterpolationTestAnalysis, self).FinalizeSolutionStep()
        node = [node for node in self.spheres_model_part.Nodes][0]

if __name__ == "__main__":
    # Setting parameters

    with open('ProjectParameters.json','r') as parameter_file:
        parameters = Parameters(parameter_file.read())

    # Create Model
    model = Model()

    # To avoid too many prints
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.DETAIL)

    test = InterpolationTestAnalysis(model, parameters)
    test.Run()
