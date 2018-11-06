from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import co_simulation_tools as tools
try :
    import KratosMultiphysics
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS is not available ! Please ensure that Kratos is available for usage !'+ tools.bcolors.ENDC)
    exit()

try :
    import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
except ModuleNotFoundError:
    print(tools.bcolors.FAIL + 'KRATOS_STRUCTURAL_MECHANICS is not available ! Please ensure that it is available for usage !'+ tools.bcolors.ENDC)
    exit()
# Importing the base class
from . import kratos_base_field_solver

# Other imports
from structural_mechanics_analysis import StructuralMechanicsAnalysis

def Create(solver_name, cosim_solver_settings):
    return KratosStructuralSolver(solver_name, cosim_solver_settings)

class KratosStructuralSolver(kratos_base_field_solver.KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return StructuralMechanicsAnalysis(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__