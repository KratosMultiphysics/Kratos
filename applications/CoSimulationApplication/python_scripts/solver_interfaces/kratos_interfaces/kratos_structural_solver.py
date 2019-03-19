from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from . import kratos_base_field_solver

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

def CreateSolver(model, cosim_solver_settings):
    return KratosStructuralSolver(model, cosim_solver_settings)

class KratosStructuralSolver(kratos_base_field_solver.KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return StructuralMechanicsAnalysis(self.model, self.project_parameters)

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__