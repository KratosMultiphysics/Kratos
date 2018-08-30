from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from structural_mechanics_analysis import StructuralMechanicsAnalysis

class StructuralMechanicsAnalysisNLSensitivity(StructuralMechanicsAnalysis):
    """
    This class is the special-script of the StructuralMechanicsApplication put in a class

    It is used to specifiy the non-linear (over- resp. under-linear ) behaviour of state results (e.g. displacements or stresses).
    """
    def __init__(self, model, project_parameters):
        #solver_settings = project_parameters["solver_settings"]

        super(StructuralMechanicsAnalysisNLSensitivity, self).__init__(model, project_parameters)

    def Initialize(self):
        ##here we initialize user-provided processes
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()

        if self.is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")


