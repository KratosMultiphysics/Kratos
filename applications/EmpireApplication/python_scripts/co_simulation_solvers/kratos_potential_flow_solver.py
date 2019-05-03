
from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
try:
    import KratosMultiphysics.MeshMovingApplication
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "succesfully imported")
except ImportError:
    KratosMultiphysics.Logger.PrintInfo("MeshMovingApplication", "not imported")

# Importing the base class
from kratos_base_field_solver import KratosBaseFieldSolver

# Other imports
from KratosMultiphysics.CompressiblePotentialFlowApplication.potential_flow_analysis import PotentialFlowAnalysis
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_forces_on_nodes_process import ComputeForcesOnNodesProcess
from KratosMultiphysics.CompressiblePotentialFlowApplication.define_wake_process_2d import DefineWakeProcess2D

def CreateSolver(cosim_solver_settings, level):
    return KratosPotentialFlowSolver(cosim_solver_settings, level)

class KratosPotentialFlowSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        #print('debug: project_param. from KratosPotentialFlowSolver:\n',  self.project_parameters)
        KratosMultiphysics.Logger.PrintInfo("Creation of Model for KratosPotentialFlowSolver is used")
        return PotentialFlowAnalysis(self.model, self.project_parameters)

    def InitializeSolutionStep(self):
        super(KratosPotentialFlowSolver, self).InitializeSolutionStep()

        sub_project_parameters = self.project_parameters["processes"]["boundary_conditions_process_list"]

        self.wake_process = DefineWakeProcess2D(self.model, sub_project_parameters[1]["Parameters"])

        if not hasattr(self, "wake_process"):
            raise Exception("potential flow requires specification of a process for the wake (currently specifically using 'define_wake_process_2d')")
        self.conversion_process = ComputeForcesOnNodesProcess(self.model, sub_project_parameters[4]["Parameters"])

        # KratosMultiphysics.Logger.PrintInfo("REACTIONS for this case are == ", KratosMultiphysics.REACTION)
        # self.conversion_process.Execute()

    def SolveSolutionStep(self):
        self.wake_process.FindWakeElements()
        KratosMultiphysics.Logger.PrintInfo("IIIIIIIIIII", self.wake_process)

        super(KratosPotentialFlowSolver, self).SolveSolutionStep()

        self.conversion_process.ExecuteFinalizeSolutionStep()
        self.wake_process.CleanMarking()

    def FinalizeSolutionStep(self):
        # self.wake_process.FindWakeElements()
        # self.conversion_process.ExecuteFinalizeSolutionStep()
        super(KratosPotentialFlowSolver, self).FinalizeSolutionStep()


        # self.wake_process.CleanMarking()

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__

