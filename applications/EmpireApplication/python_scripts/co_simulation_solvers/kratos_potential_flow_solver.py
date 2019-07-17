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
from KratosMultiphysics.CompressiblePotentialFlowApplication.compute_lift_process import ComputeLiftProcess
from gid_output_process import Factory as GiDFactory

def CreateSolver(cosim_solver_settings, level):
    return KratosPotentialFlowSolver(cosim_solver_settings, level)

class KratosPotentialFlowSolver(KratosBaseFieldSolver):
    def _CreateAnalysisStage(self):
        return PotentialFlowAnalysis(self.model, self.project_parameters)

    def Initialize(self):

        super(KratosPotentialFlowSolver, self).Initialize()

        sub_project_parameters = self.project_parameters["processes"]["boundary_conditions_process_list"]

        for i in range(sub_project_parameters.size()):
            if sub_project_parameters[i]["python_module"].GetString() == "define_wake_process_2d":
                self.wake_process = DefineWakeProcess2D(self.model, sub_project_parameters[i]["Parameters"])
                if not hasattr(self, "wake_process"):
                    raise Exception("potential flow requires specification of a process for the wake (currently specifically using 'define_wake_process_2d')")

            if sub_project_parameters[i]["python_module"].GetString() == "compute_forces_on_nodes_process":
                self.conversion_process = ComputeForcesOnNodesProcess(self.model, sub_project_parameters[i]["Parameters"])
            if sub_project_parameters[i]["python_module"].GetString() == "compute_lift_process":
                self.lift_process = ComputeLiftProcess(self.model, sub_project_parameters[i]["Parameters"])

        gid_parameters = self.project_parameters["output_processes"]["gid_output"][0]
        gid_parameters["Parameters"]["output_name"].SetString('Fsi_Potential_Flow_Output_with_TimeSteps')
        self.gid_output_process = GiDFactory(gid_parameters,self.model)
        self.gid_output_process.ExecuteInitialize()
        self.gid_output_process.ExecuteBeforeSolutionLoop()
        self.time = 1

    def SolveSolutionStep(self):

        self.wake_process.ExecuteInitialize()
        self.gid_output_process.ExecuteInitializeSolutionStep()
        self._GetAnalysisStage()._GetSolver().fluid_solver.solver.Clear()
        self._GetAnalysisStage()._GetSolver().fluid_solver.solver.InitializeSolutionStep()
        super(KratosPotentialFlowSolver, self).SolveSolutionStep()
        self.model["FluidModelPart"].CloneTimeStep(self.time)
        self.model["FluidModelPart"].ProcessInfo[KratosMultiphysics.STEP] += 1

        self.conversion_process.ExecuteFinalizeSolutionStep()
        self.gid_output_process.ExecuteFinalizeSolutionStep()
        self.lift_process.ExecuteFinalizeSolutionStep()
        self.gid_output_process.PrintOutput()
        self.time +=1

    def FinalizeSolutionStep(self):
        super(KratosPotentialFlowSolver, self).FinalizeSolutionStep()
        self.gid_output_process.ExecuteFinalize()

    def _GetParallelType(self):
        return self.project_parameters["problem_data"]["parallel_type"].GetString()

    def _Name(self):
        return self.__class__.__name__